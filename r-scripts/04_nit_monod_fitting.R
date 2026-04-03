# Jason R Laurich
# April 3rd, 2026

# We're going to start fitting nitrogen Monod curves to our Chlamydomonas reinhardtii and microbial growth rate data.

# For now, we will not worry about microbial data from wells with algae.

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(rTPC)
library(nls.multstart)
library(Deriv)
library(car)
library(boot)
library(minpack.lm)

monod_fun <- function(S, r.max, K.s) {
  r.max * S / (K.s + S)
}

# Chlamy alone ------------------------------------------------------------

###### Load & examine the data ######

df.c <- read.csv("processed-data/02_chlamy_µs.csv")

head(df.c) # OK so we have unique.id, with temp, nit, salt, block, mic, plate, well, rep, and µ

df.c.n <- df.c %>% 
  filter(salt == 0,
         temp == 30)

df.c.n <- df.c.n %>%
  mutate(rep.id = paste0("b", block, ".m", mic, ".r", rep)) # So this column will capture unique id's for each nitrogen gradient. 

length(unique(df.c.n$rep.id)) # 214

df.c.n <- df.c.n %>% # Don't have data for the 4th replicates at all levels of nitrogen
  filter(rep <4) # leaves us with 204

###### Run the models ######

c.summ.df <- data.frame(    # We'll create a dataframe to store the data as we fit models.
  block = numeric(),        # Block
  mic = character(),        # Microbial inocula
  rep = numeric(),          # Replicate
  id = character(),         # rep.id
  
  K.s.mod = numeric(),      # Half saturation constant (model output)
  K.s.post = numeric(),     # Half saturation constant (posterior median)
  K.s.min = numeric(),      # Half saturation constant (lower HDPI)
  K.s.max = numeric(),      # Half saturation constant (upper HDPI)
  K.s.na = numeric(),       # % NA returns
  
  r.max.mod = numeric(),    # Maximum growth rate (model output)
  r.max.post = numeric(),   # Maximum growth rate (posterior median)
  r.max.min = numeric(),    # Maximum growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum growth rate (upper HDPI)
  r.max.na = numeric(),     # % NA returns
  
  R.post = numeric(),       # R* (m = 0.56) (posterior median)
  R.min = numeric(),        # R* (m = 0.56) (lower HDPI)
  R.max = numeric(),        # R* (m = 0.56) (upper HDPI)
  R.na = numeric(),         # % NA returns
  
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  block = numeric(),        # Block
  mic = character(),        # Microbial inocula
  rep = numeric(),          # Replicate
  id = character(),         # rep.id
  
  Parameter = character(),  # Model parameter (e.g. a, tmax, etc.)
  est = numeric(),          # Estimate
  se = numeric(),           # Standard error
  p = numeric(),            # p-value
  
  stringsAsFactors = FALSE            
)

n <-0 # progression tracker

for (i in unique(df.c.n$rep.id[df.c.n$rep.id >= 1])) { # for each replicate ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.c.n %>% 
    filter(rep.id == i)
  df.i <- droplevels(df.i)
  
  monod_nls <- nls_multstart(µ ~ monod_fun(S = nit, r.max = r.max, K.s = K.s),
                           data = df.i,
                           iter = c(10,10), 
                           start_lower = c(r.max = 0, K.s = 0),
                           start_upper = c(r.max = 5, K.s = 100),
                           lower = c(r.max = 0, K.s = 0),
                           upper = c(r.max = 10, K.s = 400),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  sum <- summary(monod_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  r.max <- df.nls[1,1] # Extract parameters
  k.s <- df.nls[2,1]
  
  monod.LM <- nlsLM(µ ~ monod_fun(S = nit, r.max = r.max, K.s = K.s),
                  data = df.i,
                  start = c(r.max = r.max, K.s = k.s),
                  lower = c(r.max = max(0, r.max - 1), K.s = max(0, k.s - 200)),
                  upper = c(r.max = r.max + 1, K.s = 1000),
                  control = nls.lm.control(maxiter=500)
  )
  
  sum <- summary(monod.LM)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  boot <- Boot(monod.LM, method = 'residual')
  
  post <- as.data.frame(boot$t) # Get the bootstrapped values
  post$R <- 0.56*post$K.s/(post$r.max - 0.56)
  
  c.summ.df <- rbind(c.summ.df, data.frame(                                     # Add summary data
    block = df.i$block[1],                                                      # Block
    mic = df.i$mic[1],                                                          # Microbial inocula
    rep = df.i$rep[1],                                                          # Replicate
    id = df.i$rep.id[1],                                                        # Replicate id
    
    K.s.mod = k.s,                                                              # Half saturation constant (model output)
    K.s.post = median(post$K.s),                                                # Half saturation constant (posterior median)
    K.s.min = hdi(post$K.s[!is.na(post$K.s)], credMass = 0.95)$CI_low,          # Half saturation constant (lower HDPI)
    K.s.max = hdi(post$K.s[!is.na(post$K.s)], credMass = 0.95)$CI_high,         # Half saturation constant (upper HDPI)
    K.s.na = mean(is.na(post$K.s)),                                             # % NA returns
    
    r.max.mod = r.max,                                                          # Maximum growth rate (model output)
    r.max.post = median(post$r.max),                                            # Maximum growth rate (posterior median)
    r.max.min = hdi(post$r.max[!is.na(post$r.max)], credMass = 0.95)$CI_low,    # Maximum growth rate (lower HDPI)
    r.max.max = hdi(post$r.max[!is.na(post$r.max)], credMass = 0.95)$CI_high,   # Maximum growth rate (upper HDPI)
    r.max.na = mean(is.na(post$r.max)),                                         # % NA returns
    
    R.post = median(post$R),                                                    # R* (m = 0.56) (posterior median)
    R.min = hdi(post$R[!is.na(post$R)], credMass = 0.95)$CI_low,                # R* (m = 0.56) (lower HDPI)
    R.max = hdi(post$R[!is.na(post$R)], credMass = 0.95)$CI_high,               # R* (m = 0.56) (upper HDPI)
    R.na = mean(is.na(post$R))                                                  # % NA returns
  ))
  
  for (j in 1:2){
    fit.df <- rbind(fit.df, data.frame(                                         # Model performance data
      block = df.i$block[1],                                                    # Block
      mic = df.i$mic[1],                                                        # Microbial inocula
      rep = df.i$rep[1],                                                        # Replicate
      id = df.i$rep.id[1],                                                      # Replicate id
      
      Parameter = df.nls$parameter[j],                                          # Model parameter (e.g. a, b, tmax etc.)
      est = df.nls$Estimate[j],                                                 # Estimate
      se = df.nls$`Std. Error`[j],                                              # Error
      p = df.nls$`Pr(>|t|)`[j],                                                 # p-values
      stringsAsFactors = FALSE            
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.c.n$rep.id))))
  
}

write.csv(c.summ.df, "processed-data/06a_chlamy_Monod_nit.csv") # 204 Monod curves!
write.csv(fit.df, "processed-data/06b_chlamy_Monod_nit_fits.csv")

# Jason R Laurich
# April 7th, 2026

# We're going to start fitting salt tolerance curves to our Chlamydomonas reinhardtii and microbial growth rate data.

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

salt_fun <- function(S, a, b, c) {
  a / (1 + exp(b * (S - c)))
}

# Chlamy alone ------------------------------------------------------------

###### Load & examine the data ######

df.c <- read.csv("processed-data/02_chlamy_µs.csv")

head(df.c) # OK so we have unique.id, with temp, nit, salt, block, mic, plate, well, rep, and µ

df.c.s <- df.c %>% 
  filter(nit == 1000,
         temp == 30)

df.c.s <- df.c.s %>%
  mutate(rep.id = paste0("b", block, ".m", mic, ".r", rep)) # So this column will capture unique id's for each salt gradient. 

length(unique(df.c.s$rep.id)) # 207

df.c.s <- df.c.s %>% # Don't have data for the 4th replicates at all levels of salt
  filter(rep <4) # leaves us with 204

###### Run the models ######

c.summ.df <- data.frame(    # We'll create a dataframe to store the data as we fit models.
  block = numeric(),        # Block
  mic = character(),        # Microbial inocula
  rep = numeric(),          # Replicate
  id = character(),         # rep.id
  
  r.max.mod = numeric(),    # Maximum population growth rate (alpha) (model output)
  r.max.post = numeric(),   # Maximum population growth rate (posterior median)
  r.max.min = numeric(),    # Maximum population growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum population growth rate (upper HDPI)
  r.max.na = numeric(),     # % NA returns
  
  c.mod = numeric(),        # salt concentration at which r is half of alpha (extracted from model)
  c.post = numeric(),       # salt concentration at which r is half of alpha (posterior median)
  c.post.min = numeric(),   # salt concentration at which r is half of alpha (lower HDPI)
  c.post.max = numeric(),   # salt concentration at which r is half of alpha (upper HDPI)
  c.post.na = numeric(),    # % NA returns
  
  b.mod = numeric(),        # decline rate (extracted from model)
  b.post = numeric(),       # decline rate (posterior median)
  b.post.min = numeric(),   # decline rate (lower HDPI)
  b.post.max = numeric(),   # decline rate (upper HDPI)
  b.post.na = numeric(),    # % NA returns
  
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

for (i in unique(df.c.s$rep.id[df.c.s$rep.id >= 1])) { # for each replicate ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.c.s %>% 
    filter(rep.id == i)
  df.i <- droplevels(df.i)
  
  mu.thresh <- max(df.i$µ[df.i$salt < 6], na.rm = TRUE) / 2
  
  df.i <- df.i %>%
    filter(!(salt >= 6 & µ > mu.thresh))
    
  
  salt_nls <- nls_multstart(µ ~ salt_fun(S = salt, a = a, b = b, c = c),
                             data = df.i,
                             iter = c(10, 10, 10), 
                             start_lower =  c(a = 0.5, b = 0,  c = 0),
                             start_upper = c(a = 2,   b = 10, c = 10),
                             lower = c(a = 0.5, b = 0,  c = 0),
                             upper = c(a = 2,   b = 10, c = 12),
                             supp_errors = 'Y',
                             convergence_count = FALSE
  )
  
  sum <- summary(salt_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  a.mod <- df.nls[1,1] # Extract parameters
  b.mod <- df.nls[2,1]
  c.mod <- df.nls[3,1]
  
  salt.LM <- nlsLM(µ ~ salt_fun(S = salt, a = a, b = b, c = c),
                    data = df.i,
                    start = c(a = a.mod, b = b.mod, c = c.mod),
                    lower = c(a = max(0, a.mod - 0.5), b = max(0, b.mod - 2), c = max (0, c.mod - 2)),
                    upper = c(a = a.mod + 0.5, b = b.mod + 2, c = c.mod +2),
                    control = nls.lm.control(maxiter=500)
  )
  
  sum <- summary(salt.LM)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  boot <- Boot(salt.LM, method = 'residual')
  
  post <- as.data.frame(boot$t) # Get the bootstrapped values
  
  c.summ.df <- rbind(c.summ.df, data.frame(                                     # Add summary data
    block = df.i$block[1],                                                      # Block
    mic = df.i$mic[1],                                                          # Microbial inocula
    rep = df.i$rep[1],                                                          # Replicate
    id = df.i$rep.id[1],                                                        # Replicate id
    
    r.max.mod = a.mod,                                                          # Maximum growth rate (model output)
    r.max.post = median(post$a),                                                # Maximum growth rate (posterior median)
    r.max.min = hdi(post$a[!is.na(post$a)], credMass = 0.95)$CI_low,            # Maximum growth rate (lower HDPI)
    r.max.max = hdi(post$a[!is.na(post$a)], credMass = 0.95)$CI_high,           # Maximum growth rate (upper HDPI)
    r.max.na = mean(is.na(post$a)),                                             # % NA returns
    
    c.mod = c.mod,                                                              # salt concentration at which r is half of alpha (extracted from model)
    c.post = median(post$c),                                                    # salt concentration at which r is half of alpha (posterior median)
    c.post.min = hdi(post$c[!is.na(post$c)], credMass = 0.95)$CI_low,           # salt concentration at which r is half of alpha (lower HDPI)
    c.post.max = hdi(post$c[!is.na(post$c)], credMass = 0.95)$CI_high,          # salt concentration at which r is half of alpha (upper HDPI)
    c.post.na = mean(is.na(post$c)),                                            # % NA returns
    
    b.mod = b.mod,                                                              # decline rate (extracted from model)
    b.post = median(post$b),                                                    # decline rate (posterior median)
    b.post.min = hdi(post$b[!is.na(post$b)], credMass = 0.95)$CI_low,           # decline rate (lower HDPI)
    b.post.max = hdi(post$b[!is.na(post$b)], credMass = 0.95)$CI_high,          # decline rate (upper HDPI)
    b.post.na = mean(is.na(post$b))                                             # % NA returns
  ))
  
  for (j in 1:3){
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
  
  print(paste("Done", n, "of ", length(unique(df.c.s$rep.id))))
  
}

write.csv(c.summ.df, "processed-data/07a_chlamy_salt_tol.csv") # 204 Monod curves!
write.csv(fit.df, "processed-data/07b_chlamy_salt_tol_fits.csv")

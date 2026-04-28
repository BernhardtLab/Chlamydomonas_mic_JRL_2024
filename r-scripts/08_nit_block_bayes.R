# Jason R Laurich
# April 16th, 2026

# We're going to play around with Monod nit curves, using R2jags. We'll fit models to each block, not replicate

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(R2jags)
library(mcmcplots)
library(bayestestR)

###### Load & examine the data ######

df.c <- read.csv("processed-data/02_chlamy_µs.csv")

head(df.c) # OK so we have unique.id, with temp, nit, salt, block, mic, plate, well, rep, and µ

df.c.n <- df.c %>% 
  filter(salt == 0,
         temp == 30)

df.c.n <- df.c.n %>%
  mutate(rep.id = paste0("b", block, ".m", mic)) # So this column will capture unique id's for each nitrogen gradient. 

length(unique(df.c.n$rep.id)) # 214

df.c.n <- df.c.n %>% # Don't have data for the 4th replicates at all levels of nitrogen
  filter(rep <4,
         block != 2) # leaves us with 204

ggplot(df.c.n, aes(x = nit, y = µ)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ rep.id, ncol = 6) +
  theme_classic()

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
  
  aff.mod = numeric(),      # affinity (1/Ks) (model output)
  aff.post = numeric(),     # affinity (1/Ks) (posterior median)
  aff.min = numeric(),      # affinity (1/Ks) (lower HDPI)
  aff.max = numeric(),      # affinity (1/Ks) (upper HDPI)
  aff.na = numeric(),       # % NA returns
  
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

# Set generous MCMC settings still for our models. 
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model. 

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 3), # Initial guess for r_max
    K_s = runif(1, 0.1, 10),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

S.pred <- seq(0, 1000, 0.5) # Nitrogen gradient we are interested in here (concentration)
N.S.pred <-length(S.pred)

n <-0 # progression tracker

for (i in unique(df.c.n$rep.id[df.c.n$rep.id >= 26])) { # for each replicate ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.c.n %>% 
    filter(rep.id == i)
  df.i <- droplevels(df.i)
  
  trait <- df.i$µ    # format the data for jags
  N.obs <- length(trait)
  
  nit <- df.i$nit
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = nit, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod.jag <- jags( # Run the nitrogen Monod function. 
    data = jag.data,
    inits = inits.monod,
    parameters.to.save = parameters.monod,
    model.file = "monod.nit.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  save(monod.jag, file = paste0("R2jags-models/rep_", i, "_nit_monod.RData")) # save the monod model
  # This folder is listed in gitignore, because the objects are too big to load
  
  post <- as.data.frame(monod.jag$BUGSoutput$sims.matrix) # The posteriors
  
  post <- post %>% 
    select(K_s, r_max)
  
  post$aff <- 1/post$K_s
  
  c.summ.df <- rbind(c.summ.df, data.frame(                                     # Add summary data
    block = df.i$block[1],                                                      # Block
    mic = df.i$mic[1],                                                          # Microbial inocula
    rep = df.i$rep[1],                                                          # Replicate
    id = df.i$rep.id[1],                                                        # Replicate id
    
    K.s.mod = monod.jag$BUGSoutput$summary[1,1],                                # Half saturation constant (model output)
    K.s.post = median(post$K_s, na.rm = T),                                     # Half saturation constant (posterior median)
    K.s.min = hdi(post$K_s, ci = 0.95)$CI_low,                                  # Half saturation constant (lower HDPI)
    K.s.max = hdi(post$K_s, ci = 0.95)$CI_high,                                 # Half saturation constant (upper HDPI)
    K.s.na = mean(is.na(post$K_s)),                                             # % NA returns
    
    r.max.mod = monod.jag$BUGSoutput$summary[3,1],                              # Maximum growth rate (model output)
    r.max.post = median(post$r_max, na.rm = T),                                 # Maximum growth rate (posterior median)
    r.max.min = hdi(post$r_max, ci = 0.95)$CI_low,                              # Maximum growth rate (lower HDPI)
    r.max.max = hdi(post$r_max, ci = 0.95)$CI_high,                             # Maximum growth rate (upper HDPI)
    r.max.na = mean(is.na(post$r_max)),                                         # % NA returns
    
    aff.mod = 1/monod.jag$BUGSoutput$summary[1,1],                              # affinity (1/Ks) (model output)
    aff.post = median(post$aff, na.rm = T),                                     # affinity (1/Ks) (posterior median)
    aff.min = hdi(post$aff, ci = 0.95)$CI_low,                                  # affinity (1/Ks) (lower HDPI)
    aff.max = hdi(post$aff, ci = 0.95)$CI_high,                                 # affinity (1/Ks) (upper HDPI)
    aff.na = mean(is.na(post$aff))                                              # % NA returns
    
  ))
  
  nit_sum <- monod.jag$BUGSoutput$summary[c(1:3),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:3){
    fit.df <- rbind(fit.df, data.frame(                                         # Model performance data
      block = df.i$block[1],                                                    # Block
      mic = df.i$mic[1],                                                        # Microbial inocula
      rep = df.i$rep[1],                                                        # Replicate
      id = df.i$rep.id[1],                                                      # Replicate id
      
      Parameter = rownames(nit_sum)[j],          # Model parameter (e.g. K_s, r_max, etc.)
      mean = nit_sum[j,1],                       # Posterior mean
      Rhat = nit_sum[j,8],                       # Rhat values
      n.eff = nit_sum[j,9]                       # Sample size estimates (should be ~6000)
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.c.n$rep.id))))
  
}

write.csv(c.summ.df, "processed-data/11a_chlamy_Monod_nit_bayes.csv") 
write.csv(fit.df, "processed-data/11b_chlamy_Monod_nit_fits_bayes.csv")

# Load and compile if needed ----------------------------------------------

n <-0 # progression tracker

for (i in unique(df.c.n$rep.id[df.c.n$rep.id >= 1])) { # for each replicate ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.c.n %>% 
    filter(rep.id == i)
  df.i <- droplevels(df.i)
  
  load(paste0("R2jags-models/rep_", i, "_nit_monod.RData")) # load the lactin2 model
  
  post <- as.data.frame(monod.jag$BUGSoutput$sims.matrix) # The posteriors
  
  post <- post %>% 
    select(K_s, r_max)
  
  post$aff <- 1/post$K_s
  
  c.summ.df <- rbind(c.summ.df, data.frame(                                     # Add summary data
    block = df.i$block[1],                                                      # Block
    mic = df.i$mic[1],                                                          # Microbial inocula
    rep = df.i$rep[1],                                                          # Replicate
    id = df.i$rep.id[1],                                                        # Replicate id
    
    K.s.mod = monod.jag$BUGSoutput$summary[1,1],                                # Half saturation constant (model output)
    K.s.post = median(post$K_s, na.rm = T),                                     # Half saturation constant (posterior median)
    K.s.min = hdi(post$K_s, ci = 0.95)$CI_low,                                  # Half saturation constant (lower HDPI)
    K.s.max = hdi(post$K_s, ci = 0.95)$CI_high,                                 # Half saturation constant (upper HDPI)
    K.s.na = mean(is.na(post$K_s)),                                             # % NA returns
    
    r.max.mod = monod.jag$BUGSoutput$summary[3,1],                              # Maximum growth rate (model output)
    r.max.post = median(post$r_max, na.rm = T),                                 # Maximum growth rate (posterior median)
    r.max.min = hdi(post$r_max, ci = 0.95)$CI_low,                              # Maximum growth rate (lower HDPI)
    r.max.max = hdi(post$r_max, ci = 0.95)$CI_high,                             # Maximum growth rate (upper HDPI)
    r.max.na = mean(is.na(post$r_max)),                                         # % NA returns
    
    aff.mod = 1/monod.jag$BUGSoutput$summary[1,1],                              # affinity (1/Ks) (model output)
    aff.post = median(post$aff, na.rm = T),                                     # affinity (1/Ks) (posterior median)
    aff.min = hdi(post$aff, ci = 0.95)$CI_low,                                  # affinity (1/Ks) (lower HDPI)
    aff.max = hdi(post$aff, ci = 0.95)$CI_high,                                 # affinity (1/Ks) (upper HDPI)
    aff.na = mean(is.na(post$aff))                                              # % NA returns
    
  ))
  
  nit_sum <- monod.jag$BUGSoutput$summary[c(1:3),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:3){
    fit.df <- rbind(fit.df, data.frame(                                         # Model performance data
      block = df.i$block[1],                                                    # Block
      mic = df.i$mic[1],                                                        # Microbial inocula
      rep = df.i$rep[1],                                                        # Replicate
      id = df.i$rep.id[1],                                                      # Replicate id
      
      Parameter = rownames(nit_sum)[j],          # Model parameter (e.g. K_s, r_max, etc.)
      mean = nit_sum[j,1],                       # Posterior mean
      Rhat = nit_sum[j,8],                       # Rhat values
      n.eff = nit_sum[j,9]                       # Sample size estimates (should be ~6000)
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.c.n$rep.id))))
  
}

write.csv(c.summ.df, "processed-data/100a_chlamy_Monod_nit_bayes.csv") 
write.csv(fit.df, "processed-data/100b_chlamy_Monod_nit_fits_bayes.csv")

# Jason R Laurich
# April 17th, 2026

# We're going to play around with salt tolerance curves, using R2jags. We'll fit models to each block, not replicate

# Packages & functions ----------------------------------------------------

library(tidyverse)

library(R2jags)
library(mcmcplots)
library(bayestestR)

###### Load & examine the data ######

df.c <- read.csv("processed-data/02_chlamy_µs.csv")

head(df.c) # OK so we have unique.id, with temp, nit, salt, block, mic, plate, well, rep, and µ

df.c.s <- df.c %>% 
  filter(nit == 1000,
         temp == 30)

df.c.s <- df.c.s %>%
  mutate(rep.id = paste0("b", block, ".m", mic)) # So this column will capture unique id's for each salt gradient. 

length(unique(df.c.s$rep.id)) # 68

df.c.s <- df.c.s %>% # Don't have data for the 4th replicates at all levels of salt
  filter(rep <4,
         block !=2) # leaves us with 51

ggplot(df.c.s, aes(x = salt, y = µ)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ rep.id, ncol = 6) +
  theme_classic() # Two weird data points we'll trim. salt >10, mu>1.5 

df.c.s <- df.c.s %>%
  filter(!(salt >= 8 & µ > 1.5))

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

inits.salt <- function() { # Smaller a (prior 0.5-> 2)
  list(
    a = runif(1, 0.7, 1.8),  # Smaller window for a
    b = runif(1, 0.1, 5),  # Initial guess for b
    c = runif(1, 0.1, 9),  # Initial guess for c
    sigma = runif(1, 0.1, 2)  # Initial guess for error
  )
}

parameters.salt <- c("a", "b", "c", "sigma", "r_pred_new") # Repeats

S.pred <- seq(0, 10, 0.005) # Repeats
N.S.pred <-length(S.pred)

# Repeated here
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model.

n <-0 # progression tracker

for (i in unique(df.c.s$rep.id[df.c.s$rep.id >= 1])) { # for each replicate ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.c.s %>% 
    filter(rep.id == i)
  df.i <- droplevels(df.i)

  trait <- df.i$µ    # format the data for jags
  N.obs <- length(trait)
  
  salt <- df.i$salt
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)
  
  salt.jag <- jags( # Run the salt logistic growth curve function. 
    data = jag.data,
    inits = inits.salt,
    parameters.to.save = parameters.salt,
    model.file = "salt.tolerance.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  save(salt.jag, file = paste0("R2jags-models/pop_", i, "_salt_tol.RData")) # save the monod model
  # This folder is listed in gitignore, because the objects are too big to load
  
  post <- as.data.frame(salt.jag$BUGSoutput$sims.matrix) # The posteriors
  
  post <- post %>% 
    select(a, b, c)
  
  c.summ.df <- rbind(c.summ.df, data.frame(                                     # Add summary data
    block = df.i$block[1],                                                      # Block
    mic = df.i$mic[1],                                                          # Microbial inocula
    rep = df.i$rep[1],                                                          # Replicate
    id = df.i$rep.id[1],                                                        # Replicate id
    
    r.max.mod = salt.jag$BUGSoutput$summary[1,1],                               # Maximum growth rate (model output)
    r.max.post = median(post$a, na.rm = T),                                     # Maximum growth rate (posterior median)
    r.max.min = hdi(post$a, ci = 0.95)$CI_low,                                  # Maximum growth rate (lower HDPI)
    r.max.max = hdi(post$a, ci = 0.95)$CI_high,                                 # Maximum growth rate (upper HDPI)
    r.max.na = mean(is.na(post$a)),                                             # % NA returns
    
    c.mod = salt.jag$BUGSoutput$summary[3,1],                                   # salt concentration at which r is half of alpha (extracted from model)
    c.post = median(post$c, na.rm = T),                                         # salt concentration at which r is half of alpha (posterior median)
    c.post.min = hdi(post$c, ci = 0.95)$CI_low,                                 # salt concentration at which r is half of alpha (lower HDPI)
    c.post.max = hdi(post$c, ci = 0.95)$CI_high,                                # salt concentration at which r is half of alpha (upper HDPI)
    c.post.na = mean(is.na(post$c)),                                            # % NA returns
    
    b.mod = salt.jag$BUGSoutput$summary[2,1],                                   # decline rate (extracted from model)
    b.post = median(post$b, na.rm = T),                                         # decline rate (posterior median)
    b.post.min = hdi(post$b, ci = 0.95)$CI_low,                                 # decline rate (lower HDPI)
    b.post.max = hdi(post$b, ci = 0.95)$CI_low,                                 # decline rate (upper HDPI)
    b.post.na = mean(is.na(post$b))                                             # % NA returns
  ))
  
  salt_sum <- salt.jag$BUGSoutput$summary[1:4,] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:3){
    fit.df <- rbind(fit.df, data.frame(                                         # Model performance data
      block = df.i$block[1],                                                    # Block
      mic = df.i$mic[1],                                                        # Microbial inocula
      rep = df.i$rep[1],                                                        # Replicate
      id = df.i$rep.id[1],                                                      # Replicate id
      
      Parameter = rownames(salt_sum)[j],         # Model parameter (e.g. K_s, r_max, etc.)
      mean = salt_sum[j,1],                      # Posterior mean
      Rhat = salt_sum[j,8],                      # Rhat values
      n.eff = salt_sum[j,9]                      # Sample size estimates (should be ~6000)      
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.c.s$rep.id))))
  
}

write.csv(c.summ.df, "processed-data/12a_chlamy_salt_tol_bayes.csv") # 204 salt tolerance curves!
write.csv(fit.df, "processed-data/12b_chlamy_salt_tol_fits_bayes.csv")

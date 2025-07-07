# Jason R Laurich

# July 7th, 2025

# Going to fit TPCs, Monod curves and salt tolerance curves to data from block 1 (carbon added) and 2 (no carbon)

# Load packages -----------------------------------------------------------

library(tidyr)
library(cowplot)
library(tidyverse)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(rTPC)
library(nls.multstart)

# Load and examine data ---------------------------------------------------

### Temperature data ### 

df.t1 <- read.csv("processed-data/2a_blk1_temp_mus.csv") # block 1, temp
df.t1$block <- 1

df.t2 <- read.csv("processed-data/4a_blk2_temp_mus.csv") # block 2, temp
df.t2$block <- 2

df.t <- rbind(df.t1, df.t2)

df.t$mic <- as.factor(df.t$mic)

df.t <- df.t %>% 
  filter(mic != 'BLANK') %>%
  droplevels() %>%
  mutate(mic = fct_relevel(mic, 
                               "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                               "9", "10", "11", "12", "13", "14", "15", "all")) %>% 
  mutate(blk.mic = paste0(block, mic))

### Nitrogen data ### 

df.n1 <- read.csv("processed-data/2b_blk1_nit_mus.csv") # block 1, nit
df.n1$block <- 1

df.n2 <- read.csv("processed-data/4b_blk2_nit_mus.csv") # block 2, nit
df.n2$block <- 2

df.n <- rbind(df.n1, df.n2)

df.n$mic <- as.factor(df.n$mic)

df.n <- df.n %>% 
  filter(mic != 'BLANK') %>%
  droplevels() %>%
  mutate(mic = fct_relevel(mic, 
                           "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                           "9", "10", "11", "12", "13", "14", "15", "all")) %>% 
  mutate(blk.mic = paste0(block, mic))

### Salt data ###

df.s1 <- read.csv("processed-data/2c_blk1_salt_mus.csv") # block 1, salt
df.s1$block <- 1

df.s2 <- read.csv("processed-data/4c_blk2_salt_mus.csv") # block 2, salt
df.s2$block <- 2

df.s <- rbind(df.s1, df.s2)

df.s$mic <- as.factor(df.s$mic)

df.s <- df.s %>% 
  filter(mic != 'BLANK') %>%
  droplevels() %>%
  mutate(mic = fct_relevel(mic, 
                           "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                           "9", "10", "11", "12", "13", "14", "15", "all")) %>% 
  mutate(blk.mic = paste0(block, mic))

# Fit TPCs ----------------------------------------------------------------

blk1.2.tpc.summ.df <- data.frame(  # We'll create a dataframe to store the data as we fit models.
  Block = numeric(),               # Block
  Mic = character(),               # Mic
  T.min.raw = numeric(),           # Minimum T (Jags raw)
  T.max.raw = numeric(),           # Maximum T (Jags raw)
  T.opt.raw = numeric(),           # Optimal T (Jags raw)
  r.max.raw = numeric(),           # Maximum growth rate (Jags raw)
  T.br.raw = numeric(),            # T breadth (Jags raw)
  stringsAsFactors = FALSE         # Avoid factor conversion
)

# Let's do larger models for the final things (10 times larger)
ni.fit <- 220000    # iterations / chain
nb.fit <- 20000     # burn in periods for each chain
nt.fit <- 200       # thinning interval : (200,000 - 20,000) / 200 = 1000 posterior estimates / chain
nc.fit <- 2         # number of chains, total of 2,000 estimates for each model. 

parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred") # repeated here

inits.lactin.cust<- function() { # Pulling initial values centres from the start_vals function in rTPC
  list(
    cf.a = rnorm(1, mean = start.vals.lac[1], sd = 0.05),
    cf.tmax = rnorm(1, mean = start.vals.lac[3], sd = 1),
    cf.delta_t = rnorm(1, mean = start.vals.lac[4], sd = 1),
    cf.b = rnorm(1, mean = start.vals.lac[2], sd = 0.05),
    cf.sigma = runif(1, 0.1, 2)
  )
}

for (i in unique(df.t$blk.mic)[32:length(unique(df.t$blk.mic))]){ # For each unique block x mic combo
  
  df.i <- df.t %>%
    filter(blk.mic == i, !is.na(r.exp)) %>%
    arrange(temp)
  
  trait <- df.i$r.exp    # format the data for jags
  N.obs <- length(trait)
  
  Temp.xs <- seq(min(df.i$temp) - 3, max(df.i$temp) + 3, 0.1) # Temperature gradient we're interested in - upped the granularity here
  N.Temp.xs <-length(Temp.xs) # We'll reset this internally since the gradient varies substantially
  
  temp <- df.i$temp
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  start.vals.lac <- get_start_vals(df.i$temp, df.i$r.exp, model_name = 'lactin2_1995')
  
  lac_jag <- jags(
    data = jag.data, 
    inits = inits.lactin.cust, 
    parameters.to.save = parameters.lactin2, 
    model.file = "lactin2_flex.txt",
    n.thin = nt.fit, 
    n.chains = nc.fit, 
    n.burnin = nb.fit, 
    n.iter = ni.fit, 
    DIC = TRUE, 
    working.directory = getwd()
  ) # ~ 10 min to run?
  
  print(paste("Done", i))
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(min(df.i$temp) - 3, max(df.i$temp) + 3, 0.1)
  
  blk1.2.tpc.summ.df <- rbind(blk1.2.tpc.summ.df, data.frame(   # Add summary data
    Block = df.i$block[1],                                      # Block
    Mic = df.i$mic[1],                                          # Mic
    T.min.raw = df.jags$temp[min(which(df.jags$mean > 0))],     # Minimum T
    T.max.raw = df.jags$temp[max(which(df.jags$mean > 0))],     # Maximum T
    T.br.raw = df.jags$temp[max(which(df.jags$mean > 0.5))] - 
      df.jags$temp[min(which(df.jags$mean > 0.5))],             # T breadth
    T.opt.raw = df.jags$temp[which.max(df.jags$mean)],          # Optimal T
    r.max.raw = max(df.jags$mean)                               # Maximum growth rate   
  ))
}

write.csv(blk1.2.tpc.summ.df, "processed-data/5a_blocks1.2_TPC_stats.csv") # save the file.

# Nitrogen Monod curves ---------------------------------------------------

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

blk1.2.n.monod.summ.df <- data.frame(  # We'll create a dataframe to store the data as we fit models.
  Block = numeric(),                 # Block
  Mic = character(),                 # Mic
  K.s = numeric(),                   # Half-saturation constant
  r.max = numeric(),                 # Maximum population growth rate
  R.jag = numeric(),                 # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),                 # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
  stringsAsFactors = FALSE           # Avoid factor conversion
)

for (i in unique(df.n$blk.mic)[1:length(unique(df.n$blk.mic))]){ # For each unique block x mic combo
  
  df.i <- df.n %>%
    filter(blk.mic == i, !is.na(r.exp)) %>%
    arrange(nit)
  
  trait <- df.i$r.exp    # format the data for jags
  N.obs <- length(trait)
  
  S.pred <- seq(0, 1000, 1) # Nitrogen gradient we're interested in - upped the granularity here
  N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
  nit <- df.i$nit
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = nit, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the light Monod function. 
    data = jag.data,
    inits = inits.monod,
    parameters.to.save = parameters.monod,
    model.file = "monod.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  print(paste("Done", i))
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3, 1005),]   # generate the sequence of r.pred values
  df.jags$nit <- seq(0, 1000, 1)
  
  blk1.2.n.monod.summ.df <- rbind(blk1.2.n.monod.summ.df, data.frame(                         # Add summary data
    Block = df.i$block[1],                                                                    # Block
    Mic = df.i$mic[1],                                                                        # Mic
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                  # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                # Maximum population growth rate
    R.jag = df.jags$nit[which(df.jags$mean > 0.5)[1]],                                        # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.5*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.5)   # Minimum resource requirement for positive growth (from math)
  ))
}

write.csv(blk1.2.n.monod.summ.df, "processed-data/5b_blocks1.2_N_monod_stats") # save the file.

# Salt tolerance reversed logisitic growth curves -------------------------

inits.salt <- function() {
  list(
    a = runif(1, 0.1, 5),  # Initial guess for a
    b = runif(1, 0.1, 5),  # Initial guess for b
    c = runif(1, 0.1, max(df.i$salt)),  # Initial guess for c
    sigma = runif(1, 0.1, 2)  # Initial guess for error
  )
}

parameters.salt <- c("a", "b", "c", "sigma", "r_pred_new") # Save these

blk1.2.salt.summ.df <- data.frame(  # We'll create a dataframe to store the data as we fit models.
  Block = numeric(),                # Block
  Mic = character(),                # Mic
  r.max = numeric(),                # Maximum population growth rate (alpha)
  c.mod = numeric(),                # salt concentration at which r is half of alpha (extracted from model)
  c.pred = numeric(),               # salt concentration at which r is half of alpha (extracted from predicted values)
  stringsAsFactors = FALSE          # Avoid factor conversion
)


for (i in unique(df.s$blk.mic)[1:length(unique(df.s$blk.mic))]){ # For each unique block x mic combo
  
  df.i <- df.s %>%
    filter(blk.mic == i, !is.na(r.exp)) %>%
    arrange(salt)
  
  trait <- df.i$r.exp    # format the data for jags
  N.obs <- length(trait)
  
  S.pred <- seq(0, 12, 0.025) # Salt gradient we are interested in we'll keep N.S.pred more or less consistent across abiotic gradients for now
  N.S.pred <-length(S.pred)
  
  salt <- df.i$salt
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)
  
  salt_jag <- jags( # Run the salt logistic growth curve function. 
    data = jag.data,
    inits = inits.salt,
    parameters.to.save = parameters.salt,
    model.file = "salt.tol.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  print(paste("Done", i))
  
  df.jags <- data.frame(salt_jag$BUGSoutput$summary)[-c(1:4, 486),]   # generate the sequence of r.pred values
  df.jags$salt <- seq(0, 12, 0.025)

  blk1.2.salt.summ.df <- rbind(blk1.2.salt.summ.df, data.frame(                                         # Add summary data
    Block = df.i$block[1],                                                                                    # Block
    Mic = df.i$mic[1],                                                                                        # Mic
    r.max = salt_jag$BUGSoutput$summary[1,1],                                                                # Maximum population growth rate
    c.mod = salt_jag$BUGSoutput$summary[3,1],                                                                # salt concentration at which r is half of alpha (extracted from model)
    c.pred = df.jags$salt[which.min(abs(df.jags$mean - (monod_jag$BUGSoutput$summary[1,1] / 2)))]   # salt concentration at which r is half of alpha (extracted from predicted values)                                                   
  ))

}

write.csv(blk1.2.salt.summ.df, "processed-data/5c_blocks1.2_salt_stats") # save the file.



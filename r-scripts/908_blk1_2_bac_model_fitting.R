# Jason R Laurich

# July 14th, 2025

# Going to fit TPCs, Monod curves and salt tolerance curves to data from block 1 (carbon added) and 2 (no carbon) for microbes

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

df.t <- read.csv("processed-data/6a_blk1.2_mic_temp_mus.csv") 
df.t$mic <- as.factor(df.t$mic)

df.t <- df.t %>% 
  filter(mic != 'BLANK') %>%
  droplevels() %>%
  mutate(mic = fct_relevel(mic, 
                           "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                           "9", "10", "11", "12", "13", "14", "15", "all")) %>% 
  mutate(blk.mic = paste(block, mic, Chlamy, sep = "."))

### Nitrogen data ### 

df.n <- read.csv("processed-data/6b_blk1.2_mic_nit_mus.csv") 
df.n$mic <- as.factor(df.n$mic)

df.n <- df.n %>% 
  filter(mic != 'BLANK') %>%
  droplevels() %>%
  mutate(mic = fct_relevel(mic, 
                           "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                           "9", "10", "11", "12", "13", "14", "15", "all")) %>% 
  mutate(blk.mic = paste(block, mic, Chlamy, sep = "."))

### Salt data ###

df.s <- read.csv("processed-data/6c_blk1.2_mic_salt_mus.csv") 
df.s$mic <- as.factor(df.s$mic)

df.s <- df.s %>% 
  filter(mic != 'BLANK') %>%
  droplevels() %>%
  mutate(mic = fct_relevel(mic, 
                           "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                           "9", "10", "11", "12", "13", "14", "15", "all")) %>% 
  mutate(blk.mic = paste(block, mic, Chlamy, sep = "."))

# Fit TPCs ----------------------------------------------------------------

blk1.2.tpc.summ.df <- data.frame(  # We'll create a dataframe to store the data as we fit models.
  Block = numeric(),               # Block
  Mic = character(),               # Mic
  Chlamy = character(),            # Chlamy present?
  T.min.raw = numeric(),           # Minimum T (Jags raw)
  T.max.raw = numeric(),           # Maximum T (Jags raw)
  T.opt.raw = numeric(),           # Optimal T (Jags raw)
  r.max.raw = numeric(),           # Maximum growth rate (Jags raw)
  T.br.raw = numeric(),            # T breadth (Jags raw)
  stringsAsFactors = FALSE         # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Block = numeric(),        # Block
  Mic = character(),        # Mic
  Chlamy = character(),            # Chlamy present?
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~6000)
  stringsAsFactors = FALSE            
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

#### Start here ####

for (i in unique(df.t$blk.mic)[2:length(unique(df.t$blk.mic))]){ # For each unique block x mic combo x chlamy
  
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
  
  fit_success <- FALSE
  attempt <- 1
  
  while (!fit_success && attempt <= 5) {
    try_result <- tryCatch({
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
      )
      fit_success <- TRUE
      lac_jag
    }, error = function(e) {
      message(paste("Attempt", attempt, "failed with error:", conditionMessage(e)))
      attempt <<- attempt + 1
      NULL
    })
  }
  
  if (!fit_success) {
    warning(paste("Model fitting failed after 5 attempts for", i))
    next  # Skip to the next blk.mic if all attempts fail
  }
  
  print(paste("Done", i))
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(min(df.i$temp) - 3, max(df.i$temp) + 3, 0.1)
  
  blk1.2.tpc.summ.df <- rbind(blk1.2.tpc.summ.df, data.frame(   # Add summary data
    Block = df.i$block[1],                                      # Block
    Mic = df.i$mic[1],                                          # Mic
    Chlamy = df.i$Chlamy[1],                                    # Chlamy present?
    T.min.raw = df.jags$temp[min(which(df.jags$mean > 0))],     # Minimum T
    T.max.raw = df.jags$temp[max(which(df.jags$mean > 0))],     # Maximum T
    T.br.raw = df.jags$temp[max(which(df.jags$mean > 0.5))] - 
      df.jags$temp[min(which(df.jags$mean > 0.5))],             # T breadth
    T.opt.raw = df.jags$temp[which.max(df.jags$mean)],          # Optimal T
    r.max.raw = max(df.jags$mean)                               # Maximum growth rate   
  ))
  
  for (j in 1:6){
    fit.df <- rbind(fit.df, data.frame(                         # Model performance data
      Block = df.i$block[1],                                    # Block
      Mic = df.i$mic[1],                                        # Mic
      Chlamy = df.i$Chlamy[1],                                  # Chlamy present?
      Parameter = rownames(lac_jag$BUGSoutput$summary)[j],      # Model parameter (e.g. cf.a, cf.tmax, etc.)
      mean = lac_jag$BUGSoutput$summary[j,1],                   # Posterior mean
      Rhat = lac_jag$BUGSoutput$summary[j,8],                   # Rhat values
      n.eff = lac_jag$BUGSoutput$summary[j,9],                  # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
  
}

write.csv(blk1.2.tpc.summ.df, "processed-data/7a_blocks1.2_mic_TPC_stats.csv") # save the file.
write.csv(fit.df, "processed-data/7b_block1.2_mic_TPC_fits.csv")


# Nitrogen Monods ---------------------------------------------------------

inits.monod <- function() {
  list(
    r_max = runif(1, 0.1, 5),
    K_s = runif(1, 0.1, 5),
    sigma = runif(1, 0.1, 1)
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new")

blk1.2.n.monod.summ.df <- data.frame(
  Block = numeric(),
  Mic = character(),
  Chlamy = character(),
  K.s = numeric(),
  r.max = numeric(),
  R.jag = numeric(),
  R.mth = numeric(),
  stringsAsFactors = FALSE
)

for (i in unique(df.n$blk.mic)) {
  
  df.i <- df.n %>%
    filter(blk.mic == i, !is.na(r.exp)) %>%
    arrange(nit)
  
  # Prep data for JAGS
  trait <- df.i$r.exp
  N.obs <- length(trait)
  S.pred <- seq(0, 1000, 1)
  N.S.pred <- length(S.pred)
  nit <- df.i$nit
  jag.data <- list(trait = trait, N.obs = N.obs, S = nit, S.pred = S.pred, N.S.pred = N.S.pred)
  
  # Retry up to 5 times
  success <- FALSE
  for (attempt in 1:5) {
    tryCatch({
      monod_jag <- jags(
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
      success <- TRUE
      break
    }, error = function(e) {
      message(paste("Attempt", attempt, "failed for", i, ":", e$message))
    })
  }
  
  if (!success) {
    message(paste("All attempts failed for", i))
    next
  }
  
  # Process successful output
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3, 1005), ]
  df.jags$nit <- S.pred
  
  blk1.2.n.monod.summ.df <- rbind(blk1.2.n.monod.summ.df, data.frame(
    Block = df.i$block[1],
    Mic = df.i$mic[1],
    Chlamy = df.i$Chlamy[1],
    K.s = monod_jag$BUGSoutput$summary["K_s", "mean"],
    r.max = monod_jag$BUGSoutput$summary["r_max", "mean"],
    R.jag = df.jags$nit[which(df.jags$mean > 0.5)[1]],
    R.mth = 0.5 * monod_jag$BUGSoutput$summary["K_s", "mean"] /
      (monod_jag$BUGSoutput$summary["r_max", "mean"] - 0.5)
  ))
  
  message(paste("Done", i))
}

write.csv(blk1.2.n.monod.summ.df, "processed-data/7c_blocks1.2_N_monod_stats.csv", row.names = FALSE)


# Salt --------------------------------------------------------------------

inits.salt <- function() {
  list(
    a = runif(1, 0.1, 5),  # Initial guess for a
    b = runif(1, 0.1, 5),  # Initial guess for b
    c = runif(1, 0.1, max(df.i$salt)),  # Initial guess for c
    sigma = runif(1, 0.1, 2)  # Initial guess for error
  )
}

parameters.salt <- c("a", "b", "c", "sigma", "r_pred_new") # Save these

blk1.2.s.salt.summ.df <- data.frame(
  Block = numeric(),
  Mic = character(),
  Chlamy = character(),
  r.max = numeric(),                # Maximum population growth rate (alpha)
  c.mod = numeric(),                # salt concentration at which r is half of alpha (extracted from model)
  c.pred = numeric(),               # salt concentration at which r is half of alpha (extracted from predicted values)
  stringsAsFactors = FALSE          # Avoid factor conversion
)

fit.df <- data.frame(
  Block = numeric(),
  Mic = character(),
  Chlamy = character(),
  Parameter = character(),  # Model parameter (e.g. K_s, r_max, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~6000)
  stringsAsFactors = FALSE      
)

for (i in unique(df.s$blk.mic)) {
  
  df.i <- df.s %>%
    filter(blk.mic == i, !is.na(r.exp)) %>%
    arrange(salt)
  
  trait <- df.i$r.exp
  N.obs <- length(trait)
  
  S.pred <- seq(0, 12, 0.025) # Salt gradient we are interested in we'll keep N.S.pred more or less consistent across abiotic gradients for now
  N.S.pred <-length(S.pred)
  
  salt <- df.i$salt
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)
  
  success <- FALSE
  for (attempt in 1:5) {
    tryCatch({
      salt_jag <- jags(
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
      success <- TRUE
      break
    }, error = function(e) {
      message(paste("Attempt", attempt, "failed for", i, ":", e$message))
    })
  }
  
  if (!success) {
    message(paste("All attempts failed for", i))
    next
  }
  
  df.jags <- data.frame(salt_jag$BUGSoutput$summary)[-c(1:4, 486),]   # generate the sequence of r.pred values
  df.jags$salt <- seq(0, 12, 0.025)
  
  blk1.2.s.salt.summ.df <- rbind(blk1.2.s.salt.summ.df, data.frame(
    Block = df.i$block[1],
    Mic = df.i$mic[1],
    Chlamy = df.i$Chlamy[1],
    r.max = salt_jag$BUGSoutput$summary[1,1],                                                       # Maximum population growth rate
    c.mod = salt_jag$BUGSoutput$summary[3,1],                                                       # salt concentration at which r is half of alpha (extracted from model)
    c.pred = df.jags$salt[which.min(abs(df.jags$mean - (salt_jag$BUGSoutput$summary[1,1] / 2)))]   # salt concentration at which r is half of alpha (extracted from predicted values)                                                   
  ))
  
  salt_sum <- salt_jag$BUGSoutput$summary[c(1:4, 486),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(        # Model performance data
      Block = df.i$block[1],                   # Block
      Mic = df.i$mic[1],                       # Mic 
      Chlamy = df.i$Chlamy[1],
      Parameter = rownames(salt_sum)[j],       # Model parameter (e.g. K_s, r_max, etc.)
      mean = salt_sum[j,1],                    # Posterior mean
      Rhat = salt_sum[j,8],                    # Rhat values
      n.eff = salt_sum[j,9],                   # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  
  }
  
  message(paste("Done", i))
}

write.csv(blk1.2.s.salt.summ.df, "processed-data/7d_blocks1.2_salt_stats.csv", row.names = FALSE)
write.csv(fit.df, "processed-data/7e_blocks1.2_salt_fits.csv", row.names = FALSE)

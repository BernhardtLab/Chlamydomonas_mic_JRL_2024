# Jason R Laurich
# April 14th, 2026

# We're going to play around with TPC curves, using R2jags. We'll fit models to each block, not replicate

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(rTPC)
library(nls.multstart)
library(Deriv)
library(MuMIn)

library(R2jags)
library(mcmcplots)
library(bayestestR)

lactin2 <- function(temp, a, b, tmax, d.t) { 
  exp(a * temp) - exp(a * tmax - (tmax - temp) / d.t) + b
} # Define the Lactin II function

lactin2_deriv <- function(temp, a, b, tmax, d.t) {
  rho <- a
  T_max <- tmax
  delta_T <- d.t
  
  term1 <- rho * exp(rho * temp)
  term2 <- (1 / delta_T) * exp(rho * T_max - (T_max - temp) / delta_T)
  
  return(term1 - term2)
} # Derivative of the Lactin II function

lactin2_halfmax <- function(temp, a, b, tmax, d.t, r_half) {
  exp(a * temp) - exp(a * tmax - (tmax - temp) / d.t) + b - r_half
} # Calculate thermal breadth at 1/2 µ_max. 

###### Load & examine the data ######

df.c <- read.csv("processed-data/02_chlamy_µs.csv")

head(df.c) # OK so we have unique.id, with temp, nit, salt, block, mic, plate, well, rep, and µ

df.c.t <- df.c %>% 
  filter(salt == 0,
         nit == 1000)

df.c.t <- df.c.t %>%
  mutate(rep.id = paste0("b", block, ".m", mic)) # So this column will capture unique id's for each termperature gradient. 

length(unique(df.c.t$rep.id)) # 68

df.c.t <- df.c.t %>% # Don't have data for the 4th replicates at all levels of temperature
  filter(rep <4,
         block != 2)

ggplot(df.c.t, aes(x = temp, y = µ)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ rep.id, ncol = 6) +
  theme_classic()

###### Run the models ######

c.summ.df <- data.frame(    # We'll create a dataframe to store the data as we fit models.
  block = numeric(),        # Block
  mic = character(),        # Microbial inocula
  rep = numeric(),          # Replicate
  id = character(),         # rep.id
  
  T.min = numeric(),        # Minimum T (calculus)
  T.min.min = numeric(),    # Minimum T (lower HDPI)
  T.min.max = numeric(),    # Minimum T (upper HDPI)
  T.min.na = numeric(),     # %% NA returns
  
  T.max = numeric(),        # Maximum T (calculus)
  T.max.min = numeric(),    # Maximum T (lower HDPI)
  T.max.max = numeric(),    # Maximum T (upper HDPI)
  T.max.na = numeric(),     # %% NA returns
  
  T.opt = numeric(),        # Optimal T (calculus)
  T.opt.min = numeric(),    # Optimal T (lower HDPI)
  T.opt.max = numeric(),    # Optimal T (upper HDPI)
  
  r.max = numeric(),        # Maximum growth rate (calculus)
  r.max.min = numeric(),    # Maximum growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum growth rate (upper HDPI)
  
  T.br.min = numeric(),     # T breadth (calculus)
  T.br.min.min = numeric(), # T breadth (lower HDPI)
  T.br.min.max = numeric(), # T breadth (upper HDPI)
  T.br.min.na = numeric(),  # T breadth na's
  
  T.br.max = numeric(),     # T breadth (calculus)
  T.br.max.min = numeric(), # T breadth (lower HDPI)
  T.br.max.max = numeric(), # T breadth (upper HDPI)
  T.br.max.na = numeric(),  # T breadth na's
  
  a = numeric(),            # parameter: a
  b = numeric(),            # parameter: b
  tmax = numeric(),         # parameter: tmax
  d.t = numeric(),          # parameter: deltaT
  
  a.mod = numeric(),        # nls.LM parameter: a
  b.mod = numeric(),        # nls.LM parameter: b
  tmax.mod = numeric(),     # nls.LM parameter: tmax
  d.t.mod = numeric(),      # nls.LM parameter: deltaT
  
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

parameters.lactin2 <- c("a", "b", "tmax", "d.t", "sigma", "r.pred") # parameters to estimate

Temp.xs <- seq(0, 50, 0.05) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

inits.lactin <- function() { # The other static initial values. 
  list(
    a = runif(1, 0.05, 0.15),  # More constrained initial values
    tmax = runif(1, 37, 43),
    d.t = runif(1, 1, 5),
    b = runif(1, -2.5, -1),
    sigma = runif(1, 0.1, 2)
  )
}

rep.ids <- unique(df.c.t$rep.id) # Save the unique rep ids so we can run the foor loop in chunks if needed

n <- 0 # for tracking progress

for (i in unique(df.c.t$rep.id[df.c.t$rep.id >= 6])) { # for each replicate ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.c.t %>% 
    filter(rep.id == i)
  
  df.i <- droplevels(df.i)
  
  trait <- df.i$µ    # format the data for jags
  N.obs <- length(trait)
  temp <- df.i$temp
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs) # assemble the jag data
  
  lac.jag <- jags(  # run the model
    data = jag.data, 
    inits = inits.lactin,
    parameters.to.save = parameters.lactin2, 
    model.file = "lactin.txt",
    n.thin = nt.fit, 
    n.chains = nc.fit, 
    n.burnin = nb.fit, 
    n.iter = ni.fit, 
    DIC = TRUE, 
    working.directory = getwd()
  )
  
  save(lac.jag, file = paste0("R2jags-models/rep_", i, "_lactin2.RData")) # save the lactin2 model
  # This folder is listed in gitignore, because the objects are too big to load
  
  post <- as.data.frame(lac.jag$BUGSoutput$sims.matrix) # The posterior distributions
  
  post <- post %>% 
    select(a, b, d.t, tmax)
  
  post.a <- median(post$a, na.rm = T)    # Extract parameters
  post.b <- median(post$b, na.rm = T)
  post.tmax <- median(post$tmax, na.rm = T)
  post.d.t <- median(post$d.t, na.rm = T)
  
  # Calculate summary metrics for the curve across all posteriors
  
  # Topt
  calc_Topt <- function(a, b, tmax, d.t) {
    tryCatch(
      uniroot(
        function(temp) lactin2_deriv(temp, a, b, tmax, d.t),
        interval = c(-10, 45)
      )$root,
      error = function(e) NA
    )
  }
  
  T.opt <- mapply(
    calc_Topt,
    a = post$a,
    b = post$b,
    tmax = post$tmax,
    d.t = post$d.t
  )
  
  #rmax
  r.max <- mapply(
    lactin2,
    temp = T.opt,
    a = post$a,
    b = post$b,
    tmax = post$tmax,
    d.t = post$d.t
  )
  
  #Tmin
  Tmin.safe <- function(Topt, a, b, tmax, d.t, # Tmin equation
                        lower = -50) {
    
    f_low  <- lactin2(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, b, tmax, d.t),
      interval = c(lower, Topt)
    )$root
  }
  
  T.min <- mapply(
    Tmin.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t
  )
  
  #Tmax
  Tmax.safe <- function(Topt, a, b, tmax, d.t, # Tmax equation
                        upper = 50) {
    
    f_low  <- lactin2(temp = Topt, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = upper,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, b, tmax, d.t),
      interval = c(Topt, upper)
    )$root
  }
  
  T.max <- mapply(
    Tmax.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t
  )
  
  # Tbr
  r.half <- r.max/2
  
  Tmin.half.safe <- function(Topt, a, b, tmax, d.t, r_half, # Tmin equation
                             lower = -10) {
    
    f_low  <- lactin2_halfmax(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    f_high <- lactin2_halfmax(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2_halfmax(temp, a, b, tmax, d.t, r_half),
      interval = c(lower, Topt)
    )$root
  }
  
  T.min.half <- mapply(
    Tmin.half.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t,
    r_half = r.half
  )
  
  Tmax.half.safe <- function(Topt, a, b, tmax, d.t, r_half, # Tmax equation
                             upper = 50) {
    
    f_low  <- lactin2_halfmax(temp = Topt, a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    f_high <- lactin2_halfmax(temp = upper,  a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2_halfmax(temp, a, b, tmax, d.t, r_half),
      interval = c(Topt, upper)
    )$root
  }
  
  T.max.half <- mapply(
    Tmax.half.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t,
    r_half = r.half
  )
  
  c.summ.df <- rbind(c.summ.df, data.frame(                                     # Add summary data
    block = df.i$block[1],                                                      # Block
    mic = df.i$mic[1],                                                          # Microbial inocula
    rep = df.i$rep[1],                                                          # Replicate
    id = df.i$rep.id[1],                                                        # Replicate id
    
    T.min = median(T.min, na.rm = T),                                           # Minimum T (calculus)
    T.min.min = hdi(T.min[!is.na(T.min)], credMass = 0.95)$CI_low,              # Minimum T (lower HDPI)
    T.min.max = hdi(T.min[!is.na(T.min)], credMass = 0.95)$CI_high,             # Minimum T (upper HDPI)
    T.min.na = mean(is.na(T.min)),                                              # % NA returns
    
    T.max = median(T.max, na.rm = T),                                           # Maximum T (calculus)
    T.max.min = hdi(T.max[!is.na(T.max)], credMass = 0.95)$CI_low,              # Maximum T (lower HDPI)
    T.max.max = hdi(T.max[!is.na(T.max)], credMass = 0.95)$CI_high,             # Maximum T (upper HDPI)
    T.max.na = mean(is.na(T.max)),                                              # % NA returns
    
    T.opt = median(T.opt, na.rm = T),                                           # Optimal T (calculus)
    T.opt.min = hdi(T.opt[!is.na(T.opt)], credMass = 0.95)$CI_low,              # Optimal T (lower HDPI)
    T.opt.max = hdi(T.opt[!is.na(T.opt)], credMass = 0.95)$CI_high,             # Optimal T (upper HDPI)
    
    r.max = median(r.max, na.rm = T),                                           # Maximum growth rate  (calculus)
    r.max.min = hdi(r.max[!is.na(r.max)], credMass = 0.95)$CI_low,              # Maximum growth rate  (lower HDPI)
    r.max.max = hdi(r.max[!is.na(r.max)], credMass = 0.95)$CI_high,             # Maximum growth rate  (upper HDPI)
    
    T.br.min = median(T.min.half, na.rm = T),                                   # T breadth (calculus)
    T.br.min.min = hdi(T.min.half[!is.na(T.min.half)], credMass = 0.95)$CI_low, # T breadth (lower HDPI)
    T.br.min.max = hdi(T.min.half[!is.na(T.min.half)], credMass = 0.95)$CI_high,# T breadth (upper HDPI)
    T.br.min.na = mean(is.na(T.min.half)),                                      # % NA returns
    
    T.br.max = median(T.max.half, na.rm = T),                                   # T breadth (calculus)
    T.br.max.min = hdi(T.max.half[!is.na(T.max.half)], credMass = 0.95)$CI_low, # T breadth (lower HDPI)
    T.br.max.max = hdi(T.max.half[!is.na(T.max.half)], credMass = 0.95)$CI_high,# T breadth (upper HDPI)
    T.br.max.na = mean(is.na(T.max.half)),                                      # % NA returns
    
    a = post.a,                                                                 # parameter: a
    b = post.b,                                                                 # parameter: b
    tmax = post.tmax,                                                           # parameter: tmax
    d.t = post.d.t,                                                             # parameter: deltaT
    
    a.mod = lac.jag$BUGSoutput$summary[1,1],                                    # Jags parameter: a
    b.mod = lac.jag$BUGSoutput$summary[2,1],                                    # Jags parameter: b
    tmax.mod = lac.jag$BUGSoutput$summary[1007,1],                               # Jags parameter: tmax
    d.t.mod =  lac.jag$BUGSoutput$summary[3,1]                                  # Jags parameter: deltaT
  ))
  
  lac.sum <- as.data.frame(lac.jag$BUGSoutput$summary[c(1:3,1007),])
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(                                         # Model performance data
      block = df.i$block[1],                                                    # Block
      mic = df.i$mic[1],                                                        # Microbial inocula
      rep = df.i$rep[1],                                                        # Replicate
      id = df.i$rep.id[1],                                                      # Replicate id
      
      Parameter = rownames(lac.sum)[j],                                         # Model parameter (e.g. cf.a, cf.tmax, etc.)
      mean = lac.sum[j,1],                                                      # Posterior mean
      Rhat = lac.sum[j,8],                                                      # Rhat values
      n.eff = lac.sum[j,9],                                                     # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE                   
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.c.t$rep.id))))
  
}

write.csv(c.summ.df, "processed-data/10a_chlamy_TPCs_bayes.csv") # 204 TPCs!
write.csv(fit.df, "processed-data/10b_chlamy_TPC_bayes_fits.csv")

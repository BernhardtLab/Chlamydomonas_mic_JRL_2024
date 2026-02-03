# Jason R Laurich
# February 3rd, 2026

# We're going to start fitting thermal performance curves to our Chlamydomonas reinhardtii and microbial growth rate data.

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
library(bayestestR)

lactin2 <- function(temp, a, tmax, d.t, b) { 
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
} # OK we're going to modify the function to calculate T_breadth, based on a modified lactin.

# Chlamy alone ------------------------------------------------------------

###### Load & examine the data ######

df.c <- read.csv("processed-data/02_chlamy_µs.csv")

head(df.c) # OK so we have unique.id, with temp, nit, salt, block, mic, plate, well, rep, and µ

df.c.t <- df.c %>% 
  filter(salt == 0,
         nit == 1000)

df.c.t <- df.c.t %>%
  mutate(rep.id = paste0("b", block, ".m", mic, ".r", rep)) # So this column will capture unique id's for each termperature gradient. 

length(unique(df.c.t$rep.id)) # 205 TPCs. 

df.c.t <- df.c.t %>% # Don't have data for the 4th replicates at all levels of temperature
  filter(rep <4)

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

n <-0 # progression tracker

for (i in unique(df.c.t$rep.id[df.c.t$rep.id >= 1])) { # for each replicate ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.c.t %>% 
    filter(rep.id == i)
  df.i <- droplevels(df.i)
  
  lac_nls <- nls_multstart(µ ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = c(0.01, -2.5, 36, 1),
                           start_upper = c(0.19, -0.5, 44, 5),
                           lower = c(0, -3, 30, 0.1),
                           upper = c(0.5, -0.001, 50, 6),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  sum <- summary(lac_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  lac.a <- df.nls[1,1] # Extract parameters
  lac.b <- df.nls[2,1]
  lac.tmax <- df.nls[3,1]
  lac.d.t <- df.nls[4,1]
  
  lac.LM <- nlsLM(µ ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                  data = df.i,
                  start = c(a = lac.a, b = lac.b, tmax = lac.tmax, delta_t = lac.d.t),
                  lower = c(a = max(0.001,lac.a - 0.05), b = lac.b - 0.5, tmax = max(0, lac.tmax -3), delta_t = max(0.1, lac.d.t - 2)),
                  upper = c(a = lac.a + 0.05, b = min(lac.b + 0.5, -0.001), tmax = lac.tmax + 3, delta_t = lac.d.t + 2),
                  control = nls.lm.control(maxiter=500)
  )
  
  sum <- summary(lac.LM)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  boot <- Boot(lac.LM, method = 'residual')
  
  post <- as.data.frame(boot$t) # Get the bootstrapped values
  post$d.t <- post$delta_t
  
  boot.a <- median(post$a, na.rm = T)    # Extract parameters
  boot.b <- median(post$b, na.rm = T)
  boot.tmax <- median(post$tmax, na.rm = T)
  boot.d.t <- median(post$d.t, na.rm = T)
  
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
                        lower = -100) {
    
    f_low  <- lactin2(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, tmax, d.t, b),
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
      function(temp) lactin2(temp, a, tmax, d.t, b),
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
    
    a = boot.a,                                                                 # parameter: a
    b = boot.b,                                                                 # parameter: b
    tmax = boot.tmax,                                                           # parameter: tmax
    d.t = boot.d.t,                                                             # parameter: deltaT
    
    a.mod = df.nls[1,1],                                                        # nls.LM parameter: a
    b.mod = df.nls[2,1],                                                        # nls.LM parameter: b
    tmax.mod = df.nls[3,1],                                                     # nls.LM parameter: tmax
    d.t.mod = df.nls[4,1]                                                       # nls.LM parameter: deltaT
  ))
  
  for (j in 1:4){
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
  
  print(paste("Done", n, "of ", length(unique(df.c.t$rep.id))))
  
}

write.csv(c.summ.df, "processed-data/05a_chlamy_TPCs.csv") # 204 TPCs!
write.csv(fit.df, "processed-data/05b_chlamy_TPC_fits.csv")

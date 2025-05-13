# Jason R Laurich
# May 13th, 2025

# I'll be fitting some TPC, Monod curves, and reversed logistic growth curve models to the µ data from the first block here
# Just to explore for now, I need more replication and I still need to account for some errors in the raw data for block 1.


# Load packages -----------------------------------------------------------

library(dplyr)
library(cowplot)
library(R2jags)
library(mcmcplots) 
library(gridExtra)
library(Deriv)

# Load the data -----------------------------------------------------------

df.T <- read.csv("processed-data/2a_blk1_temp_mus.csv")

head(df.T)
str(df.T)
df.T <- df.T %>% 
  mutate(mic = fct_relevel(mic, 
                               "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                               "9", "10", "11", "12", "13", "14", "15", "all"))

df.N <- read.csv("processed-data/2b_blk1_nit_mus.csv")

head(df.N)
str(df.N)
df.N <- df.N %>% 
  mutate(mic = fct_relevel(mic, 
                           "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                           "9", "10", "11", "12", "13", "14", "15", "all"))

df.S <- read.csv("processed-data/2c_blk1_salt_mus.csv")

head(df.S)
str(df.S)
df.S <- df.S %>% 
  mutate(mic = fct_relevel(mic, 
                           "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                           "9", "10", "11", "12", "13", "14", "15", "all"))

# We're going to cut this all down for now, just to look at some variation among microbes' effects on TPCs, monods and salt tolerance quickly
# 1/10th of the normal amounts, half the chains
ni.fit <- 33000   # iterations / chain
nb.fit <- 3000    # burn in periods for each chain
nt.fit <- 30      # thinning interval : (33,000 - 3,000) / 30 = 1000 posterior estimates / chain
nc.fit <- 3       # number of chains, total of 3,000 estimates for each model. 

# TPCs --------------------------------------------------------------------

inits.lactin.final <- function() { # The final initial values set we landed on after experimenting. 
  list(
    cf.a = runif(1, 0.05, 0.15),  # More constrained initial values
    cf.tmax = runif(1, 37, 43),
    cf.delta_t = runif(1, 1, 5),
    cf.b = runif(1, -2.5, -1),
    cf.sigma = runif(1, 0.1, 2)
  )
}

parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred") # what to fit/estimate and save

Temp.xs <- seq(0, 45, 0.05) # Temperature gradient
N.Temp.xs <-length(Temp.xs)

for (i in levels(df.T$mic)){
  
  df.i <- df.T %>% 
    filter(mic == i) %>% 
    droplevels()
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  temp <- df.i$temp
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  lac_jag <- jags(
    data = jag.data, 
    inits = inits.lactin.final, 
    parameters.to.save = parameters.lactin2, 
    model.file = "lactin2_final.txt",
    n.thin = nt.fit, 
    n.chains = nc.fit, 
    n.burnin = nb.fit, 
    n.iter = ni.fit, 
    DIC = TRUE, 
    working.directory = getwd()
  )
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:6),]
  df.jags.plot$temp <- seq(0, 45, 0.05)
  assign(paste0("df.T.jags.", i), df.jags.plot)
  
  message("Done ", i)
  
}

# Now let's plot all of the data!

p.t <- ggplot(df.T, aes(x = temp, y = r.exp, colour = mic)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  
  scale_color_manual(
    name = "Microbial treatment",  # Update the legend title
    values = c("none" = "black",
               "1" = "mediumorchid4",
               "2" = "magenta2",
               "3" = "tan",
               "4" = "turquoise1",
               "5" = "dodgerblue",
               "6" = "royalblue4",  
               "7" = "blue",
               "8" = "deeppink",
               "9" = "darkorange2",
               "10" = "orange",
               "11" = "firebrick",
               "12" = "orangered2",
               "13" = "red1",  
               "14" = "forestgreen",
               "15" = "darkolivegreen2",
               "all" = "goldenrod1")
  ) +
  
  geom_line(data = df.T.jags.none, aes(x = temp, y= mean), colour = "black", linewidth = 1) +
  geom_line(data = df.T.jags.1, aes(x = temp, y= mean), colour = "mediumorchid4", linewidth = 1) +
  geom_line(data = df.T.jags.2, aes(x = temp, y= mean), colour = "magenta2", linewidth = 1) +
  geom_line(data = df.T.jags.3, aes(x = temp, y= mean), colour = "tan", linewidth = 1) +
  geom_line(data = df.T.jags.4, aes(x = temp, y= mean), colour = "turquoise1", linewidth = 1) +
  geom_line(data = df.T.jags.5, aes(x = temp, y= mean), colour = "dodgerblue", linewidth = 1) +
  geom_line(data = df.T.jags.6, aes(x = temp, y= mean), colour = "royalblue4", linewidth = 1) +
  geom_line(data = df.T.jags.7, aes(x = temp, y= mean), colour = "blue", linewidth = 1) +
  geom_line(data = df.T.jags.8, aes(x = temp, y= mean), colour = "deeppink", linewidth = 1) +
  geom_line(data = df.T.jags.9, aes(x = temp, y= mean), colour = "darkorange", linewidth = 1) +
  geom_line(data = df.T.jags.10, aes(x = temp, y= mean), colour = "orange", linewidth = 1) +
  geom_line(data = df.T.jags.11, aes(x = temp, y= mean), colour = "firebrick", linewidth = 1) +
  geom_line(data = df.T.jags.12, aes(x = temp, y= mean), colour = "orangered2", linewidth = 1) +
  geom_line(data = df.T.jags.13, aes(x = temp, y= mean), colour = "red1", linewidth = 1) +
  geom_line(data = df.T.jags.14, aes(x = temp, y= mean), colour = "forestgreen", linewidth = 1) +
  geom_line(data = df.T.jags.15, aes(x = temp, y= mean), colour = "darkolivegreen2", linewidth = 1) +
  geom_line(data = df.T.jags.all, aes(x = temp, y= mean), colour = "goldenrod1", linewidth = 1) +
  
  ylim(-0.5, 3.5) +
  xlim(0, 45) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "Thermal performance curves (Lactin II)") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.t

# Nitrogen Monod curves ---------------------------------------------------

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

S.pred <- seq(0, 1000, 0.5) # Nitrogen gradient
N.S.pred <-length(S.pred)

for (i in levels(df.N$mic)){
  
  df.i <- df.N %>% 
    filter(mic == i) %>% 
    droplevels()
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  nit <- df.i$nit
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = nit, S.pred = S.pred, N.S.pred = N.S.pred)
  
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
  
  df.jags.plot <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
  df.jags.plot$nit <- seq(0, 1000, 0.5)
  assign(paste0("df.N.jags.", i), df.jags.plot)
  
  message("Done ", i)
  
}

p.n <- ggplot(df.N, aes(x = nit, y = r.exp, colour = mic)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  
  scale_color_manual(
    name = "Microbial treatment",  # Update the legend title
    values = c("none" = "black",
               "1" = "mediumorchid4",
               "2" = "magenta2",
               "3" = "tan",
               "4" = "turquoise1",
               "5" = "dodgerblue",
               "6" = "royalblue4",  
               "7" = "blue",
               "8" = "deeppink",
               "9" = "darkorange2",
               "10" = "orange",
               "11" = "firebrick",
               "12" = "orangered2",
               "13" = "red1",  
               "14" = "forestgreen",
               "15" = "darkolivegreen2",
               "all" = "goldenrod1")
  ) +
  
  geom_line(data = df.N.jags.none, aes(x = nit, y= mean), colour = "black", linewidth = 1) +
  geom_line(data = df.N.jags.1, aes(x = nit, y= mean), colour = "mediumorchid4", linewidth = 1) +
  geom_line(data = df.N.jags.2, aes(x = nit, y= mean), colour = "magenta2", linewidth = 1) +
  geom_line(data = df.N.jags.3, aes(x = nit, y= mean), colour = "tan", linewidth = 1) +
  geom_line(data = df.N.jags.4, aes(x = nit, y= mean), colour = "turquoise1", linewidth = 1) +
  geom_line(data = df.N.jags.5, aes(x = nit, y= mean), colour = "dodgerblue", linewidth = 1) +
  geom_line(data = df.N.jags.6, aes(x = nit, y= mean), colour = "royalblue4", linewidth = 1) +
  geom_line(data = df.N.jags.7, aes(x = nit, y= mean), colour = "blue", linewidth = 1) +
  geom_line(data = df.N.jags.8, aes(x = nit, y= mean), colour = "deeppink", linewidth = 1) +
  geom_line(data = df.N.jags.9, aes(x = nit, y= mean), colour = "darkorange", linewidth = 1) +
  geom_line(data = df.N.jags.10, aes(x = nit, y= mean), colour = "orange", linewidth = 1) +
  geom_line(data = df.N.jags.11, aes(x = nit, y= mean), colour = "firebrick", linewidth = 1) +
  geom_line(data = df.N.jags.12, aes(x = nit, y= mean), colour = "orangered2", linewidth = 1) +
  geom_line(data = df.N.jags.13, aes(x = nit, y= mean), colour = "red1", linewidth = 1) +
  geom_line(data = df.N.jags.14, aes(x = nit, y= mean), colour = "forestgreen", linewidth = 1) +
  geom_line(data = df.N.jags.15, aes(x = nit, y= mean), colour = "darkolivegreen2", linewidth = 1) +
  geom_line(data = df.N.jags.all, aes(x = nit, y= mean), colour = "goldenrod1", linewidth = 1) +
  
  ylim(-0.5, 2.5) +
  xlim(-10, 1010) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Nitrate concentration (µM)", 
       y = "Exponential growth rate (µ)",
       title = "Monod curves (Nitrogen gradient)") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.n

# Salt reversed logistic growth curves ------------------------------------

inits.salt <- function() {
  list(
    a = runif(1, 0.1, 5),  # Initial guess for a
    b = runif(1, 0.1, 5),  # Initial guess for b
    c = runif(1, 0.1, max(df.i$salt)),  # Initial guess for c
    sigma = runif(1, 0.1, 2)  # Initial guess for error
  )
}

parameters.salt <- c("a", "b", "c", "sigma", "r_pred_new") # Save these

S.pred <- seq(0, 12, 0.005) # Salt gradient we are interested in we'll keep N.S.pred more or less consistent across abiotic gradients for now
N.S.pred <-length(S.pred)

for (i in levels(df.S$mic)){
  
  df.i <- df.S %>% 
    filter(mic == i) %>% 
    droplevels()
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  salt <- df.i$salt
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the salt logistic growth curve function. 
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
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:4,2406),]
  df.jags.plot$salt <- seq(0, 12, 0.005)
  assign(paste0("df.S.jags.", i), df.jags.plot)
  
  message("Done ", i)
  
}

# Now let's plot all of the data!

p.s <- ggplot(df.S, aes(x = salt, y = r.exp, colour = mic)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  
  scale_color_manual(
    name = "Microbial treatment",  # Update the legend title
    values = c("none" = "black",
               "1" = "mediumorchid4",
               "2" = "magenta2",
               "3" = "tan",
               "4" = "turquoise1",
               "5" = "dodgerblue",
               "6" = "royalblue4",  
               "7" = "blue",
               "8" = "deeppink",
               "9" = "darkorange2",
               "10" = "orange",
               "11" = "firebrick",
               "12" = "orangered2",
               "13" = "red1",  
               "14" = "forestgreen",
               "15" = "darkolivegreen2",
               "all" = "goldenrod1")
  ) +
  
  geom_line(data = df.S.jags.none, aes(x = salt, y= mean), colour = "black", linewidth = 1) +
  geom_line(data = df.S.jags.1, aes(x = salt, y= mean), colour = "mediumorchid4", linewidth = 1) +
  geom_line(data = df.S.jags.2, aes(x = salt, y= mean), colour = "magenta2", linewidth = 1) +
  geom_line(data = df.S.jags.3, aes(x = salt, y= mean), colour = "tan", linewidth = 1) +
  geom_line(data = df.S.jags.4, aes(x = salt, y= mean), colour = "turquoise1", linewidth = 1) +
  geom_line(data = df.S.jags.5, aes(x = salt, y= mean), colour = "dodgerblue", linewidth = 1) +
  geom_line(data = df.S.jags.6, aes(x = salt, y= mean), colour = "royalblue4", linewidth = 1) +
  geom_line(data = df.S.jags.7, aes(x = salt, y= mean), colour = "blue", linewidth = 1) +
  geom_line(data = df.S.jags.8, aes(x = salt, y= mean), colour = "deeppink", linewidth = 1) +
  geom_line(data = df.S.jags.9, aes(x = salt, y= mean), colour = "darkorange", linewidth = 1) +
  geom_line(data = df.S.jags.10, aes(x = salt, y= mean), colour = "orange", linewidth = 1) +
  geom_line(data = df.S.jags.11, aes(x = salt, y= mean), colour = "firebrick", linewidth = 1) +
  geom_line(data = df.S.jags.12, aes(x = salt, y= mean), colour = "orangered2", linewidth = 1) +
  geom_line(data = df.S.jags.13, aes(x = salt, y= mean), colour = "red1", linewidth = 1) +
  geom_line(data = df.S.jags.14, aes(x = salt, y= mean), colour = "forestgreen", linewidth = 1) +
  geom_line(data = df.S.jags.15, aes(x = salt, y= mean), colour = "darkolivegreen2", linewidth = 1) +
  geom_line(data = df.S.jags.all, aes(x = salt, y= mean), colour = "goldenrod1", linewidth = 1) +
  
  ylim(-0.5, 2.5) +
  xlim(0, 12) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Salt (g/L)", 
       y = "Exponential growth rate (µ)",
       title = "Reversed logistic growth curves (Salt tolerance)") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.s


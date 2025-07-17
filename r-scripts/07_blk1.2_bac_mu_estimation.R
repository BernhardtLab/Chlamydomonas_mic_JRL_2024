# Jason R Laurich
# July 13th, 2025

# Going to estimate exponential growth rate (µ) for each well of the first block of my Chlamydomonas x microbes x global change
# experiment (Spring 2025, blocks 1 and 2 - no carbon added to the wells)

# Estimating this for bacteria-only wells and bacterial performance alongside Chlamy

# Load packages -----------------------------------------------------------

library(nls.multstart)
library(tidyr)
library(cowplot)
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)

# Load and explore the data -----------------------------------------------------------

df.1 <- read.csv("processed-data/01_blk1_rawdata.csv") # Get the raw data for block 1
df.1$block <- 1

df.2 <- read.csv("processed-data/03_blk2_rawdata.csv") # Get the raw data for block 2
df.2$block <- 2

df <- rbind (df.1, df.2)

# We will be working with OD600 data here, but this needs to be transformed into cell counts. First, we need to backcalculate
# Chlamy's contribution to OD600 in algae-only wells and then remove that portion from shared wells.

df %>%
  filter(Chlamy.y.n == "y", Microbe == "none") %>%
  ggplot(aes(x = RFU, y = OD600)) +
  geom_point(alpha = 0.5) +
  labs(
    title = "Relationship between RFU and OD600",
    x = "RFU",
    y = "OD600"
  ) +
  theme_classic()

# How much is day explaining?

df %>%
  filter(Chlamy.y.n == "y", Microbe == "none") %>%
  ggplot(aes(x = RFU, y = OD600, colour = days)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_viridis_c(option = "D", direction = -1) +
  labs(
    title = "Relationship between RFU and OD600 over time",
    x = "RFU",
    y = "OD600",
    colour = "Days"
  ) +
  theme_classic()

# Let's remove the early points (inoculation)?

df %>%
  filter(Chlamy.y.n == "y", Microbe == "none", days >= 1) %>%
  ggplot(aes(x = RFU, y = OD600, colour = days)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_viridis_c(option = "D", direction = -1) +
  labs(
    title = "Relationship between RFU and OD600 over time",
    x = "RFU",
    y = "OD600",
    colour = "Days"
  ) +
  theme_classic()

# How much are salt and nitrogen explaining?

df %>%
  filter(Chlamy.y.n == "y", Microbe == "none", days >= 1) %>%
  ggplot(aes(x = RFU, y = OD600, colour = Nitrogen.conc.µM)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_viridis_c(option = "D", direction = -1) +
  labs(
    title = "Relationship between RFU and OD600 over time",
    x = "RFU",
    y = "OD600",
    colour = "Nitrogen"
  ) +
  theme_classic()

df %>%
  filter(Chlamy.y.n == "y", Microbe == "none", days >= 1) %>%
  ggplot(aes(x = RFU, y = OD600, colour = Salt.conc.g.l)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_viridis_c(option = "D", direction = -1) +
  labs(
    title = "Relationship between RFU and OD600 over time",
    x = "RFU",
    y = "OD600",
    colour = "Salt"
  ) +
  theme_classic()

# OK let's just add a line to this plot:

df %>%
  filter(Chlamy.y.n == "y", Microbe == "none") %>%
  ggplot(aes(x = RFU, y = OD600, colour = days)) +
  geom_smooth(method='lm') +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_viridis_c(option = "D", direction = -1) +
  labs(
    title = "Relationship between RFU and OD600 over time",
    x = "RFU",
    y = "OD600",
    colour = "Days"
  ) +
  theme_classic()

# Calculating bacterial cell density --------------------------------------

mod <- lm(OD600 ~ RFU, data= df[df$Chlamy.y.n == 'y' & df$Microbe == 'none', ])
C.OD600 <- coef(mod)["RFU"]

df <- df %>%
  mutate(
    mic.prod = case_when(
      Chlamy.y.n == "n" ~ OD600 * 8e8,
      Chlamy.y.n == "y" ~ pmax(OD600 - (RFU * C.OD600), 0) * 8e8,
      TRUE ~ NA_real_
    )
  )

df %>%
  ggplot(aes(x = Microbe, y = mic.prod, colour = Chlamy.y.n)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_classic()

df %>%
  ggplot(aes(x = Microbe, y = mic.prod, colour = Chlamy.y.n)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.5), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.5), width = 0.2) +
  labs(y = "Microbial production (cells)", colour = "Chlamy present?") +
  theme_classic()

# OK I don't trust these data - for now we will work only with the microbe only wells. 
# No let's run it with both!

df <- df %>% 
  filter(Microbe != 'BLANK') %>%
  droplevels() %>%
  mutate(Microbe = fct_relevel(Microbe, 
                               "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                               "9", "10", "11", "12", "13", "14", "15", "all")) %>% 
  mutate(well.ID = paste(Plate, Row, Column, sep = ".")) %>% 
  mutate(logRFU = log(RFU + 0.0001)) # Let's reorder and clean these out

df.T <- df %>% 
  filter(Salt.conc.g.l == 0, Nitrogen.conc.µM == 1000) 

df.N <- df %>% 
  filter(Salt.conc.g.l == 0, Temperature.C == 30)

df.S <- df %>% 
  filter(Nitrogen.conc.µM == 1000, Temperature.C == 30)

# Temperature -------------------------------------------------------------

df.bac.r.exp.t <- data.frame( # Initializing a dataframe to store the results for each well, microbe and temp level
  block = numeric(),
  Chlamy = character(),
  mic = character(),
  temp = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

for (i in levels(df.T$Microbe)){
  
  for (t in unique(df.T$Temperature.C)){
    
    df.it <- df.T %>% # focus on the microbe and temperature, drop unneeded levels
      filter(Temperature.C == t, Microbe == i) %>% 
      droplevels()
    
    for (w in unique(df.it$well.ID)){
      
      df.it.w <- subset(df.it, df.it$well.ID == w)
      
      df.it.w <- df.it.w[order(df.it.w$days), ] # Order by days just in case something weird happened
      
      if (df.it.w$RFU[2] <  df.it.w$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.w <- df.it.w[-1,]
      }
      
      df.it.w$N0 <- df.it.w$RFU[1] # Set the initial value to the first (lowest)
      
      t.series <- unique(df.it.w$days) # Re-initialize this internally - we will only save summary data for each unique pop x salt x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.w.sl <- df.it.w[df.it.w$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.w.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID at microbe x temperature level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 10%
        
        if (ln.slopes[r] != 0 & round(ln.slopes[r], digits = 5) !=0) { # We also need to account for tiny values that are basically 0 (e.g. 5 e-16, but are messing up our loops)
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] 
          
          if (percent.chg >= 1.10 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.w.th <- df.it.w[df.it.w$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.w.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.bac.r.exp.t <- rbind(df.bac.r.exp.t, data.frame(
          block = df.it.w$block[1],
          Chlamy = df.it.w$Chlamy.y.n[1],
          mic = df.it.w$Microbe[1],
          temp = df.it.w$Temperature.C[1],
          well.ID = df.it.w$well.ID[1],
          r.exp = NA
        ))
        
      }else{
        df.bac.r.exp.t <- rbind(df.bac.r.exp.t, data.frame(
          block = df.it.w$block[1],
          Chlamy = df.it.w$Chlamy.y.n[1],
          mic = df.it.w$Microbe[1],
          temp = df.it.w$Temperature.C[1],
          well.ID = df.it.w$well.ID[1],
          r.exp = summary(r_exp)$parameters[1,1]
        ))
        
        
      }
      
    }
    
  }
  
}

write.csv(df.bac.r.exp.t, "processed-data/6a_blk1.2_mic_temp_mus.csv") # let's save the file.

# Nitrogen -------------------------------------------------------------

df.bac.r.exp.n <- data.frame( # Initializing a dataframe to store the results for each well, microbe and temp level
  block = numeric(),
  Chlamy = character(),
  mic = character(),
  nit = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

for (i in levels(df.N$Microbe)){
  
  for (t in unique(df.N$Nitrogen.conc.µM)){
    
    df.it <- df.N %>% # focus on the microbe and nitrogen, drop unneeded levels
      filter(Nitrogen.conc.µM == t, Microbe == i) %>% 
      droplevels()
    
    for (w in unique(df.it$well.ID)){
      
      df.it.w <- subset(df.it, df.it$well.ID == w)
      
      df.it.w <- df.it.w[order(df.it.w$days), ] # Order by days just in case something weird happened
      
      if (df.it.w$RFU[2] <  df.it.w$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.w <- df.it.w[-1,]
      }
      
      df.it.w$N0 <- df.it.w$RFU[1] # Set the initial value to the first (lowest)
      
      t.series <- unique(df.it.w$days) # Re-initialize this internally - we will only save summary data for each unique pop x salt x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.w.sl <- df.it.w[df.it.w$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.w.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID at microbe x temperature level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 10%
        
        if (ln.slopes[r] != 0 & round(ln.slopes[r], digits = 5) !=0) { # We also need to account for tiny values that are basically 0 (e.g. 5 e-16, but are messing up our loops)
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] 
          
          if (percent.chg >= 1.10 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.w.th <- df.it.w[df.it.w$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.w.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.bac.r.exp.n <- rbind(df.bac.r.exp.n, data.frame(
          block = df.it.w$block[1],
          Chlamy = df.it.w$Chlamy.y.n[1],
          mic = df.it.w$Microbe[1],
          nit = df.it.w$Nitrogen.conc.µM[1],
          well.ID = df.it.w$well.ID[1],
          r.exp = NA
        ))
        
      }else{
        df.bac.r.exp.n <- rbind(df.bac.r.exp.n, data.frame(
          block = df.it.w$block[1],
          Chlamy = df.it.w$Chlamy.y.n[1],
          mic = df.it.w$Microbe[1],
          nit = df.it.w$Nitrogen.conc.µM[1],
          well.ID = df.it.w$well.ID[1],
          r.exp = summary(r_exp)$parameters[1,1]
        ))
        
        
      }
      
    }
    
  }
  
}

write.csv(df.bac.r.exp.n, "processed-data/6b_blk1.2_mic_nit_mus.csv") # let's save the file.

# Salt -------------------------------------------------------------

df.bac.r.exp.s <- data.frame( # Initializing a dataframe to store the results for each well, microbe and temp level
  block = numeric(),
  Chlamy = character(),
  mic = character(),
  salt = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

for (i in levels(df.S$Microbe)){
  
  for (t in unique(df.S$Salt.conc.g.l)){
    
    df.it <- df.S %>% # focus on the microbe and nitrogen, drop unneeded levels
      filter(Salt.conc.g.l == t, Microbe == i) %>% 
      droplevels()
    
    for (w in unique(df.it$well.ID)){
      
      df.it.w <- subset(df.it, df.it$well.ID == w)
      
      df.it.w <- df.it.w[order(df.it.w$days), ] # Order by days just in case something weird happened
      
      if (df.it.w$RFU[2] <  df.it.w$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.w <- df.it.w[-1,]
      }
      
      df.it.w$N0 <- df.it.w$RFU[1] # Set the initial value to the first (lowest)
      
      t.series <- unique(df.it.w$days) # Re-initialize this internally - we will only save summary data for each unique pop x salt x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.w.sl <- df.it.w[df.it.w$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.w.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID at microbe x temperature level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 10%
        
        if (ln.slopes[r] != 0 & round(ln.slopes[r], digits = 5) !=0) { # We also need to account for tiny values that are basically 0 (e.g. 5 e-16, but are messing up our loops)
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] 
          
          if (percent.chg >= 1.10 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.w.th <- df.it.w[df.it.w$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.w.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.bac.r.exp.s <- rbind(df.bac.r.exp.s, data.frame(
          block = df.it.w$block[1],
          Chlamy = df.it.w$Chlamy.y.n[1],
          mic = df.it.w$Microbe[1],
          salt = df.it.w$Salt.conc.g.l[1],
          well.ID = df.it.w$well.ID[1],
          r.exp = NA
        ))
        
      }else{
        df.bac.r.exp.s <- rbind(df.bac.r.exp.s, data.frame(
          block = df.it.w$block[1],
          Chlamy = df.it.w$Chlamy.y.n[1],
          mic = df.it.w$Microbe[1],
          salt = df.it.w$Salt.conc.g.l[1],
          well.ID = df.it.w$well.ID[1],
          r.exp = summary(r_exp)$parameters[1,1]
        ))
        
        
      }
      
    }
    
  }
  
}

write.csv(df.bac.r.exp.s, "processed-data/6c_blk1.2_mic_salt_mus.csv") # let's save the file.

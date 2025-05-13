# Jason R Laurich
# May 12th, 2025

# Going to estimate exponential growth rate (µ) for each well of the first block of my Chlamydomonas x microbes x global change
# experiment (Spring 2025)


# Load packages -----------------------------------------------------------

library(nls.multstart)
library(tidyr)
library(cowplot)
library(ggplot2)
library(dplyr)

# Load data ---------------------------------------------------------------

df <- read.csv("processed-data/01_blk1_rawdata.csv") # Get the experimental design file.

df$Microbe <- as.factor(df$Microbe)
levels(df$Microbe) # Right so we have BLANK in there as an extra factor level. That's fine.

df <- df %>% 
  filter(Microbe != 'BLANK') %>%
  droplevels() %>%
  mutate(Microbe = fct_relevel(Microbe, 
                               "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                               "9", "10", "11", "12", "13", "14", "15", "all")) %>% 
  mutate(well.ID = paste(Plate, Row, Column, sep = ".")) %>% 
  mutate(logRFU = log(RFU + 0.0001)) # Let's reorder and clean these out
 
df.T <- df %>% 
  filter(Chlamy.y.n == 'y', Salt.conc.g.l == 0, Nitrogen.conc.µM == 1000) 

df.N <- df %>% 
  filter(Chlamy.y.n == 'y', Salt.conc.g.l == 0, Temperature.C == 30)

df.S <- df %>% 
  filter(Chlamy.y.n == 'y', Nitrogen.conc.µM == 1000, Temperature.C == 30)

# Temperature -------------------------------------------------------------

df.r.exp.t <- data.frame( # Initializing a dataframe to store the results for each well, microbe and temp level
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
        
        df.r.exp.t <- rbind(df.r.exp.t, data.frame(
          mic = df.it.w$Microbe[1],      
          temp = df.it.w$Temperature.C[1],        
          well.ID = df.it.w$well.ID[1],                
          r.exp = NA        
        ))
        
      }else{
        df.r.exp.t <- rbind(df.r.exp.t, data.frame(    # Add data to our summary table
          mic = df.it.w$Microbe[1],                    # Microbes
          temp = df.it.w$Temperature.C[1],             # Temperature
          well.ID = df.it.w$well.ID[1],                # Well ID
          r.exp = summary(r_exp)$parameters[1,1]       # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
}

write.csv(df.r.exp.t, "processed-data/2a_blk1_temp_mus.csv") # let's save the file.

# Nitrogen ----------------------------------------------------------------

df.r.exp.n <- data.frame( # Initializing a dataframe to store the results for each well, microbe and nitrogen level
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
        
        df.r.exp.n <- rbind(df.r.exp.n, data.frame(
          mic = df.it.w$Microbe[1],      
          nit = df.it.w$Nitrogen.conc.µM[1],        
          well.ID = df.it.w$well.ID[1],                
          r.exp = NA        
        ))
        
      }else{
        df.r.exp.n <- rbind(df.r.exp.n, data.frame(    # Add data to our summary table
          mic = df.it.w$Microbe[1],                    # Microbes
          nit = df.it.w$Nitrogen.conc.µM[1],           # Nitrogen
          well.ID = df.it.w$well.ID[1],                # Well ID
          r.exp = summary(r_exp)$parameters[1,1]       # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
}

write.csv(df.r.exp.n, "processed-data/2b_blk1_nit_mus.csv") # let's save the file.

# Salt --------------------------------------------------------------------

df.r.exp.s <- data.frame( # Initializing a dataframe to store the results for each well, microbe and salt level
  mic = character(),
  salt = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

for (i in levels(df.S$Microbe)){
  
  for (t in unique(df.S$Salt.conc.g.l)){
    
    df.it <- df.S %>% # focus on the microbe and salt, drop unneeded levels
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
        
        df.r.exp.s <- rbind(df.r.exp.s, data.frame(
          mic = df.it.w$Microbe[1],      
          salt = df.it.w$Salt.conc.g.l[1],        
          well.ID = df.it.w$well.ID[1],                
          r.exp = NA        
        ))
        
      }else{
        df.r.exp.s <- rbind(df.r.exp.s, data.frame(    # Add data to our summary table
          mic = df.it.w$Microbe[1],                    # Microbes
          salt = df.it.w$Salt.conc.g.l[1],             # Salt
          well.ID = df.it.w$well.ID[1],                # Well ID
          r.exp = summary(r_exp)$parameters[1,1]       # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
}

write.csv(df.r.exp.s, "processed-data/2c_blk1_salt_mus.csv") # let's save the file.


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

df <- df %>% # Let's reorder and clean these out
  filter(Microbe != 'BLANK') %>%
  droplevels() %>%
  mutate(Microbe = fct_relevel(Microbe, 
                               "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                               "9", "10", "11", "12", "13", "14", "15", "all"))

df.T <- df %>% 
  filter(Chlamy.y.n == 'y', Salt.conc.g.l == 0, Nitrogen.conc.µM == 1000) 

df.N <- df %>% 
  filter(Chlamy.y.n == 'y', Salt.conc.g.l == 0, Temperature.C == 30)

df.S <- df %>% 
  filter(Chlamy.y.n == 'y', Nitrogen.conc.µM == 100, Temperature.C == 30)



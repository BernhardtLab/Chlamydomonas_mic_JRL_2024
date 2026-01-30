# Jason R Laurich
# January 29th, 2026

# We're going to upload all of the raw data (RFUs) from the CSVs and save them into a single file with all of the timeseries data for every block.

# Load packages -----------------------------------------------------------

library(tidyverse)
library(zoo)
library(readxl)
library(lubridate)
library(openxlsx)
library(cowplot)
library(dplyr)

# Upload the design files for each block --------------------------------

df.design.1 <- read.csv("raw-data/02-block-1-design.csv") # Get the experimental design file.
head(df.design.1) # This matches treatment to wells
df.design.1$Block <- 1

df.design.1 <- df.design.1 %>% 
  mutate(unique.id = paste0("b1.t", Temperature.C, ".p", Plate.at.T, ".w", Well.at.T))

df.design.1 %>% 
  filter(!is.na(Notes), trimws(Notes) != "") # Let's identify the problematic rows — e.g. where I left notes on potential errors etc. 

#Hmm, for now I will keep samples where I caught pipetting errors. For example, we will keep the treatments where I noted that wells received bacteria 5, not 4
df.design.1 <- df.design.1 %>%
  
  mutate(
    Replicate = suppressWarnings(as.numeric(Replicate))
  ) %>%
  
  filter(!unique.id %in% c("b1.t25.p1.w19",
                          "b1.t33.p1.w2",
                          "b1.t33.p1.w12",
                          "b1.t39.p2.w105")) %>% 
  
  mutate(Microbe = case_when(
    unique.id %in% c("b1.t30.p27.w1615",
                     "b1.t30.p28.w1641",
                     "b1.t30.p28.w1649",
                     "b1.t30.p28.w1663",
                     "b1.t30.p28.w1671",
                     "b1.t30.p29.w1685",
                     "b1.t30.p29.w1699",
                     "b1.t30.p29.w1716",
                     "b1.t30.p29.w1717",
                     "b1.t30.p29.w1730",
                     "b1.t30.p29.w1735") ~ "5",
    TRUE ~ Microbe)) %>% 
                                              
  mutate(Replicate = case_when(
    unique.id %in% c("b1.t30.p27.w1615",
                     "b1.t30.p28.w1641",
                     "b1.t30.p28.w1649",
                     "b1.t30.p28.w1663",
                     "b1.t30.p28.w1671",
                     "b1.t30.p29.w1685",
                     "b1.t30.p29.w1699",
                     "b1.t30.p29.w1716",
                     "b1.t30.p29.w1717",
                     "b1.t30.p29.w1730") ~ 4,
    
    unique.id == "b1.t30.p29.w1735" ~ 5,
    TRUE ~ Replicate)) 


df.design.2 <- read.csv("raw-data/04-block-2-design.csv") # Get the experimental design file.
head(df.design.2) # This matches treatment to wells
df.design.2$Block <- 2

df.design.2 <- df.design.2 %>% 
  mutate(unique.id = paste0("b2.t", Temperature.C, ".p", Plate.at.T, ".w", Well.at.T))

df.design.2 %>% 
  filter(!is.na(Notes), trimws(Notes) != "") # Let's identify the problematic rows — e.g. where I left notes on potential errors etc. 

# OK so this is just a general note — this whole block didn't get carbon in its media. We'll use this as a comparison point (maybe in the supplement)
# To explore the effects of carbon enrichment on mutualism etc. Even if the comparison isn't perfect. 

df.design.3 <- read.csv("raw-data/06-block-3-design.csv") # Get the experimental design file.
head(df.design.3) # This matches treatment to wells
df.design.3$Block <- 3

df.design.3 <- df.design.3 %>% 
  mutate(unique.id = paste0("b3.t", Temperature.C, ".p", Plate.at.T, ".w", Well.at.T))

df.design.3 %>% 
  filter(!is.na(Notes), trimws(Notes) != "") # Let's identify the problematic rows — e.g. where I left notes on potential errors etc.

# We'll change the treatment information for the 2 errors. 

df.design.3 <- df.design.3 %>%
  
  mutate(
    Replicate = suppressWarnings(as.numeric(Replicate))
  ) %>%
  
  mutate(Microbe = case_when(
    unique.id == "b3.t30.p7.w365" ~ "8",
    unique.id == "b3.t30.p26.w1542" ~ "4",
    TRUE ~ Microbe)) %>% 
  
  mutate(Replicate = case_when(
    unique.id %in% c("b3.t30.p7.w365",
                     "b3.t30.p26.w1542") ~ 4,
    TRUE ~ Replicate)) 
    
df.design.4 <- read.csv("raw-data/08-block-4-design.csv") # Get the experimental design file.
head(df.design.4) # This matches treatment to wells
df.design.4$Block <- 4

df.design.4 <- df.design.4 %>% 
  mutate(unique.id = paste0("b4.t", Temperature.C, ".p", Plate.at.T, ".w", Well.at.T))

df.design.4 %>% 
  filter(!is.na(Notes), trimws(Notes) != "") # Let's identify the problematic rows — e.g. where I left notes on potential errors etc. 

# OK we'll fix these errors. I'm going to take out the one that got some unknown number of bacteria, as well as the one that got chlamy accidentally
df.design.4 <- df.design.4 %>%
  
  mutate(
    Replicate = suppressWarnings(as.numeric(Replicate))
  ) %>%
  
  filter(!unique.id %in% c("b4.t30.p19.w1085",
                           "b4.t30.p20.w1152")) %>% 
  
  mutate(Microbe = case_when(
    unique.id == "b4.t30.p30.w1764" ~ "10",
    TRUE ~ Microbe)) %>% 
  
  mutate(Replicate = case_when(
    unique.id == "b4.t30.p30.w1764" ~ 4,
    TRUE ~ Replicate)) 

# OK these should be good to go! Now I can start loading in the raw data csvs and extracting RFU, OD, and timestamp data. 

# Block 1 -----------------------------------------------------------------



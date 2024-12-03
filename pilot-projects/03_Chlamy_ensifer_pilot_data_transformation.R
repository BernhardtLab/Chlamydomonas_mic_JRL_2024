# Jason R Laurich
# December 3, 2024

# Script to process plate reader data for Chlamydomnas-Ensifer pilot project (6 time points)

############### Packages #############################

library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)

############### Access the csv files ######################

df.design <- read.csv("pilot-projects/04_Chlamy_Ensifer_pilot_design.csv") # Get the experimental design file.

head(df.design) # This matches treatment to wells
str(df.design)  

df.design <- df.design[,-c(12:15)]
df.design[, 5:7] <- lapply(df.design[, 5:7], as.factor)

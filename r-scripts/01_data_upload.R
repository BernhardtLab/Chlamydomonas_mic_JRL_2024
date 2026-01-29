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

df.design.2 <- read.csv("raw-data/04-block-2-design.csv") # Get the experimental design file.

head(df.design.2) # This matches treatment to wells


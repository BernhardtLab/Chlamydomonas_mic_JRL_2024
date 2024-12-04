# Jason R Laurich
# December 3, 2024

# Script to process plate reader data for Chlamydomnas-Ensifer pilot project (6 time points)

############### Packages #############################

library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)
library(readxl)
library(lubridate)
library(openxlsx)

############### Access the csv files ######################

df.design <- read.csv("pilot-projects/04_Chlamy_Ensifer_pilot_design.csv") # Get the experimental design file.

head(df.design) # This matches treatment to wells
str(df.design)  

df.design <- df.design[,-c(12:15)]
df.design[, 5:7] <- lapply(df.design[, 5:7], as.factor)

df <- df.design[,1:7]
df$bac <- c(rep("Em1021", 96), rep("Em1022", 96)) # Add in a bacteria label
df$bac <- as.factor(df$bac)

df$plate <- c(rep("1", 96), rep("2", 96)) # Add in plate label

df$ID <- paste(df$Well.ID, df$plate, sep=".") # Unique ID for each well

# OK so I have 14 files - 7 for EM1021 and 1022 each, the naming is as follows: em1021_pilot_t3.csv
# I am going to loop through these, transforming the data into a dataframe and then adding it to a master dataframe

################# Test run ########################

file.path <- "pilot-projects/em1022_pilot_t1.xlsx"
df.rawd <- read_excel(file.path, col_names = T, range = "B46:N78")

df.rawd <- read.xlsx(file.path, rows = 47:78, cols = 2:14, colNames = FALSE)

file.path <- "pilot-projects/em1022_pilot_t2.xlsx"
df.rawd <- read.xlsx(file.path, rows = 47:78, cols = 2:14, colNames = FALSE)

# The following is the code we need to run for every plate...

df.rawd$...1 <- zoo::na.locf(df.rawd$...1) # Fill in NAs with correct row ID
df.rawd$meas <- rep(c("RFU", "OD600", "OD700", "OD750"), 8)
df.rawd <- as.data.frame(df.rawd)

RFU <- numeric()
OD600 <- numeric()
OD700 <- numeric()
OD750 <- numeric()
Row <- character()
Col <- numeric()
plate <- character()

# Loop through the df.rawd

for (i in 1:nrow(df.rawd)){
  
  x <- ceiling(i/4) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
  
  for (c in 2:13){
    if (i == 4*x){
      Row <- c(Row, df.rawd$...1[i]) # Only record row and column every 4 rows.
      Col <- c(Col, c - 1) 
      if (grepl("1021", file.path)) {
        plate <- c(plate, "1")
      } else {
        plate <- c(plate, "2")
      }
    }# Numerical value for the column} 
    
    if(df.rawd$meas[i] == 'RFU'){
      RFU <- c(RFU, df.rawd[i,c])
    }
    
    if(df.rawd$meas[i]=='OD600'){
      OD600 <- c(OD600, df.rawd[i,c])
    }
    
    if(df.rawd$meas[i]=='OD700'){
      OD700 <- c(OD700, df.rawd[i,c])
    }
    
    if(df.rawd$meas[i]=='OD750'){
      OD750 <- c(OD750, df.rawd[i,c])
    }
    
  }
}

df.plate <- data.frame(
  Row = Row,
  Column = Col,
  plate = plate,
  RFU = RFU,
  OD600 = OD600,
  OD700 = OD700,
  OD750 = OD750
)

df.plate$Well.ID <-  paste(df.plate$Row, df.plate$Column, sep="")
df.plate$ID <- paste(df.plate$Well.ID, df.plate$plate, sep=".") # Unique ID for each well and plate

vars <- read_excel(file.path, range = "B7:B8", col_names = FALSE) # Get the date and time associated with the measurement
vars<-as.data.frame(vars)

date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
df.plate$date <- rep(julian(date, origin = as.Date("2024-01-01")), nrow(df.plate)) # Julian date

time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
df.plate$time <- rep(format(time, "%H:%M:%S"), nrow(df.plate)) # Extracted time

df.plate$datetime <- as.POSIXct(paste(date, df.plate$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

df.plate$elapsed_time <- c(0, diff(as.numeric(df.plate$datetime))) # Time in seconds

df1 <- df.plate %>%
  left_join(df, by = c("Row", "Column", "Well.ID","plate","ID")) # This is the full data we want for each frame (don't need elapsed time)

################# Proper loop ########################

df.pilot <- data.frame() # Empty data frame for adding the df.plt data. 

for (i in 0:6){ # File name indices
  
  file.path <- paste("pilot-projects/em1021_pilot_t", i, ".xlsx", sep = "") # 1021 data
  file.path2 <- paste("pilot-projects/em1022_pilot_t", i, ".xlsx", sep = "") # 1022 data
  
  df.rawd1 <- read_excel(file.path, range = "B46:N78")
  df.rawd2 <- read_excel(file.path2, range = "B46:N78")
  
  df.rawd1$...1 <- zoo::na.locf(df.rawd1$...1) # dfrawd1 first
  df.rawd1$meas <- rep(c("RFU", "OD600", "OD700", "OD750"), 8)
  df.rawd1 <- as.data.frame(df.rawd1)
  
  RFU <- numeric()
  OD600 <- numeric()
  OD700 <- numeric()
  OD750 <- numeric()
  Row <- character()
  Col <- numeric()
  plate <- character()
  
  for (i in 1:nrow(df.rawd1)){
    
    x <- ceiling(i/4) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
    
    for (c in 2:13){
      if (i == 4*x){
        Row <- c(Row, df.rawd1$...1[i]) # Only record row and column every 4 rows.
        Col <- c(Col, c - 1) 
        if (grepl("1021", file.path)) {
          plate <- c(plate, "1")
        } else {
          plate <- c(plate, "2")
        }
      }# Numerical value for the column} 
      
      if(df.rawd1$meas[i] == 'RFU'){
        RFU <- c(RFU, df.rawd1[i,c])
      }
      
      if(df.rawd1$meas[i]=='OD600'){
        OD600 <- c(OD600, df.rawd1[i,c])
      }
      
      if(df.rawd1$meas[i]=='OD700'){
        OD700 <- c(OD700, df.rawd1[i,c])
      }
      
      if(df.rawd1$meas[i]=='OD750'){
        OD750 <- c(OD750, df.rawd1[i,c])
      }
      
    }
  }
  
  df.plt <- data.frame( # Temporary holder
    Row = Row,
    Column = Col,
    plate = plate,
    RFU = RFU,
    OD600 = OD600,
    OD700 = OD700,
    OD750 = OD750
  ) 
  
  vars1 <- read_excel(file.path, range = "B7:B8", col_names = FALSE) # Get the date and time associated with the measurement
  vars1<-as.data.frame(vars1)
  
  date <- as.Date(vars1[1, 1], format = "%Y-%m-%d") #date
  df.plt$date <- rep(julian(date, origin = as.Date("2024-01-01")), nrow(df.plt)) # Julian date
  
  time <- as.POSIXct(vars1[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
  
  df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  df.1 <- df.plt %>%
    left_join(df, by = c("Row", "Column", "plate")) # This is the full data we want for each frame (don't need elapsed time)
  
  df.pilot <- rbind(df.pilot, df.1) 
  
  df.rawd2$...1 <- zoo::na.locf(df.rawd2$...1) # dfrawd2 now
  df.rawd2$meas <- rep(c("RFU", "OD600", "OD700", "OD750"), 8)
  df.rawd2 <- as.data.frame(df.rawd2)
  
  RFU <- numeric()
  OD600 <- numeric()
  OD700 <- numeric()
  OD750 <- numeric()
  Row <- character()
  Col <- numeric()
  plate <- character()
  
  for (i in 1:nrow(df.rawd2)){
    
    x <- ceiling(i/4) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
    
    for (c in 2:13){
      if (i == 4*x){
        Row <- c(Row, df.rawd2$...1[i]) # Only record row and column every 4 rows.
        Col <- c(Col, c - 1) 
        if (grepl("1021", file.path2)) {
          plate <- c(plate, "1")
        } else {
          plate <- c(plate, "2")
        }
      }# Numerical value for the column} 
      
      if(df.rawd2$meas[i] == 'RFU'){
        RFU <- c(RFU, df.rawd2[i,c])
      }
      
      if(df.rawd2$meas[i]=='OD600'){
        OD600 <- c(OD600, df.rawd2[i,c])
      }
      
      if(df.rawd2$meas[i]=='OD700'){
        OD700 <- c(OD700, df.rawd2[i,c])
      }
      
      if(df.rawd2$meas[i]=='OD750'){
        OD750 <- c(OD750, df.rawd2[i,c])
      }
      
    }
  }
  
  df.plt <- data.frame( # Temporary holder
    Row = Row,
    Column = Col,
    plate = plate,
    RFU = RFU,
    OD600 = OD600,
    OD700 = OD700,
    OD750 = OD750
  ) 
  
  vars2 <- read_excel(file.path2, range = "B7:B8", col_names = FALSE) # Get the date and time associated with the measurement
  vars2<-as.data.frame(vars2)
  
  date <- as.Date(vars2[1, 1], format = "%Y-%m-%d") #date
  df.plt$date <- rep(julian(date, origin = as.Date("2024-01-01")), nrow(df.plt)) # Julian date
  
  time <- as.POSIXct(vars2[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
  
  df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  df.2 <- df.plt %>%
    left_join(df, by = c("Row", "Column", "plate")) # This is the full data we want for each frame (don't need elapsed time)
  
  df.pilot <- rbind(df.pilot, df.2) 
}

df.pilot$days <- as.numeric(difftime(df.pilot$datetime, min(df.pilot$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

write.csv(df.pilot, "pilot-projects/05_Chlamy_Ensifer_pilot_data.csv") # Save cleaned data



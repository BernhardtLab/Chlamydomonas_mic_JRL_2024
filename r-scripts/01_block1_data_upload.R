# Jason R Laurich
# Apr 22, 2025


# Load packages -----------------------------------------------------------

library(tidyverse)
library(zoo)
library(readxl)
library(lubridate)
library(openxlsx)

# Script to upload and organize raw data --------------------------------

# 43 C data

t<-43

for (i in 0:x){ # File name indices
  
  file.path <- paste("pilot-projects/jl_1_", t, "C_", i, ".xlsx", sep = "")
  
  df.raw <- read_excel(file.path, range = "C58:K62")
  
  df.raw <- as.data.frame(df.raw)
  df.raw$meas <- c("RFU", "OD600", "OD700", "OD750")
  
  RFU <- numeric()
  OD600 <- numeric()
  OD700 <- numeric()
  OD750 <- numeric()
  Row <- character()
  Col <- numeric()
  Read <- character()
  
  for (j in 1:nrow(df.raw)){
    
    x <- ceiling(j/4) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
    
    for (c in 1:9){
      
      if (j == x){ # Only record row and column every four rows.
        Row <- c(Row, x) # row number
        Col <- c(Col, c) # column number  
        Read <- c(Read, i) # recorde the read #
      }
      
      if(df.raw$meas[j] == 'RFU'){
        RFU <- c(RFU, df.raw[j,c])
      }
      
      if(df.raw$meas[j]=='OD600'){
        OD600 <- c(OD600, df.raw[j,c])
      }
      
      if(df.raw$meas[j]=='OD700'){
        OD700 <- c(OD700, df.raw[j,c])
      }
      
      if(df.raw$meas[j]=='OD750'){
        OD750 <- c(OD750, df.raw[j,c])
      }
      
    }
    
  }  
  
  df.plt <- data.frame( # Temporary holder
    Row = Row,
    Column = Col,
    Read = Read,
    RFU = RFU,
    OD600 = OD600,
    OD700 = OD700,
    OD750 = OD750
  ) 
  
  vars <- read_excel(file.path, range = "B7:B8", col_names = FALSE) # Get the date and time associated with the measurement
  vars<-as.data.frame(vars)
  
  date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
  df.plt$date <- rep(julian(date, origin = as.Date("2025-01-01")), nrow(df.plt)) # Julian date
  
  time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
  
  df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  df.plt$Row <- as.integer(df.plt$Row)
  df.plt$Column <- as.integer(df.plt$Column)
  
  df.proc <- df.plt %>%
    left_join(df.design.temp, by = c("Row", "Column")) # This is the full data we want for each frame (don't need elapsed time)
  
  df.43 <- rbind(df.43, df.proc)
  
}  

df.43$days <- as.numeric(difftime(df.43$datetime, min(df.43$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.
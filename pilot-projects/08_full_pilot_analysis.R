# Jason R Laurich
# April 3, 2025

# Going to upload and process the raw data and pose 3 questions:

# (1) Can we back-calculate microbial OD contributions based on the relationship between RFUs and OD in Chlamy-only wells?
  # We will only use 30 C data due to evaporation at higher temperatures. 

# (2) Does bacterial inoculation affect growth rates of Chlamydomonas?
  # Here we will focus on growth rates, not TPCs per se, because we had so much evaporation at high temperatures

# (3) How do salt and nitrogen gradients affect Chlamydomonas and microbial growth?

############### Packages #############################

library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)
library(readxl)
library(lubridate)
library(openxlsx)

############### Access the csv files ######################

df.design <- read.csv("pilot-projects/09_full_pilot_summary.csv") # Get the experimental design file.

head(df.design) # This matches treatment to wells
str(df.design) 

# OK so I now have a bunch of raw files with the OD and RFU data from the plate reader that I need to read in and process
# I am going to loop through these, transforming the data into a dataframe and then adding it to a master dataframe

# The format is this: jl_2_30C_3, where 1/2 in jl_1/2_temp_read# is the plate number at a given temperature
# For the pilot, this was only for the 30 C incubator

# Need to keep in mind that the number of reads was different across treatments (colder took longer to grow)

############## Read the plate reader output files #########################

# Start with 8C data

df.8 <- data.frame() # Empty data frame for adding the df.plt data. 

t<-8
df.design.temp <- df.design[df.design$Temperature==t,]

for (i in c(1:4,6:8)){ # File name indices
  
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
  
  df.8 <- rbind(df.8, df.proc)
  
}  

df.8$days <- as.numeric(difftime(df.8$datetime, min(df.8$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# 14 C

df.14 <- data.frame() # Empty data frame for adding the df.plt data. 

t<-14
df.design.temp <- df.design[df.design$Temperature==t,]

for (i in c(1:4,6:8)){ # File name indices
  
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
  
  df.14 <- rbind(df.14, df.proc)
  
}  

df.14$days <- as.numeric(difftime(df.14$datetime, min(df.14$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# 20 C

df.20 <- data.frame() # Empty data frame for adding the df.plt data. 

t<-20
df.design.temp <- df.design[df.design$Temperature==t,]

for (i in c(1:4,6:8)){ # File name indices
  
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
  
  df.20 <- rbind(df.20, df.proc)
  
}  

df.20$days <- as.numeric(difftime(df.20$datetime, min(df.20$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.


# 25 C

df.25 <- data.frame() # Empty data frame for adding the df.plt data. 

t<-25
df.design.temp <- df.design[df.design$Temperature==t,]

for (i in c(1:4,6:8)){ # File name indices
  
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
  
  df.25 <- rbind(df.25, df.proc)
  
}  

df.25$days <- as.numeric(difftime(df.25$datetime, min(df.25$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# 33 C

df.33 <- data.frame() # Empty data frame for adding the df.plt data. 

t<-33
df.design.temp <- df.design[df.design$Temperature==t,]

for (i in c(1:4,6:8)){ # File name indices
  
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
  
  df.33 <- rbind(df.33, df.proc)
  
}  

df.33$days <- as.numeric(difftime(df.33$datetime, min(df.33$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# 35 C

df.35 <- data.frame() # Empty data frame for adding the df.plt data. 

t<-35
df.design.temp <- df.design[df.design$Temperature==t,]

for (i in c(1:4,6:8)){ # File name indices
  
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
  
  df.35 <- rbind(df.35, df.proc)
  
}  

df.35$days <- as.numeric(difftime(df.35$datetime, min(df.35$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# 37 C

df.37 <- data.frame() # Empty data frame for adding the df.plt data. 

t<-37
df.design.temp <- df.design[df.design$Temperature==t,]

for (i in c(1:4,6:8)){ # File name indices
  
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
  
  df.37 <- rbind(df.37, df.proc)
  
}  

df.37$days <- as.numeric(difftime(df.37$datetime, min(df.37$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# 40 C

df.40 <- data.frame() # Empty data frame for adding the df.plt data. 

t<-40
df.design.temp <- df.design[df.design$Temperature==t,]

for (i in c(1:4,6:8)){ # File name indices
  
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
  
  df.40 <- rbind(df.40, df.proc)
  
}  

df.40$days <- as.numeric(difftime(df.40$datetime, min(df.40$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# 43 C

df.43 <- data.frame() # Empty data frame for adding the df.plt data. 

t<-43
df.design.temp <- df.design[df.design$Temperature==t,]

for (i in c(1:4,6:8)){ # File name indices
  
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

# Finally the 30C data. The issue here is that we have multiple plates that don't have the same number of reads. But they 
# will in the future!

df.30 <- data.frame() # Empty data frame for adding the df.plt data. 

t<-30
df.design.temp <- df.design[df.design$Temperature==t,]
df.design.temp <- df.design.temp[1:96,]

for (i in c(1:4,6:8)){ # File name indices
  
  file.path1 <- paste("pilot-projects/jl_1_", t, "C_", i, ".xlsx", sep = "")
  # file.path2 <- paste("pilot-projects/jl_2_", t, "C_", i, ".xlsx", sep = "") # For now we are going to ignore this - there were consistency issues with the plate reading
  
  df.raw <- read_excel(file.path1, range = "C58:N90")
  
  df.raw <- as.data.frame(df.raw)
  df.raw$meas <- rep(c("RFU", "OD600", "OD700", "OD750"), 8)
  
  #df.raw2 <- read_excel(file.path2, range = "C58:N90")
  
  #df.raw2 <- as.data.frame(df.raw2)
  #df.raw2$meas <- rep(c("RFU", "OD600", "OD700", "OD750"), 8)
  
  RFU <- numeric()
  OD600 <- numeric()
  OD700 <- numeric()
  OD750 <- numeric()
  Row <- character()
  Col <- numeric()
  Read <- character()
  
  for (j in 1:nrow(df.raw)){
    
    x <- ceiling(j/4) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
    
    for (c in 1:12){
      
      if (j == x){ # Only record row and column every four rows.
        Row <- c(Row, x) # row number
        Col <- c(Col, c) # column number  
        Read <- c(Read, i) # record the read #
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
  
  vars <- read_excel(file.path1, range = "B7:B8", col_names = FALSE) # Get the date and time associated with the measurement
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
  
  df.30 <- rbind(df.30, df.proc)
  
}  

df.30$days <- as.numeric(difftime(df.30$datetime, min(df.30$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# Let's combine these all!

df <- rbind(df.8, df.14, df.20, df.25, df.30, df.33, df.35, df.37, df.40, df.43)

write.csv(df, "pilot-projects/10_full_pilot_data.csv") # Save cleaned data

# Jason R Laurich
# July 2, 2025

# Upload and save the data from the second block of the experiment

# Load packages -----------------------------------------------------------

library(tidyverse)
library(zoo)
library(readxl)
library(lubridate)
library(openxlsx)
library(cowplot)
library(dplyr)

# Script to upload and organize raw data --------------------------------

df.design <- read.csv("raw-data/04-block-2-design.csv") # Get the experimental design file.

head(df.design) # This matches treatment to wells
str(df.design)  # Design sheet

###### 43 C data ######

t<-43 # Temp

df.43 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:11){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "D50:M62")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),6)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:10){
        
        if (j == 2*x){ # Only record row and column every two rows.
          Row <- c(Row, x) # row number
          Col <- c(Col, c) # column number  
          Plate <- c(Plate, p) # record the plate #
          Read <- c(Read, i) # record the read #
        }
        
        if(df.raw$meas[j] == 'RFU'){
          RFU <- c(RFU, df.raw[j,c])
        }
        
        if(df.raw$meas[j]=='OD600'){
          OD600 <- c(OD600, df.raw[j,c])
        }
        
      }
      
    }
    
    df.plt <- data.frame( # Temporary holder
      Row = Row,
      Column = Col,
      Plate = Plate,
      Read = Read,
      RFU = RFU,
      OD600 = OD600
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
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.43 <- rbind(df.43, df.proc)
    
  }  
  
}

df.43$days <- as.numeric(difftime(df.43$datetime, min(df.43$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

p43 <- ggplot(df.43, aes(x = days, y = RFU, group = Well.at.T)) +
  geom_line(alpha = 0.3) +
  theme_classic() +
  labs(title="43 C")

p43

###### 39 C ######

t<-39 # Temp

df.39 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:11){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "D50:M62")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),6)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:10){
        
        if (j == 2*x){ # Only record row and column every two rows.
          Row <- c(Row, x) # row number
          Col <- c(Col, c) # column number  
          Plate <- c(Plate, p) # record the plate #
          Read <- c(Read, i) # record the read #
        }
        
        if(df.raw$meas[j] == 'RFU'){
          RFU <- c(RFU, df.raw[j,c])
        }
        
        if(df.raw$meas[j]=='OD600'){
          OD600 <- c(OD600, df.raw[j,c])
        }
        
      }
      
    }
    
    df.plt <- data.frame( # Temporary holder
      Row = Row,
      Column = Col,
      Plate = Plate,
      Read = Read,
      RFU = RFU,
      OD600 = OD600
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
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.39 <- rbind(df.39, df.proc)
    
  }  
  
}

df.39$days <- as.numeric(difftime(df.39$datetime, min(df.39$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

p39 <- ggplot(df.39, aes(x = days, y = RFU, group = Well.at.T)) +
  geom_line(alpha = 0.3) +
  theme_classic() +
  labs(title="39 C")

p39

###### 35 C ######

t<-35 # Temp

df.35 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:11){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "D50:M62")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),6)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:10){
        
        if (j == 2*x){ # Only record row and column every two rows.
          Row <- c(Row, x) # row number
          Col <- c(Col, c) # column number  
          Plate <- c(Plate, p) # record the plate #
          Read <- c(Read, i) # record the read #
        }
        
        if(df.raw$meas[j] == 'RFU'){
          RFU <- c(RFU, df.raw[j,c])
        }
        
        if(df.raw$meas[j]=='OD600'){
          OD600 <- c(OD600, df.raw[j,c])
        }
        
      }
      
    }
    
    df.plt <- data.frame( # Temporary holder
      Row = Row,
      Column = Col,
      Plate = Plate,
      Read = Read,
      RFU = RFU,
      OD600 = OD600
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
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.35 <- rbind(df.35, df.proc)
    
  }  
  
}

df.35$days <- as.numeric(difftime(df.35$datetime, min(df.35$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

p35 <- ggplot(df.35, aes(x = days, y = RFU, group = Well.at.T)) +
  geom_line(alpha = 0.3) +
  theme_classic() +
  labs(title="35 C")

p35

###### 33 C ######

t<-33 # Temp

df.33 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:11){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "D50:M62")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),6)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:10){
        
        if (j == 2*x){ # Only record row and column every two rows.
          Row <- c(Row, x) # row number
          Col <- c(Col, c) # column number  
          Plate <- c(Plate, p) # record the plate #
          Read <- c(Read, i) # record the read #
        }
        
        if(df.raw$meas[j] == 'RFU'){
          RFU <- c(RFU, df.raw[j,c])
        }
        
        if(df.raw$meas[j]=='OD600'){
          OD600 <- c(OD600, df.raw[j,c])
        }
        
      }
      
    }
    
    df.plt <- data.frame( # Temporary holder
      Row = Row,
      Column = Col,
      Plate = Plate,
      Read = Read,
      RFU = RFU,
      OD600 = OD600
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
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.33 <- rbind(df.33, df.proc)
    
  }  
  
}

df.33$days <- as.numeric(difftime(df.33$datetime, min(df.33$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

p33 <- ggplot(df.33, aes(x = days, y = RFU, group = Well.at.T)) +
  geom_line(alpha = 0.3) +
  theme_classic() +
  labs(title="33 C")

p33

###### 25 C ######

t<-25 # Temp

df.25 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "D50:M62")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),6)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:10){
        
        if (j == 2*x){ # Only record row and column every two rows.
          Row <- c(Row, x) # row number
          Col <- c(Col, c) # column number  
          Plate <- c(Plate, p) # record the plate #
          Read <- c(Read, i) # record the read #
        }
        
        if(df.raw$meas[j] == 'RFU'){
          RFU <- c(RFU, df.raw[j,c])
        }
        
        if(df.raw$meas[j]=='OD600'){
          OD600 <- c(OD600, df.raw[j,c])
        }
        
      }
      
    }
    
    df.plt <- data.frame( # Temporary holder
      Row = Row,
      Column = Col,
      Plate = Plate,
      Read = Read,
      RFU = RFU,
      OD600 = OD600
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
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.25 <- rbind(df.25, df.proc)
    
  }  
  
}

df.25$days <- as.numeric(difftime(df.25$datetime, min(df.25$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

p25 <- ggplot(df.25, aes(x = days, y = RFU, group = Well.at.T)) +
  geom_line(alpha = 0.3) +
  theme_classic() +
  labs(title="25 C")

p25

###### 20 C ######

t<-20 # Temp

df.20 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:16){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "D50:M62")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),6)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:10){
        
        if (j == 2*x){ # Only record row and column every two rows.
          Row <- c(Row, x) # row number
          Col <- c(Col, c) # column number  
          Plate <- c(Plate, p) # record the plate #
          Read <- c(Read, i) # record the read #
        }
        
        if(df.raw$meas[j] == 'RFU'){
          RFU <- c(RFU, df.raw[j,c])
        }
        
        if(df.raw$meas[j]=='OD600'){
          OD600 <- c(OD600, df.raw[j,c])
        }
        
      }
      
    }
    
    df.plt <- data.frame( # Temporary holder
      Row = Row,
      Column = Col,
      Plate = Plate,
      Read = Read,
      RFU = RFU,
      OD600 = OD600
    ) 
    
    vars <- read_excel(file.path, range = "B7:B8", col_names = FALSE) # Get the date and time associated with the measurement
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2020-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.20 <- rbind(df.20, df.proc)
    
  }  
  
}

df.20$days <- as.numeric(difftime(df.20$datetime, min(df.20$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

p20 <- ggplot(df.20, aes(x = days, y = RFU, group = Well.at.T)) +
  geom_line(alpha = 0.3) +
  theme_classic() +
  labs(title="20 C")

p20

###### 14 C ######

t<-14 # Temp

df.14 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:17){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "D50:M62")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),6)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:10){
        
        if (j == 2*x){ # Only record row and column every two rows.
          Row <- c(Row, x) # row number
          Col <- c(Col, c) # column number  
          Plate <- c(Plate, p) # record the plate #
          Read <- c(Read, i) # record the read #
        }
        
        if(df.raw$meas[j] == 'RFU'){
          RFU <- c(RFU, df.raw[j,c])
        }
        
        if(df.raw$meas[j]=='OD600'){
          OD600 <- c(OD600, df.raw[j,c])
        }
        
      }
      
    }
    
    df.plt <- data.frame( # Temporary holder
      Row = Row,
      Column = Col,
      Plate = Plate,
      Read = Read,
      RFU = RFU,
      OD600 = OD600
    ) 
    
    vars <- read_excel(file.path, range = "B7:B8", col_names = FALSE) # Get the date and time associated with the measurement
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2014-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.14 <- rbind(df.14, df.proc)
    
  }  
  
}

df.14$days <- as.numeric(difftime(df.14$datetime, min(df.14$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

p14 <- ggplot(df.14, aes(x = days, y = RFU, group = Well.at.T)) +
  geom_line(alpha = 0.3) +
  theme_classic() +
  labs(title="14 C")

p14

###### 8 C ######

t<-8 # Temp

df.8 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:17){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "D50:M62")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),6)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:10){
        
        if (j == 2*x){ # Only record row and column every two rows.
          Row <- c(Row, x) # row number
          Col <- c(Col, c) # column number  
          Plate <- c(Plate, p) # record the plate #
          Read <- c(Read, i) # record the read #
        }
        
        if(df.raw$meas[j] == 'RFU'){
          RFU <- c(RFU, df.raw[j,c])
        }
        
        if(df.raw$meas[j]=='OD600'){
          OD600 <- c(OD600, df.raw[j,c])
        }
        
      }
      
    }
    
    df.plt <- data.frame( # Temporary holder
      Row = Row,
      Column = Col,
      Plate = Plate,
      Read = Read,
      RFU = RFU,
      OD600 = OD600
    ) 
    
    vars <- read_excel(file.path, range = "B7:B8", col_names = FALSE) # Get the date and time associated with the measurement
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("208-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.8 <- rbind(df.8, df.proc)
    
  }  
  
}

df.8$days <- as.numeric(difftime(df.8$datetime, min(df.8$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

p8 <- ggplot(df.8, aes(x = days, y = RFU, group = Well.at.T)) +
  geom_line(alpha = 0.3) +
  theme_classic() +
  labs(title="8 C")

p8

###### 30 C ######

t<-30 # Temp

df.30 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:32){ # Plates
  
  df.design.temp.plate <- df.design %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:11) { # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    
    tryCatch({
      
      df.raw <- read_excel(file.path, range = "D50:M62") %>% as.data.frame()
      
      df.raw$meas <- rep(c("RFU", "OD600"),6)
      
      RFU <- numeric()
      OD600 <- numeric()
      Row <- numeric()
      Col <- numeric()
      Plate <- numeric()
      Read <- numeric()
      
      for (j in 1:nrow(df.raw)){
        x <- ceiling(j/2)
        for (c in 1:10){
          if (j == 2*x){
            Row <- c(Row, x)
            Col <- c(Col, c)
            Plate <- c(Plate, p)
            Read <- c(Read, i)
          }
          if(df.raw$meas[j] == 'RFU'){
            RFU <- c(RFU, df.raw[j,c])
          }
          if(df.raw$meas[j]=='OD600'){
            OD600 <- c(OD600, df.raw[j,c])
          }
        }
      }
      
      df.plt <- data.frame(
        Row = Row,
        Column = Col,
        Plate = Plate,
        Read = Read,
        RFU = RFU,
        OD600 = OD600
      )
      
      vars <- read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame()
      date <- as.Date(vars[1, 1], format = "%Y-%m-%d")
      df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt))
      time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
      df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt))
      df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
      df.plt$Row <- as.integer(df.plt$Row)
      df.plt$Column <- as.integer(df.plt$Column)
      
      df.plt <- df.plt %>% 
        mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
      
      df.proc <- df.plt %>%
        left_join(df.design.temp.plate, by = "Well.at.T")
      
      df.30 <- rbind(df.30, df.proc)
      
    }, error = function(e){
      message(paste("File missing:", file.path, " - skipping"))
      # Skip to next i
    })
    
  }
  
}

df.30$days <- as.numeric(difftime(df.30$datetime, min(df.30$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

p30 <- ggplot(df.30, aes(x = days, y = RFU, group = Well.at.T)) +
  geom_line(alpha = 0.3) +
  theme_classic() +
  labs(title="30 C")

p30

# Data exploration with plots ---------------------------------------------

p.comp <- plot_grid(p8, p14, p20, p25, p30, p33, p35, p39, p43,
                    ncol = 5)

p.comp

df.blk2 <- rbind(df.8, df.14, df.20, df.25, df.30, df.33, df.35, df.39, df.43)

write.csv(df.blk2, "processed-data/03_blk2_rawdata.csv") # Save raw data. This still needs some accounting (e.g. errors in the notes column)

head(df.blk2)
str(df.blk2)

df.blk2$Microbe <- as.factor(df.blk2$Microbe)
levels(df.blk2$Microbe) # Right so we have BLANK in there as an extra factor level. That's fine.

df.blk2 <- df.blk2 %>%
  filter(Microbe != 'BLANK') %>%
  droplevels() %>%
  mutate(Microbe = fct_relevel(Microbe, 
                               "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                               "9", "10", "11", "12", "13", "14", "15", "all"))

df.blk2.C <- df.blk2 %>%
  filter(Chlamy.y.n == 'y')

p8.2 <- ggplot(df.blk2.C[df.blk2.C$Temperature.C == 8, ], 
               aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "8 C") +
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
  )

p8.2

p14.2 <- ggplot(df.blk2.C[df.blk2.C$Temperature.C == 14, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "14 C") +
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
  )

p14.2

p20.2 <- ggplot(df.blk2.C[df.blk2.C$Temperature.C == 20, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "20 C") +
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
  )

p20.2

p25.2 <- ggplot(df.blk2.C[df.blk2.C$Temperature.C == 25, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "25 C") +
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
  )

p25.2

p33.2 <- ggplot(df.blk2.C[df.blk2.C$Temperature.C == 33, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "33 C") +
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
  )

p33.2

p35.2 <- ggplot(df.blk2.C[df.blk2.C$Temperature.C == 35, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "35 C") +
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
  )

p35.2

p39.2 <- ggplot(df.blk2.C[df.blk2.C$Temperature.C == 39, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "39 C") +
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
  )

p39.2

p43.2 <- ggplot(df.blk2.C[df.blk2.C$Temperature.C == 43, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "43 C") +
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
  )

p43.2

df.blk2.C.t <- df.blk2.C %>% 
  filter(Nitrogen.conc.µM == 1000, Salt.conc.g.l == 0)

p30.2 <- ggplot(df.blk2.C.t[df.blk2.C.t$Temperature.C == 30, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "30 C") +
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
  )

p30.2

p.comp.2 <- plot_grid(p8.2, p14.2, p20.2, p25.2, p30.2, p33.2, p35.2, p39.2, p43.2,
                      ncol = 3)

p.comp.2

###### Nitrogen ######

df.blk2.n.C <- df.blk2 %>% 
  filter(Salt.conc.g.l == 0, Temperature.C == 30, Chlamy.y.n =='y')

p.n.0 <- ggplot(df.blk2.n.C[df.blk2.n.C$Nitrogen.conc.µM == 0 & df.blk2.n.C$RFU < 250, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "0 µM N") +
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
  )

p.n.0

p.n.1 <- ggplot(df.blk2.n.C[df.blk2.n.C$Nitrogen.conc.µM == 1 & df.blk2.n.C$RFU < 150,], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "1 µM N") +
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
  )

p.n.1

p.n.5 <- ggplot(df.blk2.n.C[df.blk2.n.C$Nitrogen.conc.µM == 5 & df.blk2.n.C$RFU <300, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "5 µM N") +
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
  )

p.n.5

p.n.10 <- ggplot(df.blk2.n.C[df.blk2.n.C$Nitrogen.conc.µM == 10 & df.blk2.n.C$RFU < 1000, ], 
                 aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "10 µM N") +
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
  )

p.n.10

p.n.25 <- ggplot(df.blk2.n.C[df.blk2.n.C$Nitrogen.conc.µM == 25, ], 
                 aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "25 µM N") +
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
  )

p.n.25

p.n.50 <- ggplot(df.blk2.n.C[df.blk2.n.C$Nitrogen.conc.µM == 50, ], 
                 aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "50 µM N") +
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
  )

p.n.50

p.n.100 <- ggplot(df.blk2.n.C[df.blk2.n.C$Nitrogen.conc.µM == 100, ], 
                  aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "100 µM N") +
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
  )

p.n.100

p.n.400 <- ggplot(df.blk2.n.C[df.blk2.n.C$Nitrogen.conc.µM == 400, ], 
                  aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "400 µM N") +
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
  )

p.n.400

p.n.700 <- ggplot(df.blk2.n.C[df.blk2.n.C$Nitrogen.conc.µM == 700, ], 
                  aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "700 µM N") +
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
  )

p.n.700

p.n.1000 <- ggplot(df.blk2.n.C[df.blk2.n.C$Nitrogen.conc.µM == 1000, ], 
                   aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "1000 µM N") +
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
  )

p.n.1000

plot_grid(p.n.0, p.n.1, p.n.5, p.n.10, p.n.25, p.n.50, p.n.100, p.n.400, p.n.700, p.n.1000, ncol = 4)

###### Salt ######

df.blk2.s.C <- df.blk2 %>% 
  filter(Nitrogen.conc.µM == 1000, Temperature.C == 30, Chlamy.y.n =='y')

p.s.0 <- ggplot(df.blk2.n.C[df.blk2.n.C$Salt.conc.g.l == 0, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "0 g/L Salt") +
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
  )

p.s.0

p.s.1.33 <- ggplot(df.blk2.s.C[df.blk2.s.C$Salt.conc.g.l == 1.33, ], 
                   aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "1.33 g/L Salt") +
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
  )

p.s.1.33

p.s.2.67 <- ggplot(df.blk2.s.C[df.blk2.s.C$Salt.conc.g.l == 2.67, ], 
                   aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "2.67 g/L Salt") +
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
  )

p.s.2.67

p.s.4 <- ggplot(df.blk2.s.C[df.blk2.s.C$Salt.conc.g.l == 4, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "4 g/L Salt") +
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
  )

p.s.4

p.s.5.33 <- ggplot(df.blk2.s.C[df.blk2.s.C$Salt.conc.g.l == 5.33 & df.blk2.s.C$RFU < 200, ], 
                   aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "5.33 g/L Salt") +
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
  )

p.s.5.33

p.s.6.67 <- ggplot(df.blk2.s.C[df.blk2.s.C$Salt.conc.g.l == 6.67 & df.blk2.s.C$RFU < 120, ], 
                   aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "6.67 g/L Salt") +
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
  )

p.s.6.67

p.s.8 <- ggplot(df.blk2.s.C[df.blk2.s.C$Salt.conc.g.l == 8 & df.blk2.s.C$RFU < 150, ], 
                aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "8 g/L Salt") +
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
  )

p.s.8

p.s.9.33 <- ggplot(df.blk2.s.C[df.blk2.s.C$Salt.conc.g.l == 9.33 & df.blk2.s.C$RFU < 65, ], 
                   aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "9.33 g/L Salt") +
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
  )

p.s.9.33

p.s.10.67 <- ggplot(df.blk2.s.C[df.blk2.s.C$Salt.conc.g.l == 10.67 & df.blk2.s.C$RFU < 100, ], 
                    aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "10.67 g/L Salt") +
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
  )

p.s.10.67

p.s.12 <- ggplot(df.blk2.s.C[df.blk2.s.C$Salt.conc.g.l == 12 & df.blk2.s.C$RFU < 52, ], 
                 aes(x = days, y = RFU, group = Well.at.T, colour = Microbe)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "12 g/L Salt") +
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
  )

p.s.12

plot_grid(p.s.0, p.s.1.33, p.s.2.67, p.s.4, p.s.5.33, p.s.6.67, p.s.8, p.s.9.33, p.s.10.67, p.s.12, ncol = 4)


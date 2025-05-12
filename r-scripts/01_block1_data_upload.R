# Jason R Laurich
# Apr 22, 2025

# Upload and save the data from the first block of the experiment
# Will need to add in the later reads when I'm on campus (wifi)

# Load packages -----------------------------------------------------------

library(tidyverse)
library(zoo)
library(readxl)
library(lubridate)
library(openxlsx)
library(cowplot)
library(dplyr)

# Script to upload and organize raw data --------------------------------

df.design <- read.csv("raw-data/02-block-1-design.csv") # Get the experimental design file.

head(df.design) # This matches treatment to wells
str(df.design)  # Design sheet

###### 43 C data ######

t<-43 # Temp

df.43 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:10){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
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
  
  for (i in 0:10){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "C48:N64")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),8)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:12){
        
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
      mutate(Well.at.T = (Plate-1)*96 + (Row-1)*12 + Column)
    
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
  
  for (i in 0:10){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "C48:N64")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),8)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:12){
        
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
      mutate(Well.at.T = (Plate-1)*96 + (Row-1)*12 + Column)
    
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
  
  for (i in 0:10){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "C48:N64")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),8)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:12){
        
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
      mutate(Well.at.T = (Plate-1)*96 + (Row-1)*12 + Column)
    
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
  
  for (i in 0:10){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "C48:N64")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),8)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:12){
        
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
      mutate(Well.at.T = (Plate-1)*96 + (Row-1)*12 + Column)
    
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
  
  for (i in 0:10){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "C48:N64")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),8)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:12){
        
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
      mutate(Well.at.T = (Plate-1)*96 + (Row-1)*12 + Column)
    
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
  
  for (i in 0:10){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "C48:N64")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),8)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:12){
        
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
      mutate(Well.at.T = (Plate-1)*96 + (Row-1)*12 + Column)
    
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
  
  for (i in 0:10){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "C48:N64")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),8)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:12){
        
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
      mutate(Well.at.T = (Plate-1)*96 + (Row-1)*12 + Column)
    
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
  
  for (i in 0:7){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    df.raw <- read_excel(file.path, range = "C48:N64")
    
    df.raw <- as.data.frame(df.raw)
    df.raw$meas <- rep(c("RFU", "OD600"),8)
    
    RFU <- numeric()
    OD600 <- numeric()
    
    Row <- numeric()
    Col <- numeric()
    Plate <- numeric()
    Read <- numeric()
    
    for (j in 1:nrow(df.raw)){
      
      x <- ceiling(j/2) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
      
      for (c in 1:12){
        
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
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*96 + (Row-1)*12 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.30 <- rbind(df.30, df.proc)
    
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

df.blk1 <- rbind(df.8, df.14, df.20, df.25, df.30, df.33, df.35, df.39, df.43)

write.csv(df.blk1, "processed-data/01_blk1_rawdata.csv") # Save raw data. This still needs some accounting (e.g. errors in the notes column)

head(df.blk1)
str(df.blk1)

df.blk1$Microbe <- as.factor(df.blk1$Microbe)
levels(df.blk1$Microbe) # Right so we have BLANK in there as an extra factor level. That's fine.

df.blk1 <- df.blk1 %>%
  filter(Microbe != 'BLANK') %>%
  droplevels() %>%
  mutate(Microbe = fct_relevel(Microbe, 
                               "none", "1", "2", "3", "4", "5", "6", "7", "8", 
                               "9", "10", "11", "12", "13", "14", "15", "all"))

p8.2 <- ggplot(df.blk1[df.blk1$Temperature.C == 8, ], 
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

p14.2 <- ggplot(df.blk1[df.blk1$Temperature.C == 14, ], 
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

p20.2 <- ggplot(df.blk1[df.blk1$Temperature.C == 20, ], 
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

p25.2 <- ggplot(df.blk1[df.blk1$Temperature.C == 25, ], 
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

p33.2 <- ggplot(df.blk1[df.blk1$Temperature.C == 33, ], 
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

p35.2 <- ggplot(df.blk1[df.blk1$Temperature.C == 35, ], 
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

p39.2 <- ggplot(df.blk1[df.blk1$Temperature.C == 39, ], 
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

p43.2 <- ggplot(df.blk1[df.blk1$Temperature.C == 43, ], 
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

df.blk1.t <- df.blk1 %>% 
  filter(Nitrogen.conc.µM == 1000, Salt.conc.g.l == 0)

p30.2 <- ggplot(df.blk1.t[df.blk1.t$Temperature.C == 30, ], 
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

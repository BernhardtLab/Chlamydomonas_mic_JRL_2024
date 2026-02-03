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

# We will extract data for everything and remove/deal with these at the µ estimation phase

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
    
df.design.4 <- read.csv("raw-data/08-block-4-design.csv") # Get the experimental design file.
head(df.design.4) # This matches treatment to wells
df.design.4$Block <- 4

df.design.4 <- df.design.4 %>% 
  mutate(unique.id = paste0("b4.t", Temperature.C, ".p", Plate.at.T, ".w", Well.at.T))

df.design.4 %>% 
  filter(!is.na(Notes), trimws(Notes) != "") # Let's identify the problematic rows — e.g. where I left notes on potential errors etc. 

# Block 1 -----------------------------------------------------------------

###### 8 C ######

t<-8 # Temp

df.1.8 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.1 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.1.8 <- rbind(df.1.8, df.proc)
    
  }  
  
}

df.1.8$days <- as.numeric(difftime(df.1.8$datetime, min(df.1.8$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 14 C ######

t<-14 # Temp

df.1.14 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.1 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.1.14 <- rbind(df.1.14, df.proc)
    
  }  
  
}

df.1.14$days <- as.numeric(difftime(df.1.14$datetime, min(df.1.14$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 20 C ######

t<-20 # Temp

df.1.20 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.1 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.1.20 <- rbind(df.1.20, df.proc)
    
  }  
  
}

df.1.20$days <- as.numeric(difftime(df.1.20$datetime, min(df.1.20$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 25 C ######

t<-25 # Temp

df.1.25 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.1 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.1.25 <- rbind(df.1.25, df.proc)
    
  }  
  
}

df.1.25$days <- as.numeric(difftime(df.1.25$datetime, min(df.1.25$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 30 C ######

t<-30 # Temp

df.1.30 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:32){ # Plates
  
  df.design.temp.plate <- df.design.1 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.1.30 <- rbind(df.1.30, df.proc)
    
  }  
  
}

df.1.30$days <- as.numeric(difftime(df.1.30$datetime, min(df.1.30$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 33 C ######

t<-33 # Temp

df.1.33 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.1 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.1.33 <- rbind(df.1.33, df.proc)
    
  }  
  
}

df.1.33$days <- as.numeric(difftime(df.1.33$datetime, min(df.1.33$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 35 C ######

t<-35 # Temp

df.1.35 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.1 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.1.35 <- rbind(df.1.35, df.proc)
    
  }  
  
}

df.1.35$days <- as.numeric(difftime(df.1.35$datetime, min(df.1.35$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 39 C ######

t<-39 # Temp

df.1.39 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.1 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.1.39 <- rbind(df.1.39, df.proc)
    
  }  
  
}

df.1.39$days <- as.numeric(difftime(df.1.39$datetime, min(df.1.39$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 43 C ######

t<-43 # Temp

df.1.43 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.1 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block1_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.1.43 <- rbind(df.1.43, df.proc)
    
  }  
  
}

df.1.43$days <- as.numeric(difftime(df.1.43$datetime, min(df.1.43$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# Block 2 -----------------------------------------------------------------

###### 8 C ######

t<-8 # Temp

df.2.8 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.2 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.2.8 <- rbind(df.2.8, df.proc)
    
  }  
  
}

df.2.8$days <- as.numeric(difftime(df.2.8$datetime, min(df.2.8$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 14 C ######

t<-14 # Temp

df.2.14 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.2 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.2.14 <- rbind(df.2.14, df.proc)
    
  }  
  
}

df.2.14$days <- as.numeric(difftime(df.2.14$datetime, min(df.2.14$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 20 C ######

t<-20 # Temp

df.2.20 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.2 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.2.20 <- rbind(df.2.20, df.proc)
    
  }  
  
}

df.2.20$days <- as.numeric(difftime(df.2.20$datetime, min(df.2.20$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 25 C ######

t<-25 # Temp

df.2.25 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.2 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.2.25 <- rbind(df.2.25, df.proc)
    
  }  
  
}

df.2.25$days <- as.numeric(difftime(df.2.25$datetime, min(df.2.25$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 30 C ######

t<-30 # Temp

df.2.30 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:32){ # Plates
  
  df.design.temp.plate <- df.design.2 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.2.30 <- rbind(df.2.30, df.proc)
    
  }  
  
}

df.2.30$days <- as.numeric(difftime(df.2.30$datetime, min(df.2.30$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 33 C ######

t<-33 # Temp

df.2.33 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.2 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.2.33 <- rbind(df.2.33, df.proc)
    
  }  
  
}

df.2.33$days <- as.numeric(difftime(df.2.33$datetime, min(df.2.33$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 35 C ######

t<-35 # Temp

df.2.35 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.2 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.2.35 <- rbind(df.2.35, df.proc)
    
  }  
  
}

df.2.35$days <- as.numeric(difftime(df.2.35$datetime, min(df.2.35$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 39 C ######

t<-39 # Temp

df.2.39 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.2 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.2.39 <- rbind(df.2.39, df.proc)
    
  }  
  
}

df.2.39$days <- as.numeric(difftime(df.2.39$datetime, min(df.2.39$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 43 C ######

t<-43 # Temp

df.2.43 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.2 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block2_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.2.43 <- rbind(df.2.43, df.proc)
    
  }  
  
}

df.2.43$days <- as.numeric(difftime(df.2.43$datetime, min(df.2.43$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# Block 3 -----------------------------------------------------------------

###### 8 C ######

t<-8 # Temp

df.3.8 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.3 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block3_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.3.8 <- rbind(df.3.8, df.proc)
    
  }  
  
}

df.3.8$days <- as.numeric(difftime(df.3.8$datetime, min(df.3.8$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 14 C ######

t<-14 # Temp

df.3.14 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.3 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block3_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.3.14 <- rbind(df.3.14, df.proc)
    
  }  
  
}

df.3.14$days <- as.numeric(difftime(df.3.14$datetime, min(df.3.14$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 20 C ######

t<-20 # Temp

df.3.20 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.3 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block3_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.3.20 <- rbind(df.3.20, df.proc)
    
  }  
  
}

df.3.20$days <- as.numeric(difftime(df.3.20$datetime, min(df.3.20$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 25 C ######

t<-25 # Temp

df.3.25 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.3 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block3_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.3.25 <- rbind(df.3.25, df.proc)
    
  }  
  
}

df.3.25$days <- as.numeric(difftime(df.3.25$datetime, min(df.3.25$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 30 C ######

t<-30 # Temp

df.3.30 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:32){ # Plates
  
  df.design.temp.plate <- df.design.3 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block3_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.3.30 <- rbind(df.3.30, df.proc)
    
  }  
  
}

df.3.30$days <- as.numeric(difftime(df.3.30$datetime, min(df.3.30$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 33 C ######

t<-33 # Temp

df.3.33 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.3 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:12){ # Reads
    
    file.path <- paste("raw-data/JRL_block3_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.3.33 <- rbind(df.3.33, df.proc)
    
  }  
  
}

df.3.33$days <- as.numeric(difftime(df.3.33$datetime, min(df.3.33$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 35 C ######

t<-35 # Temp

df.3.35 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.3 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:12){ # Reads
    
    file.path <- paste("raw-data/JRL_block3_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.3.35 <- rbind(df.3.35, df.proc)
    
  }  
  
}

df.3.35$days <- as.numeric(difftime(df.3.35$datetime, min(df.3.35$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 39 C ######

t<-39 # Temp

df.3.39 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.3 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:12){ # Reads
    
    file.path <- paste("raw-data/JRL_block3_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.3.39 <- rbind(df.3.39, df.proc)
    
  }  
  
}

df.3.39$days <- as.numeric(difftime(df.3.39$datetime, min(df.3.39$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 43 C ######

t<-43 # Temp

df.3.43 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.3 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:12){ # Reads
    
    file.path <- paste("raw-data/JRL_block3_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.3.43 <- rbind(df.3.43, df.proc)
    
  }  
  
}

df.3.43$days <- as.numeric(difftime(df.3.43$datetime, min(df.3.43$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

# Block 4 -----------------------------------------------------------------

###### 8 C ######

t<-8 # Temp

df.4.8 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.4 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block4_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.4.8 <- rbind(df.4.8, df.proc)
    
  }  
  
}

df.4.8$days <- as.numeric(difftime(df.4.8$datetime, min(df.4.8$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 14 C ######

t<-14 # Temp

df.4.14 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.4 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block4_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.4.14 <- rbind(df.4.14, df.proc)
    
  }  
  
}

df.4.14$days <- as.numeric(difftime(df.4.14$datetime, min(df.4.14$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 20 C ######

t<-20 # Temp

df.4.20 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.4 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block4_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.4.20 <- rbind(df.4.20, df.proc)
    
  }  
  
}

df.4.20$days <- as.numeric(difftime(df.4.20$datetime, min(df.4.20$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 25 C ######

t<-25 # Temp

df.4.25 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.4 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:20){ # Reads
    
    file.path <- paste("raw-data/JRL_block4_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.4.25 <- rbind(df.4.25, df.proc)
    
  }  
  
}

df.4.25$days <- as.numeric(difftime(df.4.25$datetime, min(df.4.25$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 30 C ######

t<-30 # Temp

df.4.30 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:32){ # Plates
  
  df.design.temp.plate <- df.design.4 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block4_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.4.30 <- rbind(df.4.30, df.proc)
    
  }  
  
}

df.4.30$days <- as.numeric(difftime(df.4.30$datetime, min(df.4.30$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 33 C ######

t<-33 # Temp

df.4.33 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.4 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block4_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.4.33 <- rbind(df.4.33, df.proc)
    
  }  
  
}

df.4.33$days <- as.numeric(difftime(df.4.33$datetime, min(df.4.33$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 35 C ######

t<-35 # Temp

df.4.35 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.4 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block4_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.4.35 <- rbind(df.4.35, df.proc)
    
  }  
  
}

df.4.35$days <- as.numeric(difftime(df.4.35$datetime, min(df.4.35$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 39 C ######

t<-39 # Temp

df.4.39 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.4 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block4_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.4.39 <- rbind(df.4.39, df.proc)
    
  }  
  
}

df.4.39$days <- as.numeric(difftime(df.4.39$datetime, min(df.4.39$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

###### 43 C ######

t<-43 # Temp

df.4.43 <- data.frame() # Empty data frame for adding the df.plt data. 

for (p in 1:2){ # Plates
  
  df.design.temp.plate <- df.design.4 %>%
    filter(Temperature.C == t, Plate.at.T == p)
  
  for (i in 0:15){ # Reads
    
    file.path <- paste("raw-data/JRL_block4_", t, "_", p, "_", i, ".xlsx", sep = "")
    # t is the temp, p is the plate, i is the read
    
    # Try to read the file; if error, skip to next i
    df.raw <- tryCatch(
      read_excel(file.path, range = "D50:M62") %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(df.raw)) {
      message("File missing: ", file.path, " - skipping")
      next
    }
    
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
    
    vars <- tryCatch(
      read_excel(file.path, range = "B7:B8", col_names = FALSE) %>% as.data.frame(),
      error = function(e) NULL
    )
    
    if (is.null(vars)) {
      message("Vars missing: ", file.path, " - skipping")
      next
    }
    
    vars<-as.data.frame(vars)
    
    date <- as.Date(vars[1, 1], format = "%Y-%m-%d") #date
    df.plt$date <- rep(julian(date, origin = as.Date("2030-01-01")), nrow(df.plt)) # Julian date
    
    time <- as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    df.plt$time <- rep(format(time, "%H:%M:%S"), nrow(df.plt)) # Extracted time
    
    df.plt$datetime <- as.POSIXct(paste(date, df.plt$time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    df.plt$Row <- as.integer(df.plt$Row)
    df.plt$Column <- as.integer(df.plt$Column)
    
    df.plt <- df.plt %>% 
      mutate(Well.at.T = (Plate-1)*60 + (Row-1)*10 + Column)
    
    df.proc <- df.plt %>%
      left_join(df.design.temp.plate, by = "Well.at.T") # This is the full data we want for each frame (don't need elapsed time)
    
    df.4.43 <- rbind(df.4.43, df.proc)
    
  }  
  
}

df.4.43$days <- as.numeric(difftime(df.4.43$datetime, min(df.4.43$datetime, na.rm = TRUE), units = "days")) # Calculate time passed since the start of the experiment.

df.sum <- rbind(df.1.8, df.1.14, df.1.20, df.1.25, df.1.30, df.1.33, df.1.35, df.1.39, df.1.43,
                df.2.8, df.2.14, df.2.20, df.2.25, df.2.30, df.2.33, df.2.35, df.2.39, df.2.43,
                df.3.8, df.3.14, df.3.20, df.3.25, df.3.30, df.3.33, df.3.35, df.3.39, df.3.43,
                df.4.8, df.4.14, df.4.20, df.4.25, df.4.30, df.4.33, df.4.35, df.4.39, df.4.43)

write.csv(df.sum, "processed-data/01_raw_timeseries_data.csv") # Save raw data. 133,320 observations!

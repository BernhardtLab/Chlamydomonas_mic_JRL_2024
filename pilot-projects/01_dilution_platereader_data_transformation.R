# Jason R Laurich
# November 28, 2024

# Script to process plate reader data for dilution experimental pilot
# Save readable, interpretable csv file

############### Packages #############################

library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)

############### Access the csv file ######################

df.design <- read.csv("pilot-projects/01_dilution_exp_design.csv") # all concentrations are µL, and are relative to pure cultures as a fraction of 1
# Also all 0.2 µL entries are actually 0 - was not able to pipette those volumes.

head(df.design) # This matches treatment to wells
str(df.design)  

df.rawd <- read.csv("pilot-projects/02.5_dilution_data.csv")
head(df.rawd)
str(df.rawd)

############### Process and transform data #####################

df.rawd$X[df.rawd$X == ""] <- NA # These 2 steps will fill in the empty spaces beneath A, B....
df.rawd$X <- zoo::na.locf(df.rawd$X)
df.rawd$meas <- rep(c("RFU", "OD700", "OD600", "OD750"), 8)

RFU <- numeric()
OD700 <- numeric()
OD600 <- numeric()
OD750 <- numeric()
Row <- character()
Col <- numeric()

# Loop through the df.rawd

for (i in 1:nrow(df.rawd)){
  
  x <- ceiling(i/4) # This function round all values above an integer to the next. I'll use this to only record one line for each row of wells.
  
  for (c in 2:13){
    if (i == 4*x){
      Row <- c(Row, df.rawd$X[x]) # Only record row and column every 4 rows.
      Col <- c(Col, c - 1) 
      }# Numerical value for the column} 
    
    if(df.rawd$meas[i] == 'RFU'){
      RFU <- c(RFU, df.rawd[i,c])
    }
    
    if(df.rawd$meas[i]=='OD700'){
      OD700 <- c(OD700, df.rawd[i,c])
    }
    
    if(df.rawd$meas[i]=='OD600'){
      OD600 <- c(OD600, df.rawd[i,c])
    }
    
    if(df.rawd$meas[i]=='OD750'){
      OD750 <- c(OD750, df.rawd[i,c])
    }
  }
}

df <- data.frame(
  Row = Row,
  Column = Col,
  RFU = RFU,
  OD700 = OD700,
  OD600 = OD600,
  OD750 = OD750
)

df <- df %>%
  left_join(df.design, by = c("Row", "Column")) # Merge the data

write.csv(df, "pilot-projects/03_dilution_data_clean.csv") # Save cleaned data

# OK, now we want to examine the microbial optical density stuff

blankOD600<- mean(df$OD600[1:3])
blankRFU<- mean(df$RFU[1:3])

ggplot(data = subset(df, V.Chlamy == 0 & V.Em1021 != 0 & V.Em1021 != 0.2), aes(x = C.Em1021, y = OD600 - blankOD600)) +
  geom_point() +
  labs(
    title = "OD600 vs. C.Em1021",
    x = "C.Em1021",
    y = "OD600"
  ) +
  ylim(-0.01,0.2) +
  theme_classic()

ggplot(data = subset(df, V.Chlamy == 0 & V.Em1022 != 0 & V.Em1022 != 0.2), aes(x = C.Em1022, y = OD600 - blankOD600)) +
  geom_point() +
  labs(
    title = "OD600 vs. C.Em1022",
    x = "C.Em1022",
    y = "OD600"
  ) +
  ylim(-0.01,0.2) +
  theme_classic()

# It looks like I the relationship is not quite linear here. I may need to proof it for each bacteria. The machine is definetly capable of measuring 
# These cells at 1% of their concentration yesterday! Might be worth just running a plate quick, seeing what it can pick up on.

ggplot(data = subset(df, V.Chlamy != 0 & V.Em1021 == 0 & V.Em1022 == 0), aes(x = OD600 - blankOD600, y = RFU - blankRFU)) +
  geom_point() +
  labs(
    title = "OD600 vs. RFU for Chlamy only",
    x = "RFU",
    y = "OD600"
  ) +
  theme_classic()




# Jason R Laurich
# November 28, 2024

# Script to analyse pilot data from the dilution experiments.
# For now, want to initially explore the relationship between Ensifer dilution and OD readings (600, 700, 750)
# Get a sense of the minimum detection threshold of our plate reader. 

############## Packages ##########################################

library(ggplot2)
library(gridExtra)
library(ggrepel)

############## Get and examine the data #######################

df <- read.csv("pilot-projects/03_dilution_data_clean.csv")

str(df)
head(df)

############## Exploration and plotting #######################

blankOD600 <- mean(df$OD600[1:3]) # Let's establish the baseline for all measurements in blank wells
blankOD700 <- mean(df$OD700[1:3])
blankOD750 <- mean(df$OD750[1:3])
blankRFU <- mean(df$RFU[1:3])

# First we want to know: what is the relationship between Chlamy RFU and optical density?

C.600 <- ggplot(data = subset(df, V.Chlamy != 0 & V.Em1021 == 0 & V.Em1022 == 0), aes(x = OD600 - blankOD600, y = RFU - blankRFU)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") + # Add linear regression line +
  labs(
    title = "OD600 vs. RFU for Chlamy only",
    x = "OD600",
    y = "RFU"
  ) +
  theme_classic()

C.700 <- ggplot(data = subset(df, V.Chlamy != 0 & V.Em1021 == 0 & V.Em1022 == 0), aes(x = OD700 - blankOD700, y = RFU - blankRFU)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") + # Add linear regression line +
  labs(
    title = "OD700 vs. RFU for Chlamy only",
    x = "OD700",
    y = "RFU"
  ) +
  theme_classic()

C.750 <- ggplot(data = subset(df, V.Chlamy != 0 & V.Em1021 == 0 & V.Em1022 == 0), aes(x = OD750 - blankOD750, y = RFU - blankRFU)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") + # Add linear regression line +
  labs(
    title = "OD750 vs. RFU for Chlamy only",
    x = "OD750",
    y = "RFU"
  ) +
  theme_classic()

grid.arrange(C.600, C.700, C.750) # Couple points are off for all (high OD, low RFU)

C.600 <- ggplot(data = subset(df, V.Chlamy != 0 & V.Em1021 == 0 & V.Em1022 == 0), 
                aes(x = OD600 - blankOD600, y = RFU - blankRFU, label = Well.ID)) +
  geom_point() +
  geom_text_repel(size = 3) + # Add labels with Well.ID, adjust size as needed
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") + # Add linear regression line
  labs(
    title = "OD600 vs. RFU for Chlamy only",
    x = "OD600",
    y = "RFU"
  ) +
  theme_classic()

C.600 <- ggplot(data = subset(df, V.Chlamy != 0 & V.Em1021 == 0 & V.Em1022 == 0, Well.ID != c('B5','A10','B9')), 
                aes(x = OD600 - blankOD600, y = RFU - blankRFU, label = Well.ID)) +
  geom_point() +
  geom_text_repel(size = 3) + # Add labels with Well.ID, adjust size as needed
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") + # Add linear regression line
  labs(
    title = "OD600 vs. RFU for Chlamy only",
    x = "OD600",
    y = "RFU"
  ) +
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


############# More sophisticated model exploration #####################
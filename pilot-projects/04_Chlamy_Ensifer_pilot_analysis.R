# Jason R Laurich
# December 4, 2024

# Going to look at and model the data for the Chlamydomonas x Enisfer pilot.

############### Packages #############################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(gridExtra)

############### Access the csv file ######################

df <- read.csv("pilot-projects/05_Chlamy_Ensifer_pilot_data.csv")

str(df)
head(df)

names(df)[13:16] <- c('chlamy','em1021','em1022','nut.lvl')
df$bac.lvl <- pmax(df$em1021, df$em1022, na.rm = TRUE) # Get value of bacterial inoculation

df[, 13:19] <- lapply(df[, 13:19], as.factor)

############## Run a quick model #########################

mod <- lmer(RFU~bac.lvl*nut.lvl*days + (1|bac), data=df)

summary(mod)

p <- ggplot(df, aes(x = as.numeric(as.character(days)), y = RFU, color = as.factor(bac.lvl))) +
  geom_line(aes(group = interaction(Row, Column, plate)), alpha = 0.7) + # Connect points for each well
  geom_point(size = 2) +
  facet_grid(bac ~ nut.lvl) + # Rows for bacterial treatments, columns for nutrient levels
  scale_color_manual(
    values = c(
      "0" = "springgreen4",
      "0.01" = "black",
      "0.1" = "magenta3",
      "1" = "goldenrod1"
    ),
    name = "Bacterial\nInoculation Level"
  ) +
  labs(
    title = "Fluorescence Over Time",
    x = "Elapsed Days",
    y = "RFU"
  ) +
  theme_minimal()

p

ggsave("pilot-projects/01_RFU_time_bac.pdf", p, width = 30, height = 20, units = "cm")



# Jason R Laurich
# February 2nd, 2026

# We're going to work with the raw timeseries data and estimate µ's from exponential growth curves
# We'll do this first for Chlamydomonas using RFUs, then look into calibrating a OD600~RFU relationship
# we can use to roughly back-calculate variation in microbial abundance when they are associated with Chlamy
# For the microbe-only wells, we'll just use raw OD600 data. 

# Load packages -----------------------------------------------------------

library(tidyverse)
library(cowplot)
library(nls.multstart)

# Upload the raw time series data --------------------------------

df <- read.csv("processed-data/01_raw_timeseries_data.csv") # Get the raw data.
head(df) # This matches treatment to wells

df %>% 
  filter(Chlamy.y.n == "BLANK") %>% 
  summarize(mean(OD600), mean(RFU)) # Mean OD for blanks is 1.328268, RFU 18.2394

df <- df %>% 
  filter (Chlamy.y.n != "BLANK") # 122808 observations after removing the blanks. 

### Let's identify the errors and decide what to do with them. 

#Hmm, for now I will keep samples where I caught pipetting errors. For example, we will keep the treatments where I noted that wells received bacteria 5, not 4
df.design.1 <- read.csv("raw-data/02-block-1-design.csv") # Get the experimental design file.
head(df.design.1) # This matches treatment to wells
df.design.1$Block <- 1

df.design.1 <- df.design.1 %>% 
  mutate(unique.id = paste0("b1.t", Temperature.C, ".p", Plate.at.T, ".w", Well.at.T))

df.design.1 %>% 
  filter(!is.na(Notes), trimws(Notes) != "") # Let's identify the problematic rows — e.g. where I left notes on poten

df <- df %>%
  
  mutate(
    Replicate = suppressWarnings(as.numeric(Replicate))
  ) %>%
  
  filter(!unique.id %in% c("b1.t25.p1.w19",
                           "b1.t33.p1.w2",
                           "b1.t33.p1.w12",
                           "b1.t39.p2.w105")) %>% 
  
  mutate(Microbe = case_when(
    unique.id %in% c("b1.t30.p27.w1615",
                     "b1.t30.p28.w1641",
                     "b1.t30.p28.w1649",
                     "b1.t30.p28.w1663",
                     "b1.t30.p28.w1671",
                     "b1.t30.p29.w1685",
                     "b1.t30.p29.w1699",
                     "b1.t30.p29.w1716",
                     "b1.t30.p29.w1717",
                     "b1.t30.p29.w1730",
                     "b1.t30.p29.w1735") ~ "5",
    TRUE ~ Microbe)) %>% 
  
  mutate(Replicate = case_when(
    unique.id %in% c("b1.t30.p27.w1615",
                     "b1.t30.p28.w1641",
                     "b1.t30.p28.w1649",
                     "b1.t30.p28.w1663",
                     "b1.t30.p28.w1671",
                     "b1.t30.p29.w1685",
                     "b1.t30.p29.w1699",
                     "b1.t30.p29.w1716",
                     "b1.t30.p29.w1717",
                     "b1.t30.p29.w1730") ~ 4,
    
    unique.id == "b1.t30.p29.w1735" ~ 5,
    TRUE ~ Replicate)) 


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

# We'll change the treatment information for the 2 errors. 

df <- df %>%
  
  mutate(
    Replicate = suppressWarnings(as.numeric(Replicate))
  ) %>%
  
  mutate(Microbe = case_when(
    unique.id == "b3.t30.p7.w365" ~ "8",
    unique.id == "b3.t30.p26.w1542" ~ "4",
    TRUE ~ Microbe)) %>% 
  
  mutate(Replicate = case_when(
    unique.id %in% c("b3.t30.p7.w365",
                     "b3.t30.p26.w1542") ~ 4,
    TRUE ~ Replicate)) 

df.design.4 <- read.csv("raw-data/08-block-4-design.csv") # Get the experimental design file.
head(df.design.4) # This matches treatment to wells
df.design.4$Block <- 4

df.design.4 <- df.design.4 %>% 
  mutate(unique.id = paste0("b4.t", Temperature.C, ".p", Plate.at.T, ".w", Well.at.T))

df.design.4 %>% 
  filter(!is.na(Notes), trimws(Notes) != "") # Let's identify the problematic rows — e.g. where I left notes on potential errors etc. 

# OK we'll fix these errors. I'm going to take out the one that got some unknown number of bacteria, as well as the one that got chlamy accidentally
df <- df %>%
  
  mutate(
    Replicate = suppressWarnings(as.numeric(Replicate))
  ) %>%
  
  filter(!unique.id %in% c("b4.t30.p19.w1085",
                           "b4.t30.p20.w1152")) %>% 
  
  mutate(Microbe = case_when(
    unique.id == "b4.t30.p30.w1764" ~ "10",
    TRUE ~ Microbe)) %>% 
  
  mutate(Replicate = case_when(
    unique.id == "b4.t30.p30.w1764" ~ 4,
    TRUE ~ Replicate)) 

# So in the end, we have 122,808 time series data after removing errors and blanks.

# Estimate µ for Chlamy! --------------------------------------------------

df.chlamy <- df %>% 
  filter(Chlamy.y.n == 'y') %>% 
  mutate(log.RFU = log(RFU + 0.001))

# So we have unique replicate ids in unique.id, temp in Temperature.C, salt in Salt.conc.g.l, and nitrogen in Nitrogen.conc.µM...
# We'll need to record each of these for each replicate, as well the exponential growth rate. We also want block, and microbe, plate, and well

df.µ <- data.frame(                    # Summary dataframe for µ estimates
  unique.id = character(),             # Unique ID
  temp = numeric(),                    # Temperature
  nit = numeric(),                     # Nitrogen
  salt = numeric(),                    # Salt
  block = numeric(),                   # Block
  mic = character(),                   # Microbe
  plate = numeric(),                   # Plate at Temp
  well = numeric(),                    # Well at Temp
  rep = numeric(),                    # Replicate
  µ = numeric()                        # Thresholded µ
)

# Each unique well ID has its own singular combination of treatments.

for (i in unique(df.chlamy$unique.id)){ # for every ID (unique to block, temp, plate, and well) {
  
  df.i <- df.chlamy %>% 
    filter(unique.id == i) 
  
  df.i<- df.i[order(df.i$days), ]
  
  if (df.i$RFU[2] < df.i$RFU[1]) {
    df.i <- df.i[-1, ]
  } # So if there is a drop in RFUs between the first and second read, we will treat the second data point as N0.
  
  df.i <- df.i %>% 
    mutate(N0 = RFU[1])
  
  t.series <- unique(df.i$days) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
  t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
  
  ln.slopes <- c() # Re-initialize this too!
  
  for (z in t.series){
    
    df.i.sl <- df.i[df.i$days <= z, ] # Subset the data to exclude time points above our window
    
    ln_slope <- lm(log.RFU~days, data = df.i.sl)
    
    ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
    
  } # So now we have our slopes fit to the logged RFU data
  
  s <- max(2, which.max(ln.slopes))  # We need at least 3 data points
  
  df.i.th <- df.i[df.i$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
  # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
  
  if (length(unique(na.omit(df.i.th$RFU))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
    µ.est <- 0
    
  } else {
  
    µ.mod <- tryCatch(
      nls_multstart(
        RFU ~ N0 * exp(r * days),
        data = df.i.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      ),
      error = function(e) NULL
    )
    
    if (is.null(µ.mod)) {
      µ.est <- NA_real_
    } else {
      µ.est <- coef(µ.mod)[["r"]]  
    }
  }
  
  df.µ <- rbind(df.µ, data.frame(
    unique.id = df.i$unique.id[1],
    temp = df.i$Temperature.C[1],
    nit = df.i$Nitrogen.conc.µM[1],
    salt = df.i$Salt.conc.g.l[1],
    block = df.i$Block[1],
    mic = df.i$Microbe[1],
    plate = df.i$Plate.at.T[1],
    well = df.i$Well.at.T[1],
    rep = df.i$Replicate[1],
    µ = µ.est
  ))
  
} # 5306 estimates of µ

df.µ %>% filter(is.na(µ)) # How did this work? Are there any NA values that should be fixed?
# Looks good! No NA's

write.csv(df.µ, "processed-data/02_chlamy_µs.csv") #5306 measurements

# Estimate µ for microbes alone --------------------------------------------------

df.mic <- df %>% 
  filter(Chlamy.y.n == 'n') %>% 
  mutate(log.OD = log(OD600 + 0.001))

# So we have unique replicate ids in unique.id, temp in Temperature.C, salt in Salt.conc.g.l, and nitrogen in Nitrogen.conc.µM...
# We'll need to record each of these for each replicate, as well the exponential growth rate. We also want block, and microbe, plate, and well

df.µ <- data.frame(                    # Summary dataframe for µ estimates
  unique.id = character(),             # Unique ID
  temp = numeric(),                    # Temperature
  nit = numeric(),                     # Nitrogen
  salt = numeric(),                    # Salt
  block = numeric(),                   # Block
  mic = character(),                   # Microbe
  plate = numeric(),                   # Plate at Temp
  well = numeric(),                    # Well at Temp
  rep = numeric(),                     # Replicate
  µ = numeric()                        # Thresholded µ
)

# Each unique well ID has its own singular combination of treatments.

for (i in unique(df.mic$unique.id)){ # for every ID (unique to block, temp, plate, and well) {
  
  df.i <- df.mic %>% 
    filter(unique.id == i) 
  
  df.i<- df.i[order(df.i$days), ]
  
  if (df.i$OD600[2] < df.i$OD600[1]) {
    df.i <- df.i[-1, ]
  } # So if there is a drop in OD600 between the first and second read, we will treat the second data point as N0.
  
  df.i <- df.i %>% 
    mutate(N0 = OD600[1])
  
  t.series <- unique(df.i$days) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
  t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
  
  ln.slopes <- c() # Re-initialize this too!
  
  for (z in t.series){
    
    df.i.sl <- df.i[df.i$days <= z, ] # Subset the data to exclude time points above our window
    
    ln_slope <- lm(log.OD~days, data = df.i.sl)
    
    ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
    
  } # So now we have our slopes fit to the logged RFU data
  
  s <- max(2, which.max(ln.slopes))  # We need at least 3 data points
  
  df.i.th <- df.i[df.i$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
  # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
  
  if (length(unique(na.omit(df.i.th$OD600))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
    µ.est <- 0
    
  } else {
    
    µ.mod <- tryCatch(
      nls_multstart(
        OD600 ~ N0 * exp(r * days),
        data = df.i.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      ),
      error = function(e) NULL
    )
    
    if (is.null(µ.mod)) {
      µ.est <- NA_real_
    } else {
      µ.est <- coef(µ.mod)[["r"]]  
    }
  }
  
  df.µ <- rbind(df.µ, data.frame(
    unique.id = df.i$unique.id[1],
    temp = df.i$Temperature.C[1],
    nit = df.i$Nitrogen.conc.µM[1],
    salt = df.i$Salt.conc.g.l[1],
    block = df.i$Block[1],
    mic = df.i$Microbe[1],
    plate = df.i$Plate.at.T[1],
    well = df.i$Well.at.T[1],
    rep = df.i$Replicate[1],
    µ = µ.est
  ))
  
} # 4979 estimates of µ

df.µ %>% filter(is.na(µ)) # How did this work? Are there any NA values that should be fixed?
# Looks good! No NA's

write.csv(df.µ, "processed-data/03_microbe_alone_µs.csv") #5306 measurements

# Estimating microbial µ with Chlamy --------------------------------------

df.ch.alone <- df %>% 
  filter(Chlamy.y.n =='y', Microbe == 'none')

df.ch.alone %>% 
  filter(days >= 1) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = days)) +
  geom_point() +
  facet_wrap(~ Block, labeller = labeller(Block = function(x) paste("Block", x))) +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of day?

df.ch.alone %>% 
  filter(days >= 1) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = Salt.conc.g.l)) +
  geom_point() +
  facet_wrap(~ Block, labeller = labeller(Block = function(x) paste("Block", x))) +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of salt?

df.ch.alone %>% 
  filter(days >= 1) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = Nitrogen.conc.µM)) +
  geom_point() +
  facet_wrap(~ Block, labeller = labeller(Block = function(x) paste("Block", x))) +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of nitrogen?

df.ch.alone %>% 
  filter(days >= 1) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = Temperature.C)) +
  geom_point() +
  facet_wrap(~ Block, labeller = labeller(Block = function(x) paste("Block", x))) +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of temperature?

df.ch.alone %>% 
  filter(days >= 1) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = Temperature.C)) +
  geom_point() +
  facet_wrap(~ Block, labeller = labeller(Block = function(x) paste("Block", x))) +
  geom_smooth(aes(group = Temperature.C),
            method = "lm",
            se = FALSE) +
  theme_classic() # Effects of temperature?

df.ch.alone %>% 
  filter(days >= 1) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = Block)) +
  geom_point() +
  facet_wrap(~ Temperature.C) +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of temperature?


# Jason R Laurich
# April 25th, 2026

# We're going to work with the raw timeseries data and estimate µ's from exponential growth curves
# for microbes on their own first, then we will explore fitting µ to microbes in mixed wells with chlamy.

# Load packages -----------------------------------------------------------

library(tidyverse)
library(cowplot)
library(nls.multstart)
library(lme4)

# Upload the raw time series data --------------------------------

df <- read.csv("processed-data/01_raw_timeseries_data.csv") # Get the raw data.
head(df) # This matches treatment to wells

df %>% 
  filter(Chlamy.y.n == "BLANK") %>% 
  summarize(mean(OD600), mean(RFU)) # Mean OD for blanks is 1.328268, RFU 18.2394

df <- df %>% 
  filter (Chlamy.y.n != "BLANK") # 122894 observations after removing the blanks. 

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
  
  µ = numeric(),                       # Thresholded µ
  k = numeric(),                       # Maximum density (ish, proxy for carrying capacity?)
  time = numeric()                     # How many days are needed to estimate the exponential growth phase?
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
    
    µ = µ.est,
    k = df.i %>% filter(days > 1) %>% summarize(k = max(OD600, na.rm = T)) %>% pull(k),
    time = max(df.i.th$days)                  
  ))
  
} # 5306 estimates of µ

df.µ %>% filter(is.na(µ)) # How did this work? Are there any NA values that should be fixed?
# Looks good! No NA's

write.csv(df.µ, "processed-data/13_mics_alone_µs.csv") #5306 measurements

# Examine the data! -------------------------------------------------------

df.µ.m <- read.csv("processed-data/13_mics_alone_µs.csv")

ggplot(df.µ.m %>% filter(salt == 0, nit == 1000), aes(x = temp, y = µ)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ mic, ncol = 6) +
  theme_classic()

ggplot(df.µ.m %>% filter(temp == 30, nit == 1000), aes(x = salt, y = µ)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ mic, ncol = 6) +
  theme_classic()

ggplot(df.µ.m %>% filter(temp == 30, salt == 0), aes(x = nit, y = µ)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ mic, ncol = 6) +
  theme_classic()

# Look at variation in ~ k

ggplot(df.µ.m %>% filter(salt == 0, nit == 1000), aes(x = temp, y = k)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ mic, ncol = 6) +
  geom_smooth(method = 'lm') +
  theme_classic()

ggplot(df.µ.m %>% filter(temp == 30, nit == 1000), aes(x = salt, y = k)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ mic, ncol = 6) +
  geom_smooth(method = 'lm') +
  theme_classic()

ggplot(df.µ.m %>% filter(temp == 30, salt == 0), aes(x = nit, y = k)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ mic, ncol = 6) +
  geom_smooth(method = 'lm') +
  theme_classic()

# Estimate µ for microbes with chlamy --------------------------------------------------

# We need to explore the relationship b/w RFU and OD600 in chlamy-only wells.

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

# OK so I think we will estimate each block and each t separately.
# And we'll match up the data with the time period used for estimating µ in the underlying data. 

df.µ.c <- read.csv("processed-data/02_chlamy_µs.csv")
head(df.µ.c)

# Block 1

df.ch.1 <- df.ch.alone %>% 
  filter(Block == 1) %>% 
  left_join(
    df.µ.c %>% 
      select(unique.id, time),
    by = "unique.id"
  )

df.ch.1 %>% 
  filter(days >= 0.5) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = days)) +
  geom_point() +
  facet_wrap(~ Temperature.C, scales = "free_x") +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of temperature?

df.ch.1 %>% 
  filter(abs(days - time.y) <= 2) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = days)) +
  geom_point() +
  facet_wrap(~ Temperature.C, scales = "free_x") +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of temperature?

# Block 3

df.ch.3 <- df.ch.alone %>% 
  filter(Block == 3) %>% 
  left_join(
    df.µ.c %>% 
      select(unique.id, time),
    by = "unique.id"
  )

df.ch.3 %>% 
  filter(days >= 0.5) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = days)) +
  geom_point() +
  facet_wrap(~ Temperature.C, scales = "free_x") +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of temperature?

df.ch.3 %>% 
  filter(abs(days - time.y) <= 2) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = days)) +
  geom_point() +
  facet_wrap(~ Temperature.C, scales = "free_x") +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of temperature?

# Block 4

df.ch.4 <- df.ch.alone %>% 
  filter(Block == 4) %>% 
  left_join(
    df.µ.c %>% 
      select(unique.id, time),
    by = "unique.id"
  )

df.ch.4 %>% 
  filter(days >= 0.5) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = days)) +
  geom_point() +
  facet_wrap(~ Temperature.C, scales = "free_x") +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of temperature?

df.ch.4 %>% 
  filter(abs(days - time.y) <= 2) %>% 
  ggplot(aes(y = OD600, x = RFU, colour = days)) +
  geom_point() +
  facet_wrap(~ Temperature.C, scales = "free_x") +
  geom_smooth(method = 'lm') +
  theme_classic() # Effects of temperature?

###### Fitting OD600 ~ RFU relationships in chlamy-only wells ######

df.ch.1.cal <- df.ch.1 %>%
  
  filter(days >= 0.5) %>% 
  
  group_by(Block, Temperature.C) %>%
  
  summarise(
    n = n(),
    slope = coef(lm(OD600 ~ RFU))[2],
    intercept = coef(lm(OD600 ~ RFU))[1],
    r2 = summary(lm(OD600 ~ RFU))$r.squared,
    .groups = "drop"
  )

df.ch.1.cal

mod.cal <- lmer(
  OD600 ~ RFU +
    (RFU | Temperature.C) +
    (1 | Block),
  data = df.ch.alone
)

summary(mod.cal)

mod.cal2 <- lmer(
  OD600 ~ RFU +
    (1 | Temperature.C) +
    (1 | Block),
  data = df.ch.alone
)

summary(mod.cal2)

mod.cal3 <- lmer(
  OD600 ~ RFU +
    (1 | Temperature.C) +
    (0 + RFU | Temperature.C) +
    (1 | Block),
  data = df.ch.alone
)

summary(mod.cal3)

###### Model selection ######

df.ch.th <- df.ch.alone %>% 
  filter(days >= 0.5)

mod.1 <- lmer(
  OD600 ~ RFU +
    (1 | Block),
  data = df.ch.th
)

mod.2 <- lmer(
  OD600 ~ RFU*Temperature.C +
    (1 | Block),
  data = df.ch.th
)

mod.3 <- lmer(
  OD600 ~ RFU*factor(Temperature.C) +
    (1 | Block),
  data = df.ch.th
)

mod.4 <- lmer(
  OD600 ~ RFU*poly(Temperature.C,2) +
    (1 | Block),
  data = df.ch.th
)

AIC(mod.1, mod.2, mod.3, mod.4)
BIC(mod.1, mod.2, mod.3, mod.4)

newdat <- expand.grid(
  RFU = seq(min(df.ch.th$RFU),
            max(df.ch.th$RFU),
            length.out = 100),
  Temperature.C = sort(unique(df.ch.th$Temperature.C)),
  Block = NA
)

newdat$pred4 <- predict(mod.4, newdata = newdat, re.form = NA)
newdat$pred3 <- predict(mod.3, newdata = newdat, re.form = NA)

ggplot(newdat, aes(RFU, pred4)) +
  geom_line(colour = "blue", linewidth = 1) +
  facet_wrap(~Temperature.C, scales = "free_x") +
  theme_classic()

ggplot(newdat) +
  geom_line(aes(RFU, pred4), colour = "blue") +
  geom_line(aes(RFU, pred3), colour = "red", linetype = "dashed") +
  facet_wrap(~Temperature.C, scales = "free_x") +
  theme_classic()

mod.5 <- lmer(
  OD600 ~ RFU + factor(Temperature.C) +
    (1 | Block),
  data = df.ch.th
)

AIC(mod.1, mod.2, mod.3, mod.4, mod.5)
BIC(mod.1, mod.2, mod.3, mod.4, mod.5)

# So model 5, where temp only affects the intercept, strongly wins out. We'll use that for now...

summary(mod.5)

df.ch.mic <- df %>% 
  filter(Chlamy.y.n == "y", 
         Microbe != "none")

head(df.ch.mic)

df.ch.mic <- df.ch.mic %>%
  mutate(
    algal.OD.pred = predict(mod.5, newdata = ., re.form = NULL)
  )

df.ch.mic <- df.ch.mic %>%
  mutate(
    microbe.OD.est = OD600- algal.OD.pred
  )

plot(microbe.OD.est ~ days, data = df.ch.mic)

###### Let's try estimating this seperately for each temp in each block! ######

df.cal <- df.ch.alone %>%
  filter(days >= 0.5)

mods <- df.cal %>%
  group_by(Block, Temperature.C) %>%
  group_split()

mod.list <- list()

for(i in seq_along(mods)) {
  
  df.i <- mods[[i]]
  
  key <- paste0("B", unique(df.i$Block),
                "_T", unique(df.i$Temperature.C))
  
  mod.list[[key]] <- lm(OD600 ~ RFU, data = df.i)
}

df.ch.mic$algal.OD.local <- NA_real_

for(nm in names(mod.list)) {
  
  mod <- mod.list[[nm]]
  
  b <- as.numeric(sub("B(\\d+)_.*", "\\1", nm))
  t <- as.numeric(sub(".*_T(\\d+)", "\\1", nm))
  
  idx <- df.ch.mic$Block == b &
    df.ch.mic$Temperature.C == t
  
  df.ch.mic$algal.OD.local[idx] <-
    predict(mod, newdata = df.ch.mic[idx, ])
}

df.ch.mic <- df.ch.mic %>%
  mutate(
    microbe.local = OD600 - algal.OD.local
  )

test <- lmer(
  OD600 ~ RFU + Microbe * Temperature.C + (1|Block),
  data = df.ch.mic
)

summary(test)

# OK so after all this exploration, I'm going to work with (for now) the data from model 5, shifting values up so that all are positive based on the minimum corrected OD
# We will likely have to validate this with microscopy, but it's the only way we're going to be able to assess growth curves for microbes...

# Get the data ready for growth rate estimation --------------------------

# df.ch.mic

df.ch.mic <- df.ch.mic %>% 
  mutate(OD.corr = microbe.OD.est - min(microbe.OD.est) + 0.001)

df.ch.mic <- df.ch.mic %>% 
  mutate(log.OD = log(OD.corr))

# Fit the curves ----------------------------------------------------------

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
  
  µ = numeric(),                       # Thresholded µ
  k = numeric(),                       # Maximum density (ish, proxy for carrying capacity?)
  time = numeric()                     # How many days are needed to estimate the exponential growth phase?
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
    
    µ = µ.est,
    k = df.i %>% filter(days > 1) %>% summarize(k = max(OD600, na.rm = T)) %>% pull(k),
    time = max(df.i.th$days)                  
  ))
  
} # 5306 estimates of µ

df.µ %>% filter(is.na(µ)) # How did this work? Are there any NA values that should be fixed?
# Looks good! No NA's

write.csv(df.µ, "processed-data/14_mics_with_chlamy_µs.csv") #4976 measurements




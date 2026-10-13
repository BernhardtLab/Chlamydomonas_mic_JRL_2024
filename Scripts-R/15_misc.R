
# summary on the whole peak abd thing?

library(tidyverse)

ts <- read.csv("Data-processed/01_timeseries_data.csv")

decline <- ts %>%
  filter(carbon.added, algae == "present") %>%   # main-text wells with Chlamy
  group_by(well_id) %>%                          # unique well identifier
  arrange(days, .by_group = TRUE) %>%
  summarise(
    peak          = max(RFU),
    final         = last(RFU),
    peak_interior = which.max(RFU) < n(),        # peak before the last read
    grew          = max(RFU) > 2 * first(RFU)    # define a "grower"
  ) %>%
  filter(grew) %>%
  mutate(drop = 1 - final / peak)

median(decline$drop)        # the "~74%" figure
mean(decline$drop < 0.10)   # share that roughly plateaued (<10% decline)
mean(decline$peak_interior) # the "92.5%" figure
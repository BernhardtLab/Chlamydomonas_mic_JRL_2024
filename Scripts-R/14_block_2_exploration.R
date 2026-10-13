# Description -------------------------------------------------------------
#
# 14_block_2_exploration.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# October 11, 2026
#
# An exploration of the consequences of omitting organic carbon in the COMBO
# media for block 2. 
#
# What this script produces:
#
# Figures-main:
#   07_concept_fig.png            conceptual figure - microbial effects on TPC and interactions with stress gradients
#   
# Load packages -----------------------------------------------------------

library(tidyverse)

# Set paths ---------------------------------------------------------------

fig.dir <- "Figures-misc"
proc.dir     <- "Data-processed"

mu <- read.csv(file.path(proc.dir, "02_chlamy_mus.csv")) %>%
  filter(!is.na(mu)) %>%
  mutate(carbon = ifelse(carbon.added, "carbon added", "no carbon"),
         grp = case_when(mic == "none" ~ "alone",
                         mic == "all"  ~ "community",
                         TRUE          ~ "single partner"))

gmean <- function(d) d %>% group_by(carbon, grp) %>%
  summarise(m = mean(mu), se = sd(mu)/sqrt(n()), n = n(), .groups = "drop")

# Effects of carbon under benign conditions -------------------------------

effect_of <- function(tab) {
  alone <- tab %>% filter(grp == "alone") %>% select(carbon, m0 = m, se0 = se)
  tab %>% filter(grp != "alone") %>% left_join(alone, by = "carbon") %>%
    mutate(effect = m - m0, se_eff = sqrt(se^2 + se0^2))
}

eff.ref <- effect_of(gmean(filter(mu, temp == 30, nit == 1000, salt == 0)))
cat("\n== microbial growth effect at benign reference (mu, with - without) ==\n")
print(as.data.frame(eff.ref |> mutate(across(where(is.numeric), ~round(.x, 3)))), row.names = FALSE)

fig1 <- ggplot(eff.ref, aes(effect, grp, colour = carbon)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = effect - se_eff, xmax = effect + se_eff),
                  position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("carbon added" = "firebrick3", "no carbon" = "steelblue")) +
  labs(x = expression("Effect on algal growth rate ("*mu*", with - without microbe, day"^-1*")"),
       y = NULL, colour = NULL,
       title = "Microbial growth effect at the benign reference") +
  theme_bw()

fig1

ggsave(file.path(fig.dir, "20_carbon_effect_summ.png"), fig1, width = 7.5, height = 3, dpi = 300)

# Gradient effects --------------------------------------------------------

grad_eff <- function(d, xvar) {
  s <- d %>% rename(level = all_of(xvar)) %>%
    group_by(carbon, grp, level) %>%
    summarise(m = mean(mu), se = sd(mu)/sqrt(n()), .groups = "drop")
  alone <- s %>% filter(grp == "alone") %>% select(carbon, level, m0 = m)
  s %>% filter(grp != "alone") %>% left_join(alone, by = c("carbon", "level")) %>%
    mutate(effect = m - m0)
}

eg <- bind_rows(
  grad_eff(filter(mu, salt == 0,  nit == 1000), "temp") %>% mutate(gradient = "Temperature (C)"),
  grad_eff(filter(mu, temp == 30, salt == 0),   "nit")  %>% mutate(gradient = "Nitrogen (uM)"),
  grad_eff(filter(mu, temp == 30, nit == 1000), "salt") %>% mutate(gradient = "Salt (g/L)")) %>%
  mutate(gradient = factor(gradient, levels = c("Temperature (C)", "Nitrogen (uM)", "Salt (g/L)")))

fig2 <- ggplot(eg, aes(level, effect, colour = carbon, linetype = grp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() + geom_point(size = 1) +
  facet_wrap(~ gradient, scales = "free_x") +
  scale_colour_manual(values = c("carbon added" = "firebrick3", "no carbon" = "steelblue")) +
  scale_linetype_manual(values = c("community" = "solid", "single partner" = "22")) +
  labs(x = NULL, y = expression("Microbial growth effect ("*mu*", day"^-1*")"),
       colour = NULL, linetype = NULL,
       title = "Microbial effect on algal growth rate across gradients: carbon vs no carbon") +
  theme_bw() + theme(legend.position = "top")

fig2

ggsave(file.path(fig.dir, "21_carbon_effect_gradients.png"), fig2, width = 12, height = 4.2, dpi = 300)

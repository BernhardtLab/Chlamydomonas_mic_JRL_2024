# Description -------------------------------------------------------------
#
# 03_growth_validation.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# September 15, 2026
#
# Validate the per-well growth-rate estimate (mu) from 02 and decide what, if
# anything, to use as a yield / "carrying capacity" trait. This is a diagnostic
# script: it reads the outputs of 01 + 02, produces validation figures, and
# prints the numbers behind each check. It does NOT change mu.
#
# What this script produces:
#
# Figures-misc/
#   01_Fig_mu_estimation_check.png   mu (nls) vs a plain log-linear slope
#   02_Fig_mu_gradient_check.png     mu of the alga alone across T / N / salt
#   03_Fig_RFU_curve_shapes.png      example trajectories (boom-bust vs plateau)
#   04_Fig_peak_RFU_capture.png      is the peak actually caught, or right-censored?
# Data-processed/
#   03_estimator_check.csv           per-well mu vs log-linear slope (record)
#
# The two questions this answers:
#
# (1) Is mu suitable? Exponential growth is a straight line on a log scale, so
#     mu should equal the slope of log(RFU) vs time over the growth phase. We
#     (a) recompute that slope with a plain lm() and check it agrees with the
#     nls estimate 02 uses, and (b) confirm mu traces sensible biology across
#     the three gradients. We KEEP the nls estimate as primary; the log-linear
#     slope is computed only to document how much the estimator choice matters.
# (2) Can we use a yield trait? Only if the curves plateau (a real carrying
#     capacity K). They do not — they boom then bust — so there is no K. What
#     is defensible is peak RFU (a transient PEAK YIELD), provided the peak is
#     actually captured rather than still rising at the last read. We validate
#     that here and carry peak RFU (already in 02) as the yield trait.
#
# Load packages -----------------------------------------------------------

library(tidyverse)

# Set paths ---------------------------------------------------------------
proc.dir <- "Data-processed"
fig.dir  <- "Figures-misc"
raw.file <- "01_timeseries_data.csv"
mu.file  <- "02_chlamy_mus.csv"

# Load the data -----------------------------------------------------------

raw <- read.csv(file.path(proc.dir, raw.file))
mu  <- read.csv(file.path(proc.dir, mu.file))

# Chlamy trajectories only, sorted within well.
ch <- raw %>%
  filter(Chlamy.y.n == "y") %>%
  arrange(unique.id, days)

# 03 validates 02's estimates, so restrict to the wells 02 kept (it drops a few
# pipetting-error wells that still exist in the raw 01 data).
ch <- ch %>% filter(unique.id %in% mu$unique.id)

# Same first-read trim 02 uses: if read 2 < read 1 (settling), drop read 1.
trim_first <- function(sub) {
  sub <- sub[order(sub$days), ]
  if (nrow(sub) > 1 && !is.na(sub$RFU[2]) && sub$RFU[2] < sub$RFU[1]) sub <- sub[-1, ]
  sub
}

# Check 1: mu (nls method) v log-linear slope over the same window --------

# For each well: take the reads up to the window endpoint 02 selected (mu$time),
# and fit a straight line to log(RFU) vs days. Its slope is the exponential rate
# by the textbook definition. We compare it to the nls mu 02 recorded.

time.by.id <- mu %>% select(unique.id, window_end = time)
ch.win <- ch %>% left_join(time.by.id, by = "unique.id")
wells  <- split(ch.win, ch.win$unique.id)

loglin_slope <- function(sub) {
  sub <- trim_first(sub)
  te  <- sub$window_end[1]                    # window endpoint from 02 (renamed on the join)
  w   <- sub[sub$days <= te, ]
  if (nrow(w) < 2) return(NA_real_)
  unname(coef(lm(log(RFU + 0.001) ~ days, data = w))[2])
}

check <- data.frame(
  unique.id    = names(wells),
  loglin_slope = vapply(wells, loglin_slope, numeric(1)),
  row.names    = NULL) %>%
  left_join(select(mu, unique.id, mu, temp), by = "unique.id") %>%
  filter(!is.na(loglin_slope))

r.est <- cor(check$mu, check$loglin_slope)

write.csv(check, file.path(proc.dir, "03_estimator_check.csv"), row.names = FALSE)

fig1 <- ggplot(check, aes(mu, loglin_slope)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(alpha = 0.20, colour = "orchid3") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = sprintf("r = %.3f\nmedian |diff| = %.3f",
                           r.est, median(abs(check$mu - check$loglin_slope)))) +
  labs(x = "mu from nls_multstart (primary, raw scale)",
       y = "slope of lm(log(RFU) ~ days)\n(log-scale)",
       title = "Estimator check: nls mu vs a plain log-linear slope",
       subtitle = "nls kept as primary; log-linear shown to document sensitivity") +
  theme_classic()

ggsave(file.path(fig.dir, "01_Fig_mu_estimation_check.png"), fig1,
       width = 6.5, height = 5.5, dpi = 300)  # save the figure


# Check 2: does mu reasonably trace expected biology? ---------------------
# Algae alone, all gradients

# If mu is real, the alga-alone curves should be a thermal hump, a rise-and-
# saturate over nitrogen, and a decline over salt. Reference levels: N = 1000
# uM, salt = 0 g/L; the N and salt gradients are run at 30 C.

aa <- mu %>% filter(mic == "none")
bio <- bind_rows(
  aa %>% filter(nit == 1000, salt == 0)  %>% transmute(axis = "Temperature (C)", level = temp, mu),
  aa %>% filter(temp == 30, salt == 0)   %>% transmute(axis = "Nitrogen (uM)",   level = nit,  mu),
  aa %>% filter(temp == 30, nit == 1000) %>% transmute(axis = "Salt (g/L)",      level = salt, mu)
) %>%
  group_by(axis, level) %>%
  summarise(mean_mu = mean(mu),
            se      = sd(mu) / sqrt(n()), .groups = "drop")

fig2 <- ggplot(bio, aes(level, mean_mu)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(colour = "orchid3") +
  geom_point(colour = "orchid3") +
  geom_errorbar(aes(ymin = mean_mu - se, ymax = mean_mu + se),
                width = 0, colour = "orchid3") +
  facet_wrap(~ axis, scales = "free_x") +
  labs(x = "gradient level", y = "mean mu (alga alone)",
       title = "Does mu trace biology? Alga-alone growth across the three gradients") +
  theme_bw()

ggsave(file.path(fig.dir, "02_Fig_mu_gradient_check.png"), fig2,
       width = 9, height = 3.6, dpi = 300) # Save the figure. 


# Check 3: Do growth curves plateau (real K) or boom-bust (peak ab --------

# Example alga-alone trajectories (block 1, benign N & salt) across all
# temperatures. If they never flatten, there is no carrying capacity.

ex <- raw %>%
  filter(Chlamy.y.n == "y", Microbe == "none", Block == 1,
         Nitrogen.conc.uM == 1000, Salt.conc.g.l == 0,
         Temperature.C %in% c(8, 14, 20, 25, 30, 33, 35, 39, 43)) %>%
  group_by(Temperature.C) %>%
  filter(Well.at.T == min(Well.at.T)) %>%       # one representative well per temp
  ungroup() %>%
  mutate(Temperature.C = factor(Temperature.C))

fig3 <- ggplot(ex, aes(days, RFU, colour = Temperature.C)) +
  geom_line() + geom_point(size = 1) +
  scale_colour_viridis_d(name = "Temp (C)", end = 0.9) +
  labs(x = "days", y = "RFU (alga alone)",
       title = "Curve shapes: boom-bust, not plateau (so no carrying capacity K)") +
  theme_classic()

ggsave(file.path(fig.dir, "03_Fig_RFU_curve_shapes.png"), fig3,
       width = 6.5, height = 4.5, dpi = 300)

# Quantify the crash: fractional drop from peak to final read.
crash_frac <- function(sub) {
  sub <- trim_first(sub)
  if (nrow(sub) < 3) return(NA_real_)
  pk <- max(sub$RFU); (pk - sub$RFU[nrow(sub)]) / pk
}

growers <- mu$unique.id[mu$mu > 0.1]
cf <- vapply(wells[names(wells) %in% growers], crash_frac, numeric(1))
cf <- cf[!is.na(cf)]
cat(sprintf("[Check 3] growers: median crash from peak = %.0f%%; plateau (<10%% drop) = %.0f%%\n",
            100 * median(cf), 100 * mean(cf < 0.1)))

# Check 4: validating peak RFU as a real measure? -------------------------

# peak RFU is only a trustworthy yield trait if the maximum falls in the
# interior of the series. A max at the LAST read means growth had not finished
# (right-censored -> peak underestimated); a max at the FIRST read means the
# well never grew.

peak_position <- function(sub) {
  sub <- trim_first(sub)
  i <- which.max(sub$RFU); n <- nrow(sub)
  if (i == 1) "first read (no growth)"
  else if (i == n) "last read (still rising)"
  else "interior (peak captured)"
}

pk <- data.frame(
  unique.id  = names(wells),
  peak_where = vapply(wells, peak_position, character(1)),
  row.names  = NULL) %>%
  left_join(select(mu, unique.id, mu), by = "unique.id") %>%
  mutate(group = ifelse(mu > 0.1, "growers (mu > 0.1)", "non-growers"))

cat("[Check 4] peak position among growers:\n")
print(round(100 * prop.table(table(
  pk$peak_where[pk$group == "growers (mu > 0.1)"])), 1)) # So this works for 92.5% of curves.

fig4 <- ggplot(pk, aes(peak_where, fill = peak_where)) +
  geom_bar() +
  facet_wrap(~ group) +
  scale_fill_manual(values = c("first read (no growth)"   = "#d95f02",
                               "interior (peak captured)" = "#1b9e77",
                               "last read (still rising)" = "#7570b3"),
                    guide = "none") +
  labs(x = NULL, y = "wells",
       title = "Is peak RFU trustworthy? Where the maximum falls in each series") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(fig.dir, "04_Fig_peak_RFU_capture.png"), fig4,
       width = 8, height = 4, dpi = 300)

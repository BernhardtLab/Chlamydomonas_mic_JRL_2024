# Description -------------------------------------------------------------
#
# 10_methods_fig_insets.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# September 29, 2026
#
# Two small inset panels for the methods (BioRender) figure, element 5:
#   (i)  one control well's RFU trajectory + fitted exponential (mu extraction)
#   (ii) the control (none) thermal performance curve (mu_max vs temp + fit)
# Faithful to 02 (same trim + mu/time) and to 04 (the saved control TPC fit).
#
# What this script produces:
#
# Figures-main:
#   04_fig_mu_plot.png            fitted curve over the raw growth timeseries data
#   05_fig_TPC_plot.png           TPC model fit for one (control probably)
#
# Load packages -----------------------------------------------------------

library(tidyverse)
library(bayestestR)

# Set paths ---------------------------------------------------------------

proc.dir <- "Data-processed"
mod.dir  <- "Models"
ind.t    <- file.path(mod.dir, "individual")
fig.dir  <- "Figures-main"
raw <- read.csv(file.path(proc.dir, "01_timeseries_data.csv"))
mu  <- read.csv(file.path(proc.dir, "02_chlamy_mus.csv"))

# Helper functions --------------------------------------------------------

lactin2 <- function(temp, a, b, tmax, dt) exp(a*temp) - exp(a*tmax - (tmax-temp)/dt) + b
bavg <- function(D, p) {
  ic <- D[, paste0("b_", p, "_Intercept")]
  bl <- grep(paste0("^b_", p, "_block"), colnames(D), value = TRUE)
  if (length(bl)) ic + rowSums(D[, bl, drop = FALSE]) / (length(bl) + 1) else ic
}

curve_draws <- function(fit, grid, pars, fun) {
  D <- as.matrix(fit); P <- lapply(pars, function(p) bavg(D, p))
  sapply(grid, function(x) do.call(fun, c(list(x), P)))
}

summ <- function(M, grid) {
  h <- apply(M, 2, function(x) { d <- hdi(x, ci = 0.95); c(d$CI_low, d$CI_high) })
  tibble(x = grid, med = apply(M, 2, median), lo = h[1, ], hi = h[2, ])
}
  
theme_inset <- theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = 9),
        axis.text  = element_text(size = 8),
        plot.margin = margin(4, 4, 4, 4))

# Inset 1: one control well, RFUs and fitted exponential ------------------

# pick a clean benign control well: positive growth, most reads in the fit window
well_id <- mu |>
  filter(mic == "none", nit == 1000, salt == 0, mu > 0.3) |>
  arrange(desc(n.points), desc(mu)) |>
  slice(1) |> pull(unique.id)
cat("growth-fit inset well:", well_id, "\n")

w <- raw |> filter(Chlamy.y.n == "y", unique.id == well_id) |> arrange(days)
w_fit <- if (nrow(w) > 1 && !is.na(w$RFU[2]) && w$RFU[2] < w$RFU[1]) w[-1, ] else w  # 02's trim
N0   <- w_fit$RFU[1]
r    <- mu$mu[mu$unique.id == well_id]     # 02's estimate
wend <- mu$time[mu$unique.id == well_id]   # 02's window endpoint (days)
d0   <- w_fit$days[1]

fit_curve <- tibble(days = seq(d0, wend, length.out = 100),
                    RFU  = N0 * exp(r * days))          # RFU = N0 * exp(r*days), N0 fixed

p_grow <- ggplot() +
  geom_point(data = w, aes(days, RFU), colour = "grey30", size = 1.4) +      # all reads
  geom_line(data = fit_curve, aes(days, RFU), colour = "black", linewidth = 1.2) +
  labs(x = "Time (days)", y = "Relative fluorescence units") +
  xlim(0,10)+
  theme_inset 

p_grow

# Inset 2: control (no microbe) TPC ---------------------------------------

fit_none_t <- readRDS(file.path(ind.t, "TPC_mic_none.rds"))
tpc.pars   <- c("a", "b", "tmax", "dt")
temp.grid  <- with(fit_none_t$data, seq(min(temp), max(temp), length.out = 200))
Fn_sum <- summ(curve_draws(fit_none_t, temp.grid, tpc.pars, lactin2), temp.grid)
topt   <- temp.grid[which.max(Fn_sum$med)]
raw_pts <- fit_none_t$data   

p_tpc <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_ribbon(data = Fn_sum, aes(x, ymin = lo, ymax = hi), fill = "grey70", alpha = 0.3) +
  geom_point(data = raw_pts, aes(temp, mu), colour = "grey30", size = 1.2, alpha = 0.55) +
  geom_line(data = Fn_sum, aes(x, med), colour = "black", linewidth = 0.9) +
  labs(x = "Temperature (\u00b0C)",
       y = expression(italic('\u03bc')[max]~"("*day^-1*")")) +
  theme_inset

p_tpc

# Save the figures --------------------------------------------------------

ggsave(file.path(fig.dir, "04_fig_mu_plot.png"), p_grow,
       width = 2, height = 2, dpi = 1000)

ggsave(file.path(fig.dir, "05_fig_TPC_plot.png"), p_tpc,
       width = 3, height = 3, dpi = 1000)


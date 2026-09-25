# Description -------------------------------------------------------------
#
# 09_effects_gradients.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# September 24, 2026
#
# This script will assess evidence for the stress-gradient hypothesis:
# does the benefit of microbial inoculation to C. reinhardtii vary as 
# a function of the intensity of abiotic stress? 
#
# We will measure the the interactions between benefits and stress by
# asssessing the difference in growth rates (mu) between inoculated and
# control (un-inoculated) C. reinhardtii. 
#
# Primary approach: structured-curve differences. We already fit each
# microbial treatment's gradient responses to temperature, nitrogen, and 
# salt variation (scripts 04-06, individual models). We'll use posterior
# draws (mu_treatment(X) - mu_nomicrobe(X)) over a grid of environmental 
# variables, identifying areas where the 95% HDI does not overlap 0. 
# We will also identify if and where relationships switch from beneficial 
# to harmful, and whether facilitation increases with stress. 
#
# Each axis will keep its raw environmental variable, not absolute distance.
# Distance from the optimum is anchored to the control's optimum.
#
# What this script produces:
#
# Figures-main:
#   03_fig3_SGH_effects.png                    benefits of microbial inoculation~gradients
#
# Load packages -----------------------------------------------------------

library(tidyverse)
library(brms)          # to load the fits + as.matrix() their draws (no refit)
library(bayestestR)    # hdi()
library(patchwork)

# Set paths ---------------------------------------------------------------

mod.dir      <- "Models"
fig.main.dir <- "Figures-main"

# Curve functions (same forms as 04-06) 
lactin2  <- function(temp, a, b, tmax, dt) exp(a * temp) - exp(a * tmax - (tmax - temp) / dt) + b
monod    <- function(nit, mumax, ks)       mumax * nit / (ks + nit)
expdecay <- function(salt, mmin, A, k)     mmin + A * exp(-k * salt)

# block-average a parameter's fixed part: intercept + mean(block coefs) if present
bavg <- function(D, p) {
  ic <- D[, paste0("b_", p, "_Intercept")]
  bl <- grep(paste0("^b_", p, "_block"), colnames(D), value = TRUE)
  if (length(bl)) ic + rowSums(D[, bl, drop = FALSE]) / (length(bl) + 1) else ic
}

# reconstruct a treatment's fitted curve over a grid, per posterior draw
#   fit  : an individual brms fit; pars : nlpar names; fun : the curve function
#   returns a (draws x length(grid)) matrix of block-averaged predictions
curve_draws <- function(fit, grid, pars, fun) {
  D <- as.matrix(fit)
  P <- lapply(pars, function(p) bavg(D, p))          # one draw-vector per parameter
  sapply(grid, function(x) do.call(fun, c(list(x), P)))
}

# interaction = treatment curve - control curve, per draw (independent fits -> align by index)
interaction_draws <- function(fit_trt, fit_ref, grid, pars, fun) {
  A <- curve_draws(fit_trt, grid, pars, fun)
  R <- curve_draws(fit_ref, grid, pars, fun)
  n <- min(nrow(A), nrow(R))
  A[seq_len(n), ] - R[seq_len(n), ]
}

# per-grid-point summary: median + 95% HDI + does the HDI clear 0?
summarise_int <- function(M, grid) {
  h <- apply(M, 2, function(x) { d <- hdi(x, ci = 0.95); c(d$CI_low, d$CI_high) })
  tibble(x = grid, med = apply(M, 2, median), lo = h[1, ], hi = h[2, ],
         sig = h[1, ] > 0 | h[2, ] < 0)
}

# Test run: temperature: community (all) interaction -------------------------------
tpc.pars <- c("a", "b", "tmax", "dt")
ind.t    <- file.path(mod.dir, "individual")            # individual fits
fit_none_t <- readRDS(file.path(ind.t, "TPC_mic_none.rds"))
fit_all  <- readRDS(file.path(ind.t, "TPC_mic_all.rds"))

temp.grid <- with(fit_none$data, seq(min(temp), max(temp), length.out = 1000))

# control's thermal optimum (peak of the no-microbe curve) -> reference line
Fn   <- curve_draws(fit_none, temp.grid, tpc.pars, lactin2)
topt <- temp.grid[which.max(apply(Fn, 2, median))]

int_all <- interaction_draws(fit_all, fit_none, temp.grid, tpc.pars, lactin2)
si_all  <- summarise_int(int_all, temp.grid)

cat(sprintf("Control Topt ~ %.1f C | interaction clears 0 over %.0f%% of the range\n",
            topt, 100 * mean(si_all$sig)))

# plot: community interaction band vs temperature
pT <- ggplot(si_all, aes(x, med)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = topt, linetype = "dotted", colour = "grey40", linewidth = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "darkgreen", alpha = 0.25) +
  geom_line(colour = "darkgreen", linewidth = 0.8) +
  geom_rug(data = filter(si_all, sig), aes(x = x), sides = "b",
           colour = "darkgreen", inherit.aes = FALSE) +          # where HDI clears 0
  labs(x = "Temperature (°C)",
       y = expression(atop("Microbial effect on growth",
                           mu[with] - mu[without] ~ (day^-1))),
       title = "A) Temperature") +
  theme_classic() +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))
pT

# Determine which microbes have significant effects across all 3 --------

## Temperature -------------------------------------------------------------

trts <- c("all", as.character(1:15))

int_t <- purrr::map_dfr(trts, function(m) {
  fit_m <- readRDS(file.path(ind.t, paste0("TPC_mic_", m, ".rds")))
  summarise_int(interaction_draws(fit_m, fit_none, temp.grid, tpc.pars, lactin2),
                temp.grid) |> mutate(mic = m)
})

# contiguous significant regions (where the 95% HDI clears 0)
step <- diff(temp.grid)[1]

sig_regions_t <- int_t |>
  filter(sig) |>
  mutate(dir = ifelse(med > 0, "facilitation", "antagonism")) |>
  arrange(mic, dir, x) |>
  group_by(mic, dir) |>
  mutate(run = cumsum(c(TRUE, diff(x) > 1.5 * step))) |>     # break into contiguous runs
  group_by(mic, dir, run) |>
  summarise(from = round(min(x), 1), to = round(max(x), 1),
            peak = round(med[which.max(abs(med))], 3), .groups = "drop") |>
  arrange(match(mic, trts), from) |>
  select(mic, dir, from, to, peak)

cat("\n== temperature: significant interaction regions (°C) ==\n")
print(as.data.frame(sig_regions_t), row.names = FALSE)

cat("\nmicrobes with any significant region:",
    paste(setdiff(unique(sig_regions_t$mic), "all"), collapse = ", "), "\n")

temp_sig <- setdiff(unique(sig_regions_t$mic), "all")   # temp

# consistent microbe colours (by ID, so they match across panels) + italic names
mic_cols <- c("2" = "goldenrod2", "11" = "brown4")   
sp <- c("1"="E. meliloti","2"="E. meliloti","3"="S. cerevisiae","4"="P. psychrotolerans",
        "5"="N. halotolerans","6"="P. macmurdoensis","7"="P. sulfinovorans","8"="B. aerius",
        "9"="C. braakii","10"="R. rosettiformans","11"="P. protegens","12"="S. pituitosa",
        "13"="R. qingshengii","14"="R. phycosphaerae","15"="P. putida")
lab_it <- function(b) parse(text = paste0("italic('", sp[b], "')~'(", b, ")'"))

# add community + "Others" to the palette
mic_cols[c("all", "Others")] <- c("mediumblue", "grey45") 

sig_any_t <- setdiff(unique(sig_regions_t$mic), "all")   # = c("2","11") for temperature

## Nitrogen ----------------------------------------------------------------

nit.pars <- c("mumax", "ks")

fit_none_n <- readRDS(file.path(ind.t, "nit_mic_none.rds"))
nit.grid <- with(fit_none$data, seq(min(nit), max(nit), length.out = 1000))

int_n <- purrr::map_dfr(trts, function(m) {
  fit_m <- readRDS(file.path(ind.t, paste0("nit_mic_", m, ".rds")))
  summarise_int(interaction_draws(fit_m, fit_none, nit.grid, nit.pars, monod),
                nit.grid) |> mutate(mic = m)
})

# contiguous significant regions (where the 95% HDI clears 0)
step <- diff(nit.grid)[1]

sig_regions_n <- int_n |>
  filter(sig) |>
  mutate(dir = ifelse(med > 0, "facilitation", "antagonism")) |>
  arrange(mic, dir, x) |>
  group_by(mic, dir) |>
  mutate(run = cumsum(c(TRUE, diff(x) > 1.5 * step))) |>
  group_by(mic, dir, run) |>
  summarise(from = round(min(x), 1), to = round(max(x), 1),
            peak = round(med[which.max(abs(med))], 3), .groups = "drop") |>
  arrange(match(mic, trts), from) |>
  select(mic, dir, from, to, peak)

cat("\n== nitrogen: significant interaction regions (µM) ==\n")
print(as.data.frame(sig_regions_n), row.names = FALSE)

cat("\nmicrobes with any significant region:",
    paste(setdiff(unique(sig_regions_n$mic), "all"), collapse = ", "), "\n")

sig_any_n <- setdiff(unique(sig_regions_n$mic), "all") # 11 and 13

mic_cols["13"] <- "plum3" # Add a colour for microbe 13

nit_sig  <- setdiff(unique(sig_regions_n$mic), "all")   # after nitrogen

## Salt --------------------------------------------------------------------

salt.pars <- c("mmin", "A", "k")

ind.s     <- file.path(mod.dir, "individual")     # salt fits live here
fit_none_s  <- readRDS(file.path(ind.s, "salt_mic_none.rds"))
salt.grid <- with(fit_none$data, seq(min(salt), max(salt), length.out = 1000))

int_s <- purrr::map_dfr(trts, function(m) {
  fit_m <- readRDS(file.path(ind.s, paste0("salt_mic_", m, ".rds")))
  summarise_int(interaction_draws(fit_m, fit_none, salt.grid, salt.pars, expdecay),
                salt.grid) |> mutate(mic = m)
})

# contiguous significant regions (where the 95% HDI clears 0)
step <- diff(salt.grid)[1]

sig_regions_s <- int_s |>
  filter(sig) |>
  mutate(dir = ifelse(med > 0, "facilitation", "antagonism")) |>
  arrange(mic, dir, x) |>
  group_by(mic, dir) |>
  mutate(run = cumsum(c(TRUE, diff(x) > 1.5 * step))) |>
  group_by(mic, dir, run) |>
  summarise(from = round(min(x), 2), to = round(max(x), 2),
            peak = round(med[which.max(abs(med))], 3), .groups = "drop") |>
  arrange(match(mic, trts), from) |>
  select(mic, dir, from, to, peak)

cat("\n== salt: significant interaction regions (g/L) ==\n")
print(as.data.frame(sig_regions_s), row.names = FALSE)

cat("\nmicrobes with any significant region:",
    paste(setdiff(unique(sig_regions_s$mic), "all"), collapse = ", "), "\n")

# Salt panel highlighting significant microbes

sig_any_s <- setdiff(unique(sig_regions_s$mic), "all") # 1, 2, 3, 4, 9, 11, 13

mic_cols["1"] <- "chocolate3" # Add a colour for microbe 1
mic_cols["3"] <- "skyblue" # Add a colour for microbe 3
mic_cols["4"] <- "olivedrab4" # Add a colour for microbe 4
mic_cols["9"] <- "darkorchid3" # Add a colour for microbe 9

salt_sig <- setdiff(unique(sig_regions_s$mic), "all")   # after salt

sig_union <- as.character(sort(as.integer(union(union(temp_sig, nit_sig), salt_sig))))

# Temperature ----------------

# curve layers
bg <- filter(int_t, !mic %in% c(sig_any_t, "all"))  # thin grey background medians
hl <- filter(int_t,  mic %in% sig_any_t)            # coloured significant medians

# significance bars: facilitation top / antagonism bottom 
yr <- range(c(si_all$lo, si_all$hi, int_t$med)); spn <- diff(yr)
place_t <- function(d, top) {
  d <- droplevels(mutate(d, mic = factor(mic, levels = c("all", sig_any_t))))
  d$row <- as.integer(d$mic)
  d$y   <- if (top) yr[2] + 0.03 * spn * d$row else yr[1] - 0.03 * spn * d$row
  d
}

fac <- place_t(filter(sig_regions_t, dir == "facilitation"), TRUE)
ant <- place_t(filter(sig_regions_t, dir == "antagonism"),  FALSE)

brks     <- c("all", sig_union, "Others")
labs_vec <- c(expression("Microbial community"),
              do.call(c, lapply(sig_union, lab_it)),
              expression("Others"))

## Panel D: microbial effects~temp -----------------------------------------

pD <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = topt, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
  geom_ribbon(data = si_all, aes(x, ymin = lo, ymax = hi), fill = "mediumblue", alpha = 0.2) +
  geom_line(data = bg,     aes(x, med, group = mic, colour = "Others"), linewidth = 0.35) +
  geom_line(data = hl,     aes(x, med, group = mic, colour = mic),      linewidth = 0.8) +
  geom_line(data = si_all, aes(x, med, colour = "all"),                 linewidth = 1) +
  # microbe bars — coloured but NOT in legend (the lines above already carry them)
  geom_segment(data = filter(fac, mic != "all"),
               aes(from, y, xend = to, yend = y, colour = mic),
               linewidth = 1.2, lineend = "butt", show.legend = FALSE) +
  geom_segment(data = filter(ant, mic != "all"),
               aes(from, y, xend = to, yend = y, colour = mic),
               linewidth = 1.2, lineend = "butt", show.legend = FALSE) +
  # community bars — fixed colour, no legend
  geom_segment(data = filter(fac, mic == "all"),
               aes(from, y, xend = to, yend = y), colour = "mediumblue",
               linewidth = 1.2, lineend = "butt") +
  geom_segment(data = filter(ant, mic == "all"),
               aes(from, y, xend = to, yend = y), colour = "mediumblue",
               linewidth = 1.2, lineend = "butt") +
  scale_colour_manual(values = mic_cols, breaks = brks, labels = labs_vec,
                      limits = brks, drop = FALSE, name = "Microbial treatment") +
  labs(x = "Temperature (°C)",
       y = expression("Microbial growth effect (" * day^-1 * ")"),
       title = "D) Temperature") +
  theme_classic() +
  theme(axis.title  = element_text(size = 10),
        axis.text   = element_text(size = 10),
        plot.title  = element_text(size = 10, face = "bold", hjust = 0.03),
        plot.margin  = margin(3, 3, 3, 3))
pD

## Panel A: control TPC ----------------------------------------------------

# control (none) TPC: median + 95% HDI over the grid
Fn     <- curve_draws(fit_none_t, temp.grid, tpc.pars, lactin2) # draws × grid
Fn_sum <- summarise_int(Fn, temp.grid)                          # x, med, lo, hi (+sig, unused)

topt   <- temp.grid[which.max(apply(Fn, 2, median))]

# raw control-treatment data (points): grab the response column name
yvar <- as.character(fit_none_t$formula$formula[[2]]) 
raw  <- fit_none_t$data

pA <- ggplot() +
  geom_vline(xintercept = topt, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_ribbon(data = Fn_sum, aes(x, ymin = lo, ymax = hi), fill = "grey40", alpha = 0.2) +
  geom_point(data = raw, aes(.data[["temp"]], .data[[yvar]]),
             colour = "grey35", size = 1, alpha = 0.55) +
  geom_line(data = Fn_sum, aes(x, med), colour = "black", linewidth = 1) +
  labs(x = "Temperature (°C)",
       y = expression("Growth rate (" * day^-1 * ")"),
       title = "A) Temperature") +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03),
        plot.margin = margin(3, 3, 3, 3))
pA

# Nitrogen ----------------------------------------------------------------

# curve layers
si_all <- filter(int_n, mic == "all")
bg <- filter(int_n, !mic %in% c(sig_any_n, "all"))
hl <- filter(int_n,  mic %in% sig_any_n)

# significance bars: facilitation top / antagonism bottom
yr <- range(c(si_all$lo, si_all$hi, int_n$med)); spn <- diff(yr)
place_n <- function(d, top) {
  d <- droplevels(mutate(d, mic = factor(mic, levels = c("all", sig_any_n))))
  d$row <- as.integer(d$mic)
  d$y   <- if (top) yr[2] + 0.03 * spn * d$row else yr[1] - 0.03 * spn * d$row
  d
}

fac <- place_n(filter(sig_regions_n, dir == "facilitation"), TRUE)
ant <- place_n(filter(sig_regions_n, dir == "antagonism"),  FALSE)

brks     <- c("all", sig_union, "Others")
labs_vec <- c(expression("Microbial community"),
              do.call(c, lapply(sig_union, lab_it)),
              expression("Others"))

## Panel E: microbial effects~nitrogen -------------------------------------

pE <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 1000, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
  geom_ribbon(data = si_all, aes(x, ymin = lo, ymax = hi), fill = "mediumblue", alpha = 0.2) +
  geom_line(data = bg,     aes(x, med, group = mic, colour = "Others"), linewidth = 0.35) +
  geom_line(data = hl,     aes(x, med, group = mic, colour = mic),      linewidth = 0.8) +
  geom_line(data = si_all, aes(x, med, colour = "all"),                 linewidth = 1) +
  geom_segment(data = filter(fac, mic != "all"),
               aes(from, y, xend = to, yend = y, colour = mic),
               linewidth = 1.2, lineend = "butt", show.legend = FALSE) +
  geom_segment(data = filter(ant, mic != "all"),
               aes(from, y, xend = to, yend = y, colour = mic),
               linewidth = 1.2, lineend = "butt", show.legend = FALSE) +
  geom_segment(data = filter(fac, mic == "all"),
               aes(from, y, xend = to, yend = y), colour = "mediumblue",
               linewidth = 1.2, lineend = "butt") +
  geom_segment(data = filter(ant, mic == "all"),
               aes(from, y, xend = to, yend = y), colour = "mediumblue",
               linewidth = 1.2, lineend = "butt") +
  scale_colour_manual(values = mic_cols, breaks = brks, labels = labs_vec,
                      limits = brks, drop = FALSE, name = "Microbial treatment") +
  labs(x = "Nitrogen (µM)",
       y = expression("Microbial growth effect (" * day^-1 * ")"),
       title = "E) Nitrogen") +
  theme_classic() +
  theme(axis.title  = element_text(size = 10),
        axis.text   = element_text(size = 10),
        plot.title  = element_text(size = 10, face = "bold", hjust = 0.03),
        plot.margin  = margin(3, 3, 3, 3))
pE

# Which microbes are almost significant at boosting growth under N-limitation?

lowN_idx <- which(nit.grid <= 250)

pd_lowN <- purrr::map_dfr(setdiff(trts, "all"), function(m) {
  M  <- interaction_draws(readRDS(file.path(ind.t, paste0("nit_mic_", m, ".rds"))),
                          fit_none_n, nit.grid, nit.pars, monod)[, lowN_idx, drop = FALSE]
  pd <- apply(M, 2, function(col) mean(col > 0))     # P(facilitation) at each low-N point
  tibble(mic = m,
         pd_max = round(max(pd), 3),
         N_at   = round(nit.grid[lowN_idx][which.max(pd)], 0),
         med_at = round(apply(M, 2, median)[which.max(pd)], 3))
}) |> arrange(desc(pd_max))

cat("\n== nitrogen: closeness to facilitation at low N (≤250 µM) ==\n")
print(as.data.frame(pd_lowN), row.names = FALSE)

## Panel B: control Monod ---------------------------------------------------

Fn     <- curve_draws(fit_none_n, nit.grid, nit.pars, monod)
Fn_sum <- summarise_int(Fn, nit.grid)

yvar <- as.character(fit_none_n$formula$formula[[2]])
raw  <- fit_none_n$data

pB <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 1000, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
  geom_ribbon(data = Fn_sum, aes(x, ymin = lo, ymax = hi), fill = "grey40", alpha = 0.2) +
  geom_point(data = raw, aes(.data[["nit"]], .data[[yvar]]),
             colour = "grey35", size = 1, alpha = 0.55) +
  geom_line(data = Fn_sum, aes(x, med), colour = "black", linewidth = 1) +
  labs(x = "Nitrogen (µM)",
       y = expression("Growth rate (" * day^-1 * ")"),
       title = "B) Nitrogen") +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03),
        plot.margin = margin(3, 3, 3, 3))
pB

# Salt --------------------------------------------------------------------

# curve layers
si_all <- filter(int_s, mic == "all")
bg <- filter(int_s, !mic %in% c(sig_any_s, "all"))
hl <- filter(int_s,  mic %in% sig_any_s)

# significance bars: facilitation top / antagonism bottom
yr <- range(c(si_all$lo, si_all$hi, int_s$med)); spn <- diff(yr)
place_s <- function(d, top) {
  d <- droplevels(mutate(d, mic = factor(mic, levels = c("all", sig_any_s))))
  d$row <- as.integer(d$mic)
  d$y   <- if (top) yr[2] + 0.03 * spn * d$row else yr[1] - 0.03 * spn * d$row
  d
}

fac <- place_s(filter(sig_regions_s, dir == "facilitation"), TRUE)
ant <- place_s(filter(sig_regions_s, dir == "antagonism"),  FALSE)

brks     <- c("all", sig_union, "Others")
labs_vec <- c(expression("Microbial community"),
              do.call(c, lapply(sig_union, lab_it)),
              expression("Others"))

## Panel F: microbial effects~salt -----------------------------------------

pF <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
  geom_ribbon(data = si_all, aes(x, ymin = lo, ymax = hi), fill = "mediumblue", alpha = 0.2) +
  geom_line(data = bg,     aes(x, med, group = mic, colour = "Others"), linewidth = 0.35) +
  geom_line(data = hl,     aes(x, med, group = mic, colour = mic),      linewidth = 0.8) +
  geom_line(data = si_all, aes(x, med, colour = "all"),                 linewidth = 1) +
  geom_segment(data = filter(fac, mic != "all"),
               aes(from, y, xend = to, yend = y, colour = mic),
               linewidth = 1.2, lineend = "butt", show.legend = FALSE) +
  geom_segment(data = filter(ant, mic != "all"),
               aes(from, y, xend = to, yend = y, colour = mic),
               linewidth = 1.2, lineend = "butt", show.legend = FALSE) +
  geom_segment(data = filter(fac, mic == "all"),
               aes(from, y, xend = to, yend = y), colour = "mediumblue",
               linewidth = 1.2, lineend = "butt") +
  geom_segment(data = filter(ant, mic == "all"),
               aes(from, y, xend = to, yend = y), colour = "mediumblue",
               linewidth = 1.2, lineend = "butt") +
  scale_colour_manual(values = mic_cols, breaks = brks, labels = labs_vec,
                      limits = brks, drop = FALSE, name = "Microbial treatment") +
  labs(x = expression("Salt (g L"^-1*")"),
       y = expression("Microbial growth effect (" * day^-1 * ")"),
       title = "F) Salt") +
  theme_classic() +
  theme(axis.title  = element_text(size = 10),
        axis.text   = element_text(size = 10),
        plot.title  = element_text(size = 10, face = "bold", hjust = 0.03),
        plot.margin  = margin(3, 3, 3, 3))
pF

## Panel C: control salt decay ---------------------------------------------

Fn     <- curve_draws(fit_none_s, salt.grid, salt.pars, expdecay)
Fn_sum <- summarise_int(Fn, salt.grid)

yvar <- as.character(fit_none_s$formula$formula[[2]])
raw  <- fit_none_s$data

pC <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_ribbon(data = Fn_sum, aes(x, ymin = lo, ymax = hi), fill = "grey40", alpha = 0.2) +
  geom_point(data = raw, aes(.data[["salt"]], .data[[yvar]]),
             colour = "grey35", size = 1, alpha = 0.55) +
  geom_line(data = Fn_sum, aes(x, med), colour = "black", linewidth = 1) +
  labs(x = expression("Salt (g L"^-1*")"),
       y = expression("Growth rate (" * day^-1 * ")"),
       title = "C) Salt") +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03),
        plot.margin = margin(3, 3, 3, 3))
pC

# Assemble the figure -----------------------------------------------------

pD2 <- pD + guides(colour = "none")
pE2 <- pE + guides(colour = "none")

design <- "
ABC#
DEFG
"
fig3 <- pA + pB + pC + pD2 + pE2 + pF + guide_area() +
  plot_layout(design = design, guides = "collect", widths = c(1, 1, 1, 0.5)) &
  theme(legend.key.size   = unit(0.7, "lines"),
        legend.key.width  = unit(1.1, "lines"),
        legend.text       = element_text(size = 7),
        legend.title      = element_text(size = 8),
        legend.background = element_blank(),
        legend.key        = element_blank())

fig3

ggsave(file.path(fig.main.dir, "03_fig3_SGH_effects.png"),
       fig3, width = 12, height = 7, dpi = 600)

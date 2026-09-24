# Description -------------------------------------------------------------
#
# 07_trait_effect_analysis.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# September 23, 2026
#
# Analysis of microbial effects on niche-determining traits in C. reinhardtii
#
# This script will address whether microbial treatments affect C. 
# reinhardtii's niche-determining traits relative to the no-microbe
# baseline AND whether community effects exceed the mean and sum of 
# individual microbial effects.
#
# Part 1: which traits are affected by microbial inoculation?
#
# Part 2: are there any cases in which emergent community effects drive
# increases in C. reinhardtii performance?
#
# Part 3: compile and export table of statistical results
#
# Part 4: generate Figure 1 
#
# What this script produces:
#
# Data-processed:
#   10_microbial_effects_stat_sum.csv      summary of statistical test of microbial effects
#
# Figures-main:
#   01_fig1_microbial_effects.png          microbial effects on 6 key traits (2 per gradient)
#
# Figures-misc:
#   10_fig_full_microbial_effects.png      all trait effects, plotted individually
#
# Load packages -----------------------------------------------------------

library(tidyverse)     # dplyr, ggplot2, purrr, tidyr
library(brms)          # reload pooled models for the per-draw community test (P2)
library(posterior)
library(bayestestR)    # hdi()
library(ggrepel)
library(patchwork)

# Set paths ---------------------------------------------------------------
proc.dir     <- "Data-processed"
mod.dir      <- "Models"
fig.main.dir <- "Figures-main"     # publication figure(s) -- NEW directory
fig.misc.dir <- "Figures-misc"     # supporting / full-trait figures

# inputs: per-axis microbe - none contrast tables (from 04 / 05 / 06)
eff.files <- c(
  temperature = file.path(proc.dir, "05_TPC_microbe_effects.csv"),
  nitrogen    = file.path(proc.dir, "07_nit_monod_microbe_effects.csv"),
  salt        = file.path(proc.dir, "09_salt_microbe_effects.csv"))

# inputs: pooled fitted models (for the per-draw emergent-community test in P2)
mod.files <- c(
  temperature = file.path(mod.dir, "01_TPC_brms.rds"),
  nitrogen    = file.path(mod.dir, "03_monod_nit_brms.rds"),
  salt        = file.path(mod.dir, "05_salt_brms.rds"))

# inputs: full posteriors for the individual models
tr.files <- c(
  temperature = file.path(mod.dir, "02_TPC_tr_ind.rds"),
  nitrogen    = file.path(mod.dir, "04_nit_monod_tr_ind.rds"),
  salt        = file.path(mod.dir, "06_salt_tr_ind.rds"))

traits.files <- c(
  temperature = file.path(proc.dir, "04_TPC_traits_bayes.csv"),
  nitrogen    = file.path(proc.dir, "06_nit_monod_traits_bayes.csv"),
  salt        = file.path(proc.dir, "08_salt_traits_bayes.csv"))

# Part 1: which traits are affected by microbial inoculation? --------------
# Load the three microbe - none contrast tables, tag axis + community/microbe,
# and summarise where the 95% HDI clears 0 -- this decides the trait shortlist
# and which model structure (pooled vs individual) to headline.

effects <- imap_dfr(eff.files, ~ read.csv(.x) |> mutate(axis = .y)) |>
  mutate(axis   = factor(axis, levels = names(eff.files)),
         method = factor(method, levels = c("pooled", "individual")),
         kind   = ifelse(mic == "all", "community", "microbe"))
# columns: mic, trait, delta, lo, hi, excl0, pNA, method, axis, kind

# What traits do microbial treatment affect? And which model is the most informative?
signif.count <- effects |>
  group_by(axis, trait, method) |>
  summarise(n_signif = sum(excl0), .groups = "drop") |>
  pivot_wider(names_from = method, values_from = n_signif)

cat("\n== treatments (of 16) clearing 0 vs no-microbe, per trait x model ==\n")
print(signif.count, n = Inf)

# Temp: Topt (2), Tbr and mumax (1 each)
# Nit: affinity and ks (3 each), mumax (1)
# Salt: salt tolerance (4) and limit (7)

# We will proceed with the individual model, seeing that they also fit
# the raw growth data better (especially for temperature)

# Part 2: evaluate evidence of emergent community effects --------------

tr <- lapply(tr.files, readRDS)     # axis -> (treatment -> draws x traits)

emergent_axis <- function(trA) {
  mics   <- setdiff(names(trA), c("none", "all"))
  traits <- colnames(trA[["all"]])
  n      <- min(vapply(trA, nrow, integer(1)))            # common draw count
  none   <- trA[["none"]][seq_len(n), , drop = FALSE]
  allc   <- trA[["all"]][seq_len(n), , drop = FALSE]
  purrr::map_dfr(traits, function(v) {
    mic <- sapply(mics, function(m) trA[[m]][seq_len(n), v])       # n x 15
    E.mean <- allc[, v] - rowMeans(mic, na.rm = TRUE)             # all - mean(15 traits)
    E.sum  <- (allc[, v] - none[, v]) - rowSums(mic - none[, v], na.rm = TRUE)
    hm <- bayestestR::hdi(E.mean[is.finite(E.mean)], ci = 0.95)
    hs <- bayestestR::hdi(E.sum[is.finite(E.sum)],  ci = 0.95)
    tibble(trait = v,
           d_mean = median(E.mean, na.rm = TRUE), lo_mean = hm$CI_low, hi_mean = hm$CI_high,
           sig_mean = hm$CI_low > 0 | hm$CI_high < 0,
           d_sum  = median(E.sum, na.rm = TRUE),  lo_sum = hs$CI_low,  hi_sum = hs$CI_high,
           sig_sum  = hs$CI_low > 0 | hs$CI_high < 0)
  })
}

emergent <- purrr::imap_dfr(tr, ~ emergent_axis(.x) |> mutate(axis = .y, .before = 1))

cat("\n== emergent community effects (all vs mean & sum of the 15) ==\n")
print(as.data.frame(emergent |> mutate(across(where(is.numeric), ~ round(.x, 3)))),
      row.names = FALSE)

# emergent effects (greater than the mean of microbial effects):
# Temp: mumax
# Nit: mumax, ks, and affinity
# Salt: salt limit

# Part 3: compile a statistical summary table

# (A) individual-model contrasts vs the no-microbe control
A <- effects |>
  filter(method == "individual") |>
  transmute(axis, trait,
            type     = ifelse(mic == "all", "community_vs_none", "microbe_vs_none"),
            contrast = ifelse(mic == "all", "all", paste0("mic_", mic)),
            median = delta, lo, hi, sig = excl0)

# (B) emergent effects, reshaped long (one row per null)
B <- emergent |>
  rename(median_mean = d_mean, median_sum = d_sum) |>
  pivot_longer(-c(axis, trait),
               names_to = c(".value", "null"),
               names_pattern = "(median|lo|hi|sig)_(mean|sum)") |>
  transmute(axis, trait,
            type     = paste0("emergent_", null),
            contrast = ifelse(null == "mean", "all - mean(15)", "all - sum(15)"),
            median, lo, hi, sig)

# combine, round on a per-value scale (signif handles ks~40 and affinity~0.01)
lev <- c("microbe_vs_none", "community_vs_none", "emergent_mean", "emergent_sum")
stat.tbl <- bind_rows(A, B) |>
  mutate(across(c(median, lo, hi), ~ signif(.x, 3)),
         type = factor(type, levels = lev)) |>
  arrange(axis, trait, type, contrast)

write.csv(stat.tbl, file.path(proc.dir, "10_microbial_effects_stat_sum.csv"),
          row.names = FALSE)

cat("\n== significant results (individual model) ==\n")
print(as.data.frame(filter(stat.tbl, sig)), row.names = FALSE)
cat("\nfull table:", nrow(stat.tbl), "rows -> 10_microbial_effects_stat_sum.csv\n")

# Part 4: Figure 1 -- trait-space effects, two panels per gradient ----------

# species lookup (numbers match the microbe IDs)
sp <- c("1"="Ensifer meliloti Em1021", "2"="Ensifer meliloti Em1022",
        "3"="Saccharomyces cerevisiae", "4"="Pseudarthrobacter psychrotolerans",
        "5"="Nesterenkonia halotolerans", "6"="Paenisporosarcina macmurdoensis",
        "7"="Pseudarthrobacter sulfinovorans", "8"="Bacillus aerius",
        "9"="Citrobacter braakii", "10"="Rhizobium rosettiformans",
        "11"="Pseudomonas protogens", "12"="Sphingomonas pituitosa",
        "13"="Rhodococcus qingshengii", "14"="Rhizorhabdus phycosphaerae",
        "15"="Pseudomonas putida")

mic_labs <- c(none = "Control (none)", all = "Microbial community", sp)

lev <- rev(c("none", "all", as.character(1:15)))   # none at top, 15 at bottom

# one panel's data: median + HDI + significance, per treatment, for one trait
fdat <- function(ax, tr_name) {
  m <- tr.sum  |> filter(axis == ax, trait == tr_name) |> select(mic, x = median, lo, hi)
  s <- sig.tab |> filter(axis == ax, trait == tr_name) |> select(mic, sig = excl0)
  left_join(m, s, by = "mic") |>
    mutate(sig   = coalesce(sig, FALSE),
           type  = case_when(mic == "none" ~ "none", mic == "all" ~ "community",
                             TRUE ~ "microbe"),
           alpha = ifelse(type == "none" | sig, 1, 0.3),          # control + sig = solid
           mic   = factor(mic, levels = lev))
}
exp_x   <- function(ax, tr) { A <- tr[[1]]; NULL }               # placeholder guard (see note)
expect  <- function(ax, tr_name) {
  trA <- get("tr")[[ax]]; mics <- setdiff(names(trA), c("none", "all"))
  n   <- min(vapply(trA, nrow, integer(1)))
  median(rowMeans(sapply(mics, function(m) trA[[m]][seq_len(n), tr_name]), na.rm = TRUE), na.rm = TRUE)
}
emg_sig <- function(ax, tr_name)
  isTRUE(emergent$sig_mean[emergent$axis == ax & emergent$trait == tr_name])

cols <- c(none = "red2", microbe = "black", community = "mediumblue")

forest_panel <- function(ax, tr_name, xlab, title, show_y = FALSE) {
  d      <- fdat(ax, tr_name)
  none_x <- d$x[d$mic == "none"]
  allrow <- filter(d, mic == "all")
  seg    <- data.frame(x = expect(ax, tr_name), xend = allrow$x); seg$mic <- allrow$mic
  es     <- emg_sig(ax, tr_name)
  
  ggplot(d, aes(x, mic, colour = type)) +
    geom_vline(xintercept = none_x, linetype = "dashed", linewidth = 0.4, colour = "grey20") +
    geom_segment(data = seg, aes(x = x, xend = xend, y = mic, yend = mic),
                 colour = "mediumblue", linewidth = if (es) 0.7 else 0.4,
                 linetype = if (es) "solid" else "dotted", inherit.aes = FALSE) +
    geom_point(data = seg, aes(x = x, y = mic), shape = 21, fill = NA,
               colour = "mediumblue", size = 2.6, stroke = 0.6, inherit.aes = FALSE) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, alpha = alpha), height = 0, linewidth = 0.7) +
    geom_point(aes(alpha = alpha), size = 2) +
    scale_alpha_identity() +
    scale_colour_manual(values = cols, breaks = c("none", "microbe", "community"),
                        labels = c("Control", "Single microbe", "Microbial community")) +
    scale_y_discrete(labels = mic_labs) +
    labs(x = xlab, y = NULL, colour = NULL, title = title) +
    theme_classic() +
    theme(axis.title  = element_text(size = 10),
          axis.text.x = element_text(size = 9),
          axis.text.y = if (show_y) element_text(size = 8, angle = 0, hjust = 1, face = "italic") else element_blank(),
          axis.ticks.y = if (show_y) element_line() else element_blank(),
          plot.title  = element_text(size = 10, face = "bold", hjust = 0.03))
}

xlab_mu  <- expression("Maximum growth rate" ~ italic("\u03bc")[italic(max)] ~ (day^-1))
pA <- forest_panel("temperature", "r.max",   xlab_mu, "A) Temperature", show_y = TRUE)
pB <- forest_panel("nitrogen",    "mumax",   xlab_mu, "B) Nitrogen")
pC <- forest_panel("salt",        "mu_max",  xlab_mu, "C) Salt")
pD <- forest_panel("temperature", "breadth",
                   expression("Thermal breadth" ~ italic("T")[italic(br)] ~ "(°C)"),
                   "D) Temperature", show_y = TRUE)
pE <- forest_panel("nitrogen", "affinity",
                   expression("Nitrogen affinity" ~ italic("\u03bc")[italic(max)] / italic("K")[italic(S)]),
                   "E) Nitrogen")
pF <- forest_panel("salt", "s_lim",
                   expression("Critical salinity" ~ italic("S")[italic(crit)] ~ (g ~ L^-1)),
                   "F) Salt")

fig1 <- (pA | pB | pC) / (pD | pE | pF) & theme(legend.position = "none")

ggsave(file.path(fig.main.dir, "01_fig1_microbial_effects.png"),
       fig1, width = 12, height = 9, dpi = 300)

# Supplementary figure: remaining traits not shown in Fig 1 
sA <- forest_panel("temperature", "Tmin",
                   expression("Thermal minimum" ~ italic("T")[italic(min)] ~ "(°C)"),
                   "A) Temperature", show_y = TRUE)
sB <- forest_panel("temperature", "Topt",
                   expression("Thermal optimum" ~ italic("T")[italic(opt)] ~ "(°C)"),
                   "B) Temperature")
sC <- forest_panel("temperature", "Tmax",
                   expression("Thermal maximum" ~ italic("T")[italic(max)] ~ "(°C)"),
                   "C) Temperature")
sD <- forest_panel("nitrogen", "ks",
                   expression("Half-saturation constant" ~ italic("K")[italic(S)] ~ (mu * M)),
                   "D) Nitrogen", show_y = TRUE)
sE <- forest_panel("salt", "salt_tol",
                   expression("Salt tolerance" ~ italic("S")[50] ~ (g ~ L^-1)),
                   "E) Salt")

fig_supp <- (sA | sB | sC) / (sD | sE | plot_spacer()) &
  theme(legend.position = "none")

ggsave(file.path(fig.misc.dir, "10_fig_full_microbial_effects.png"),
       fig_supp, width = 12, height = 9, dpi = 300)

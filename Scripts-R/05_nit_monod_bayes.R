# Description -------------------------------------------------------------
#
# 05_nit_monod_bayes.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# September 22, 2026
#
# Nitrogen-response (Monod) curves for algal mu, fit two ways:
#
# (1) hierarchical Bayesian model (brms) across all microbial treatments
# at once. Microbe is a partially-pooled effect on the Monod curve
# parameters; block is a fixed effect. The biologically meaningful traits
# (mumax, ks, affinity = mumax/ks) and the microbe-vs-control differences
# in them are derived from the joint posterior, per draw.
#
# (2) Bayesian models fit individually to each microbial treatment;
# block as a fixed effect.
#
# NOTE: block is a fixed effect because we do not have sufficient levels
# (3), and models failed to converge when it was a random effect (see
# 04_TPC_bayes.R). Its placement in the individual fits follows the same
# logic: it is retained only on the parameter it most affects in the
# pooled model.
#
# Monod:  mu = mumax * N / (ks + N)
#   mumax = maximum growth rate (asymptote at saturating N)
#   ks    = half-saturation constant (N at which mu = mumax / 2)
#
# What this script produces:
#
#   Models/02_nit_monod_brms.rds                 the fitted, pooled brms model
#   Models/nit_individual/nit_monod_mic_*.rds    the individual fits
#
# Data-processed:
#   06_nit_monod_traits_bayes.csv       per-microbe trait posteriors (med + HDI)
#   07_nit_monod_microbe_effects.csv    microbe - none trait contrasts
#
# Figures-misc:
#   06_fig_nit_monod_fits.png           fitted curves over the data for both models
#   07_fig_nit_monod_effects.png        microbe-vs-control effect posteriors


# Load packages -----------------------------------------------------------

# One-time setup (run once per machine, NOT on every source):
#   install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
#   cmdstanr::check_cmdstan_toolchain(fix = TRUE)
#   cmdstanr::install_cmdstan(cores = 2)
# We use the cmdstanr backend (compiles Stan models through CmdStan directly),
# which avoids the rstan/Windows "sink: invalid connection" error thrown by
# brms's default backend on this setup.

library(cmdstanr)
library(tidyverse)     # dplyr, ggplot2, purrr, tidyr
library(brms)
library(posterior)     # summarise_draws()
library(bayestestR)    # hdi()

# Set paths ---------------------------------------------------------------
proc.dir <- "Data-processed"
fig.dir  <- "Figures-misc"
mod.dir  <- "Models"
ind.dir  <- file.path(mod.dir, "individual")    # per-microbe fits go here

# inputs / outputs (single source of truth)
mu.file      <- "02_chlamy_mus.csv"                     # input: per-well mu (from 02)
pooled.file  <- file.path(mod.dir, "02_nit_monod_brms")       # brms appends .rds itself
traits.file  <- file.path(proc.dir, "06_nit_monod_traits_bayes.csv")
effects.file <- file.path(proc.dir, "07_nit_monod_microbe_effects.csv")

options(mc.cores = parallel::detectCores())

# Monod functions ---------------------------------------------------------

monod <- function(nit, mumax, ks) mumax * nit / (ks + nit)   # for prediction/plotting

# Derive the trait set from ONE draw's parameters.
derive_one <- function(mumax, ks) {
  c(mumax = mumax, ks = ks, affinity = mumax / ks)           # affinity = initial slope at N->0
}

# coverage: fraction of wells inside their 95% predictive interval (~0.95 = calibrated)
coverage <- function(fit, prob = 0.95) {
  pp <- posterior_predict(fit)
  y  <- fit$data$mu
  a  <- (1 - prob) / 2
  qi <- t(apply(pp, 2, quantile, c(a, 1 - a)))
  mean(y >= qi[, 1] & y <= qi[, 2])
}

# Load the data -----------------------------------------------------------
df <- read.csv(file.path(proc.dir, mu.file))

df.n <- df %>%
  filter(salt == 0, temp == 30,    # nitrogen gradient = reference T & salt
         block != 2,               # drop carbon-free block
         rep < 4,                   # drop the few rep-4 (pipetting-corrected) wells
         !is.na(mu)) %>%            # drop wells where mu couldn't be estimated
  mutate(mic   = relevel(factor(mic), ref = "none"),
         block = factor(block))

cat("Monod data:", nrow(df.n), "wells |",
    nlevels(df.n$mic), "microbe levels |",
    "N levels:", paste(sort(unique(df.n$nit)), collapse = ", "), "|",
    "blocks", paste(levels(df.n$block), collapse = ","), "\n")

## VERIFY: quick pooled nls to ballpark the priors (ignores structure).
## Use these MLEs to centre the logmumax / logks intercept priors below
## (remember the priors are on the log scale, so use log(mumax), log(ks)).
print(nls(mu ~ mumax * nit / (ks + nit), data = df.n,
          start = list(mumax = 1.7, ks = 50)))

# Modelling: hierarchical Monod -------------------------------------------

mon.formula <- bf(
  mu ~ mumax * nit / (ks + nit),
  mumax ~ 1 + block + (1 | mic),
  ks    ~ 1 + block + (1 | mic),
  sigma ~ poly(log1p(nit), 2),        # <- heteroscedasticity in N
  nl = TRUE
)  # both parameters vary by microbe, block fixed

## VERIFY the two intercept priors against the nls ballpark above:
##   mumax intercept ~ 1.6
##   ks    intercept ~ 56

mon.priors <- c(
  prior(normal(1.56, 0.5), nlpar = "mumax", coef = "Intercept"),
  prior(normal(0,    0.3), nlpar = "mumax"),                 # block contrasts
  prior(normal(56,   20),  nlpar = "ks",    coef = "Intercept"),
  prior(normal(0,    15),  nlpar = "ks"),                    # block contrasts
  prior(normal(-1.4, 0.5), class = "Intercept", dpar = "sigma"),  # log(0.235) ≈ -1.45
  prior(normal(0, 1),      class = "b",         dpar = "sigma"),  # poly terms
  prior(normal(0, 0.3),    class = "sd", nlpar = "mumax"),   # microbe-RE SD (rate units)
  prior(normal(0, 15),     class = "sd", nlpar = "ks")       # microbe-RE SD (ks units)
)

n.mic <- length(unique(df.n$mic))

init_fun <- function() list(
  b_mumax = as.array(c(1.56, 0, 0)),
  b_ks    = as.array(c(56,   0, 0)),
  sd_1 = as.array(0.10), z_1 = matrix(0, 1, n.mic),   # mumax RE
  sd_2 = as.array(8.0),  z_2 = matrix(0, 1, n.mic)    # ks RE
)

mon.fit <- brm(
  mon.formula, data = df.n, family = gaussian(),
  prior = mon.priors, init = init_fun,
  chains = 4, iter = 3000, warmup = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 1,
  file = pooled.file,
  file_refit = "on_change", backend = "cmdstanr"
)

# Diagnostics: pooled model -----------------------------------------------
ss <- summarise_draws(mon.fit, "rhat", "ess_bulk", "ess_tail")
np <- nuts_params(mon.fit)

# convergence (covers Rhat, ESS, divergences, treedepth in one)
c(max_rhat     = max(ss$rhat, na.rm = TRUE),
  min_ess_bulk = min(ss$ess_bulk, na.rm = TRUE),
  min_ess_tail = min(ss$ess_tail, na.rm = TRUE),
  divergences  = sum(np$Value[np$Parameter == "divergent__"]),
  treedepth    = sum(np$Value[np$Parameter == "treedepth__"] >= 12))

# fit / calibration / predictive accuracy
bayes_R2(mon.fit)
coverage(mon.fit)
loo_fit <- loo(mon.fit); print(loo_fit)
mean(loo_fit$diagnostics$pareto_k > 0.7)

# posterior-predictive checks + microbe variance
pp_check(mon.fit, ndraws = 100)
pp_check(mon.fit, type = "stat_grouped", stat = "sd", group = "nit")
print(VarCorr(mon.fit))

# Modelling: individual per-microbe Monod ---------------------------------

pooled <- readRDS(file.path(mod.dir, "02_nit_monod_brms.rds"))   # if not already loaded as mon.fit
bf <- fixef(pooled)
round(bf[grepl("block", rownames(bf)), c("Estimate", "Q2.5", "Q97.5")], 3)

f_ind <- bf(
  mu ~ mumax * nit / (ks + nit),
  mumax ~ 1 + block,          # block3 is real (~10%); cheap to keep
  ks    ~ 1 + block,          # block effects are huge—mandatory 
  sigma ~ poly(log1p(nit), 2),
  nl = TRUE
)

p_ind <- c(
  prior(normal(1.56, 0.5), nlpar = "mumax", coef = "Intercept"),
  prior(normal(0,    0.3), nlpar = "mumax"),                    # block contrasts
  prior(normal(56,   20),  nlpar = "ks",    coef = "Intercept"),
  prior(normal(0,    15),  nlpar = "ks"),                      # block contrasts
  prior(normal(-1.4, 0.5), class = "Intercept", dpar = "sigma"),
  prior(normal(0, 1),      class = "b",         dpar = "sigma")
)

fit_one <- function(m) {
  d  <- droplevels(subset(df.n, mic == m))
  nb <- nlevels(d$block)
  init_ind <- function() list(
    b_mumax = as.array(c(1.56, rep(0, nb - 1))),
    b_ks    = as.array(c(56,   rep(0, nb - 1)))
  )
  
  brm(f_ind, data = d, family = gaussian(), prior = p_ind, init = init_ind,
      chains = 4, iter = 3000, warmup = 1000,
      control = list(adapt_delta = 0.99, max_treedepth = 12), seed = 1,
      file = file.path(ind.dir, paste0("nit_mic_", m)),
      file_refit = "on_change", backend = "cmdstanr")
}

mics.chr <- as.character(sort(unique(df.n$mic)))
fits.ind <- setNames(lapply(mics.chr, fit_one), mics.chr)

# convergence scan across all fits
ind.diag <- purrr::map_dfr(mics.chr, function(m) {
  f <- fits.ind[[m]]; np <- brms::nuts_params(f)
  tibble(mic         = m,
         max_rhat    = round(max(brms::rhat(f), na.rm = TRUE), 3),
         divergences = sum(np$Value[np$Parameter == "divergent__"]),
         td_hits     = sum(np$Value[np$Parameter == "treedepth__"] >= 12))
})

print(ind.diag, n = Inf)

# Diagnostics: individual fits --------------------------------------------

diag_one <- function(fit) {
  ss <- summarise_draws(fit, "rhat", "ess_bulk", "ess_tail")
  np <- nuts_params(fit); lo <- suppressWarnings(loo(fit))
  tibble(max_rhat     = max(ss$rhat, na.rm = TRUE),
         min_ess_bulk = min(ss$ess_bulk, na.rm = TRUE),
         min_ess_tail = min(ss$ess_tail, na.rm = TRUE),
         divergences  = sum(np$Value[np$Parameter == "divergent__"]),
         treedepth    = sum(np$Value[np$Parameter == "treedepth__"] >= 12),
         bayes_R2     = as.numeric(bayes_R2(fit)[1, "Estimate"]),
         coverage     = coverage(fit),
         pareto_bad   = mean(lo$diagnostics$pareto_k > 0.7))
}

ind.full.diag <- purrr::map_dfr(mics.chr, function(m)
  diag_one(fits.ind[[m]]) |> dplyr::mutate(mic = m, .before = 1))

print(dplyr::mutate(ind.full.diag, dplyr::across(where(is.numeric), ~round(.x, 3))), n = Inf)

# Derive traits + microbe-vs-none contrasts -------------------------------

set.seed(1); NSUB <- 8000          
traits.set <- c("mumax", "ks", "affinity")
traits_of  <- function(post) t(mapply(derive_one, post$mumax, post$ks))

# block-average a parameter's fixed part: intercept + mean(block coefs) if present
bavg <- function(D, p) {
  ic <- D[, paste0("b_", p, "_Intercept")]
  bl <- grep(paste0("^b_", p, "_block"), colnames(D), value = TRUE)
  if (length(bl)) ic + rowSums(D[, bl, drop = FALSE]) / (length(bl) + 1) else ic
}

draws_pooled <- function(D, m) data.frame(
  mumax = bavg(D, "mumax") + D[, sprintf("r_mic__mumax[%s,Intercept]", m)],
  ks    = bavg(D, "ks")    + D[, sprintf("r_mic__ks[%s,Intercept]",    m)])
draws_ind <- function(D) data.frame(mumax = bavg(D, "mumax"), ks = bavg(D, "ks"))

# per-draw trait matrices (pooled: SAME draws across microbes -> paired contrasts)
Dp    <- as.matrix(mon.fit)
idx.c <- if (is.infinite(NSUB)) seq_len(nrow(Dp)) else sample(nrow(Dp), NSUB)
tr_pooled <- setNames(lapply(mics.chr, function(m)
  traits_of(draws_pooled(Dp, m)[idx.c, ])), mics.chr)
tr_ind <- setNames(lapply(mics.chr, function(m) {
  post <- draws_ind(as.matrix(fits.ind[[m]]))
  if (nrow(post) > NSUB) post <- post[sample(nrow(post), NSUB), ]
  traits_of(post)
}), mics.chr)

# trait table -> 06_nit_monod_traits_bayes.csv
summ_traits <- function(tr) purrr::map_dfr(colnames(tr), function(v) {
  x <- tr[, v]; ok <- x[is.finite(x)]; h <- bayestestR::hdi(ok, ci = 0.95)
  tibble(trait = v, median = median(ok), lo = h$CI_low, hi = h$CI_high,
         pNA = mean(!is.finite(x)))
})

traits <- bind_rows(
  purrr::map_dfr(mics.chr, ~ summ_traits(tr_pooled[[.x]]) |> mutate(mic = .x, method = "pooled")),
  purrr::map_dfr(mics.chr, ~ summ_traits(tr_ind[[.x]])    |> mutate(mic = .x, method = "individual"))
)

write.csv(traits, traits.file, row.names = FALSE)

# microbe v none contrasts -> 07_nit_monod_microbe_effects.csv
contrast_tbl <- function(tr_list) {
  none <- tr_list[["none"]]
  purrr::map_dfr(setdiff(names(tr_list), "none"), function(m) {
    mm <- tr_list[[m]]
    purrr::map_dfr(traits.set, function(v) {
      d <- mm[, v] - none[, v]; ok <- d[is.finite(d)]
      h <- bayestestR::hdi(ok, ci = 0.95)
      tibble(mic = m, trait = v, delta = median(ok),
             lo = h$CI_low, hi = h$CI_high,
             excl0 = h$CI_low > 0 | h$CI_high < 0,     # 95% HDI excludes 0
             pNA = mean(!is.finite(d)))
    })
  })
}

effects <- bind_rows(
  contrast_tbl(tr_pooled) |> mutate(method = "pooled"),
  contrast_tbl(tr_ind)    |> mutate(method = "individual"))

write.csv(effects, effects.file, row.names = FALSE)

# headline: # microbes clearing 0 per trait x method
effects |> group_by(trait, method) |> summarise(n = sum(excl0), .groups = "drop") |>
  tidyr::pivot_wider(names_from = method, values_from = n) |> print()

# Figures -----------------------------------------------------------------

# fitted curves over the raw data: pooled vs individual, per microbe
nit.grid <- seq(0, max(df.n$nit), length.out = 120)

curve_from <- function(post) {                 # post: data.frame(mumax, ks)
  M  <- sapply(nit.grid, function(N) monod(N, post$mumax, post$ks))
  hp <- apply(M, 2, function(x) bayestestR::hdi(x, ci = .95) |> (\(h) c(h$CI_low, h$CI_high))())
  tibble(nit = nit.grid, med = apply(M, 2, median), lo = hp[1, ], hi = hp[2, ])
}

fit.curves <- bind_rows(
  purrr::map_dfr(mics.chr, ~ curve_from(draws_pooled(Dp, .x)[idx.c, ]) |>
                   mutate(mic = .x, method = "pooled")),
  purrr::map_dfr(mics.chr, ~ curve_from(draws_ind(as.matrix(fits.ind[[.x]]))) |>
                   mutate(mic = .x, method = "individual"))
) |> mutate(mic = factor(mic, levels = c("none", "all", as.character(1:15))))

fig.fits <- ggplot() +
  geom_point(data = mutate(df.n, mic = factor(mic, levels = levels(fit.curves$mic))),
             aes(nit, mu), alpha = 0.18, size = 0.5) +
  geom_ribbon(data = fit.curves, aes(nit, ymin = lo, ymax = hi, fill = method), alpha = 0.22) +
  geom_line(data = fit.curves, aes(nit, med, colour = method), linewidth = 0.6) +
  facet_wrap(~ mic, ncol = 5) +
  scale_colour_manual(values = c(pooled = "#1b9e77", individual = "#d95f02")) +
  scale_fill_manual(  values = c(pooled = "#1b9e77", individual = "#d95f02")) +
  labs(x = "Nitrogen (uM)", y = "Growth rate (mu)",
       title = "Monod fits over raw data: pooled vs individual") +
  theme_bw() + theme(legend.position = "top")

ggsave(file.path(fig.dir, "06_fig_nit_fits.png"), fig.fits, width = 13, height = 9, dpi = 300)

# microbe - none effects: forest plot per trait, both methods
fig.eff <- effects |>
  mutate(mic = factor(mic, levels = rev(c("all", as.character(1:15))))) |>
  ggplot(aes(delta, mic, colour = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5),
                  fatten = 1.5) +
  facet_wrap(~ trait, scales = "free_x") +
  scale_colour_manual(values = c(pooled = "#1b9e77", individual = "#d95f02")) +
  labs(x = "microbe - none (95% HDI)", y = "Microbe",
       title = "Nitrogen-response effects of microbial treatment") +
  theme_bw()

ggsave(file.path(fig.dir, "07_fig_nit_effects.png"), fig.eff, width = 10, height = 6, dpi = 300)

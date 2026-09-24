# Description -------------------------------------------------------------
#
# 06_salt_expdecay_bayes.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# September 23, 2026
#
# Salt-tolerance curves for algal mu, fit two ways:
#
# (1) hierarchical Bayesian model (brms) across all microbial treatments
# at once. Microbe is a partially-pooled effect on the curve parameters;
# block is a fixed effect. The biologically meaningful traits (max growth,
# salt tolerance = salt at half-max growth, and the lethal zero-crossing)
# and the microbe-vs-control differences in them are derived from the joint
# posterior, per draw.
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
# CURVE FORM: floored exponential decay
#   mu = mmin + A * exp(-k * salt)
#     mmin = lower asymptote (growth at high salt; may be < 0 = net decline)
#     A    = amplitude of the decline
#     k    = decay rate (steepness)
# Growth declines from salt = 0 with no low-salt plateau, so a logistic's
# upper asymptote / inflection are not identifiable here (they extrapolate
# to negative salt). The floored exponential decay fits every microbe as
# well as the logistic (checked: RSS within 7% for all 17, better for the
# few where the logistic would not converge), with three interpretable,
# well-identified parameters and no funnel-prone degenerate terms.
#
# TRAITS (all derived from the fitted curve, per draw):
#   mu_max   = mmin + A                       growth at salt = 0
#   salt_tol = log(2A / (A - mmin)) / k       salt at half of mu_max
#   s_lim    = log(-A / mmin) / k             salt where mu = 0 (lethal limit;
#                                             defined only when mmin < 0)
# What this script produces:
#
#   Models/05_salt_brms.rds                    the fitted, pooled brms model
#         /06_salt_tr_ind.rds                  the full posterior for all microbes
#   Models/salt_individual/salt_mic_*.rds      the individual fits
#
# Data-processed:
#   08_salt_traits_bayes.csv      per-microbe trait posteriors (med + HDI)
#   09_salt_microbe_effects.csv   microbe - none trait contrasts
#
# Figures-misc:
#   08_fig_salt_fits.png          fitted curves over the data for both models
#   09_fig_salt_effects.png       microbe-vs-control effect posteriors
#
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
ind.dir  <- file.path(mod.dir, "individual")   # per-microbe fits go here

# inputs
mu.file      <- "02_chlamy_mus.csv"                    # input: per-well mu (from 02)
pooled.file  <- file.path(mod.dir, "05_salt_brms")     # brms appends .rds itself
traits.file  <- file.path(proc.dir, "08_salt_traits_bayes.csv")
effects.file <- file.path(proc.dir, "09_salt_microbe_effects.csv")

options(mc.cores = parallel::detectCores())

# Salt-curve functions ----------------------------------------------------

expdecay <- function(salt, mmin, A, k) mmin + A * exp(-k * salt)   # for prediction/plotting

# Derive the trait set from ONE draw's parameters. Closed forms (algebra,
# no root-finding); guarded to return NA where a trait is undefined.
derive_one <- function(mmin, A, k) {
  mu_max   <- mmin + A                                     # growth at salt = 0
  ok       <- A > 0 & k > 0 & mu_max > 0                   # curve declines from a positive max
  salt_tol <- if (ok) log(2 * A / (A - mmin)) / k else NA_real_   # salt at mu = mu_max / 2
  s_lim    <- if (ok & mmin < 0) log(-A / mmin) / k else NA_real_ # salt at mu = 0 (needs mmin < 0)
  c(mu_max = mu_max, salt_tol = salt_tol, s_lim = s_lim)
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

df.s <- df %>%
  filter(nit == 1000, temp == 30,   # salt gradient = reference T & N
         block != 2,                 # drop carbon-free block
         rep < 4,                    # drop the few rep-4 (pipetting-corrected) wells
         !is.na(mu)) %>%             # drop wells where mu couldn't be estimated
  mutate(mic   = relevel(factor(mic), ref = "none"),
         block = droplevels(factor(block)))

cat("Salt data:", nrow(df.s), "wells |",
    nlevels(df.s$mic), "microbe levels |",
    "salt levels:", paste(sort(unique(df.s$salt)), collapse = ", "), "|",
    "blocks", paste(levels(df.s$block), collapse = ","), "\n")

## VERIFY: quick pooled nls to ballpark the priors (ignores structure).
## Use these MLEs to centre the mmin / A / k intercept priors below.
print(nls(mu ~ mmin + A * exp(-k * salt), data = df.s,
          start = list(mmin = -0.1, A = 1.6, k = 0.35)))

# Modelling: hierarchical salt decay --------------------------------------

salt.formula <- bf(
  mu ~ mmin + A * exp(-k * salt),
  mmin ~ 1 + block + (1 | mic),
  A    ~ 1 + block + (1 | mic),
  k    ~ 1 + block + (1 | mic),
  sigma ~ poly(salt, 2),        # quadratic
  nl = TRUE
) # all three parameters vary by microbe, block fixed

salt.priors <- c(
  prior(normal(-0.13, 0.15), nlpar = "mmin", coef = "Intercept"),
  prior(normal(0,     0.10), nlpar = "mmin"),                  # block contrasts
  prior(normal(1.65,  0.40), nlpar = "A",    coef = "Intercept"),
  prior(normal(0,     0.30), nlpar = "A"),                     # block contrasts
  prior(normal(0.31,  0.15), nlpar = "k",    coef = "Intercept"),
  prior(normal(0,     0.10), nlpar = "k"),                     # block contrasts
  prior(normal(-1.9, 0.5), class = "Intercept", dpar = "sigma"),
  prior(normal(0, 1),      class = "b",         dpar = "sigma"),
  prior(normal(0, 0.10), class = "sd", nlpar = "mmin"),        # microbe-RE SDs
  prior(normal(0, 0.30), class = "sd", nlpar = "A"),
  prior(normal(0, 0.10), class = "sd", nlpar = "k")
)

n.mic <- nlevels(df.s$mic)
nb    <- ncol(model.matrix(~ 1 + block, data = df.s))   # block coefs (intercept + contrasts)

init_fun <- function() list(
  b_mmin = as.array(c(-0.13, rep(0, nb - 1))),
  b_A    = as.array(c(1.65,  rep(0, nb - 1))),
  b_k    = as.array(c(0.31,  rep(0, nb - 1))),
  sd_1 = as.array(0.05), z_1 = matrix(0, 1, n.mic),   # mmin RE
  sd_2 = as.array(0.15), z_2 = matrix(0, 1, n.mic),   # A RE
  sd_3 = as.array(0.05), z_3 = matrix(0, 1, n.mic)    # k RE
)

salt.fit <- brm(
  salt.formula, data = df.s, family = gaussian(),
  prior = salt.priors, init = init_fun,
  chains = 4, iter = 3000, warmup = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 1,
  file = pooled.file,
  file_refit = "on_change", backend = "cmdstanr"
)

# Diagnostics: pooled model -----------------------------------------------
ss <- summarise_draws(salt.fit, "rhat", "ess_bulk", "ess_tail")
np <- nuts_params(salt.fit)

# convergence (covers Rhat, ESS, divergences, treedepth in one)
c(max_rhat     = max(ss$rhat, na.rm = TRUE),
  min_ess_bulk = min(ss$ess_bulk, na.rm = TRUE),
  min_ess_tail = min(ss$ess_tail, na.rm = TRUE),
  divergences  = sum(np$Value[np$Parameter == "divergent__"]),
  treedepth    = sum(np$Value[np$Parameter == "treedepth__"] >= 12))

# fit / calibration / predictive accuracy
bayes_R2(salt.fit)
coverage(salt.fit)
loo_fit <- loo(salt.fit); print(loo_fit)
mean(loo_fit$diagnostics$pareto_k > 0.7)

# posterior-predictive checks + microbe variance
## CHECK the stat_grouped panels: if within-salt spread is badly mis-captured,
## model sigma on salt (e.g. sigma ~ poly(salt, 2)), as in the TPC / Monod.
pp_check(salt.fit, ndraws = 100)
pp_check(salt.fit, type = "stat_grouped", stat = "sd", group = "salt")
print(VarCorr(salt.fit))

# Modelling: individual per-microbe salt decay ----------------------------

pooled <- readRDS(file.path(mod.dir, "05_salt_brms.rds"))   # if not already loaded as salt.fit
bf <- fixef(pooled)
round(bf[grepl("block", rownames(bf)), c("Estimate", "Q2.5", "Q97.5")], 3)

f_ind <- bf(
  mu ~ mmin + A * exp(-k * salt),
  mmin ~ 1 + block,      # all three carry a real block effect (see fixef)
  A    ~ 1 + block,
  k    ~ 1 + block,
  sigma ~ poly(salt, 2), # same heteroscedasticity as the pooled model
  nl = TRUE
)

p_ind <- c(
  prior(normal(-0.13, 0.15), nlpar = "mmin", coef = "Intercept"),
  prior(normal(0,     0.10), nlpar = "mmin"),
  prior(normal(1.65,  0.40), nlpar = "A",    coef = "Intercept"),
  prior(normal(0,     0.30), nlpar = "A"),
  prior(normal(0.31,  0.15), nlpar = "k",    coef = "Intercept"),
  prior(normal(0,     0.10), nlpar = "k"),
  prior(normal(-1.9, 0.5), class = "Intercept", dpar = "sigma"),
  prior(normal(0, 1),      class = "b",         dpar = "sigma")
)


fit_one <- function(m) {
  d   <- droplevels(subset(df.s, mic == m))
  nbi <- ncol(model.matrix(~ 1 + block, data = d))   # block coefs in this microbe
  init_ind <- function() list(
    b_mmin = as.array(c(-0.13, rep(0, nbi - 1))),
    b_A    = as.array(c(1.65,  rep(0, nbi - 1))),
    b_k    = as.array(c(0.31,  rep(0, nbi - 1)))
    # no scalar sigma init: sigma is modelled (poly), brms defaults are fine
  )
  
  brm(f_ind, data = d, family = gaussian(), prior = p_ind, init = init_ind,
      chains = 4, iter = 3000, warmup = 1000,
      control = list(adapt_delta = 0.99, max_treedepth = 12), seed = 1,
      file = file.path(ind.dir, paste0("salt_mic_", m)),
      file_refit = "on_change", backend = "cmdstanr")
}

mics.chr <- as.character(sort(unique(df.s$mic)))
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
traits.set <- c("mu_max", "salt_tol", "s_lim")
traits_of  <- function(post) t(mapply(derive_one, post$mmin, post$A, post$k))

# block-average a parameter's fixed part: intercept + mean(block coefs) if present
bavg <- function(D, p) {
  ic <- D[, paste0("b_", p, "_Intercept")]
  bl <- grep(paste0("^b_", p, "_block"), colnames(D), value = TRUE)
  if (length(bl)) ic + rowSums(D[, bl, drop = FALSE]) / (length(bl) + 1) else ic
}

draws_pooled <- function(D, m) data.frame(
  mmin = bavg(D, "mmin") + D[, sprintf("r_mic__mmin[%s,Intercept]", m)],
  A    = bavg(D, "A")    + D[, sprintf("r_mic__A[%s,Intercept]",    m)],
  k    = bavg(D, "k")    + D[, sprintf("r_mic__k[%s,Intercept]",    m)])

draws_ind <- function(D) data.frame(mmin = bavg(D, "mmin"), A = bavg(D, "A"), k = bavg(D, "k"))

# per-draw trait matrices (pooled: SAME draws across microbes -> paired contrasts)
Dp    <- as.matrix(salt.fit)
idx.c <- if (is.infinite(NSUB)) seq_len(nrow(Dp)) else sample(nrow(Dp), NSUB)
tr_pooled <- setNames(lapply(mics.chr, function(m)
  traits_of(draws_pooled(Dp, m)[idx.c, ])), mics.chr)

tr_ind <- setNames(lapply(mics.chr, function(m) {
  post <- draws_ind(as.matrix(fits.ind[[m]]))
  if (nrow(post) > NSUB) post <- post[sample(nrow(post), NSUB), ]
  traits_of(post)
}), mics.chr)

saveRDS(tr_ind, file.path(mod.dir, "06_salt_tr_ind.rds"))

# trait table -> 08_salt_traits_bayes.csv
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

# microbe v none contrasts -> 09_salt_microbe_effects.csv
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
salt.grid <- seq(0, max(df.s$salt), length.out = 120)

curve_from <- function(post) {                 # post: data.frame(mmin, A, k)
  M  <- sapply(salt.grid, function(S) expdecay(S, post$mmin, post$A, post$k))
  hp <- apply(M, 2, function(x) bayestestR::hdi(x, ci = .95) |> (\(h) c(h$CI_low, h$CI_high))())
  tibble(salt = salt.grid, med = apply(M, 2, median), lo = hp[1, ], hi = hp[2, ])
}

fit.curves <- bind_rows(
  purrr::map_dfr(mics.chr, ~ curve_from(draws_pooled(Dp, .x)[idx.c, ]) |>
                   mutate(mic = .x, method = "pooled")),
  purrr::map_dfr(mics.chr, ~ curve_from(draws_ind(as.matrix(fits.ind[[.x]]))) |>
                   mutate(mic = .x, method = "individual"))
) |> mutate(mic = factor(mic, levels = c("none", "all", as.character(1:15))))

fig.fits <- ggplot() +
  geom_point(data = mutate(df.s, mic = factor(mic, levels = levels(fit.curves$mic))),
             aes(salt, mu), alpha = 0.18, size = 0.5) +
  geom_ribbon(data = fit.curves, aes(salt, ymin = lo, ymax = hi, fill = method), alpha = 0.22) +
  geom_line(data = fit.curves, aes(salt, med, colour = method), linewidth = 0.6) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  facet_wrap(~ mic, ncol = 5) +
  scale_colour_manual(values = c(pooled = "#1b9e77", individual = "#d95f02")) +
  scale_fill_manual(  values = c(pooled = "#1b9e77", individual = "#d95f02")) +
  labs(x = "Salt (g/L)", y = "Growth rate (mu)",
       title = "Salt-decay fits over raw data: pooled vs individual") +
  theme_bw() + theme(legend.position = "top")

ggsave(file.path(fig.dir, "08_fig_salt_fits.png"), fig.fits, width = 13, height = 9, dpi = 300)

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
       title = "Salt-response effects of microbial treatment") +
  theme_bw()

ggsave(file.path(fig.dir, "09_fig_salt_effects.png"), fig.eff, width = 10, height = 6, dpi = 300)

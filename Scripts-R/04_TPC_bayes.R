# Description -------------------------------------------------------------
#
# 04_TPC_bayes.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# September 15, 2026
#
# Thermal performance curves for algal mu, fit two ways: 
#
# (1) hierarchical Bayesian model (brms) across all microbial treatments 
# at once. Microbe is a partially-pooled effect on the Lactin-II curve 
# parameters; block is a fixed effect. The biologically meaningful traits 
# (Topt, Tmin/Tmax, thermal breadth, r.max) and the microbe-vs-control
# differences in them are derived from the joint posterior, per draw.
#
# (2) Bayesian models fit individually to each microbial treatment;
# block as a fixed effect. 
#
# NOTES: we are treating block as a fixed effect because we do not have
# sufficient levels (3), and models failed to converge when it was 
# treated as a random effect
#
# In the individual microbe models, we incorporate block only for the b
# parameter (which controls the vertical offset). Based on our pooled model,
# it explained the most variation (b-18%, a-7%, dt-3%, tmax-1%). We did this
# because there are not sufficient data points to estimate block effects
# for all TPC terms.
#
# Residual variance is modeled as temperature-dependent in the pooled and 
# individual models (sigma ~ poly(temp,2)). Growth-rate scatter was smallest at 
# the thermal extremes and largest at the rising/falling flanks of the curve.
#
# What this script produces:
#
#   Models/01_TPC_brms.rds                         the fitted, pooled brms model 
#         /02_TPC_tr_ind.rds                       the full posterior for all microbes
#   Models/individual/TPC_mic_1-15/none/all.rds    the individual fits
#
# Data-processed:   
#   04_TPC_traits_bayes.csv        per-microbe trait posteriors (med + HDI)
#   05_TPC_microbe_effects.csv     microbe - none trait contrasts
#
# Figures-misc: 
#   04_fig_TPC_fits.png            fitted curves over the data for both models
#   05_fig_TPC_trait_effects.png       microbe-vs-control effect posteriors
#
# Load packages -----------------------------------------------------------

# One-time setup (run once per machine, NOT on every source):
#   install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
#   cmdstanr::check_cmdstan_toolchain(fix = TRUE)
#   cmdstanr::install_cmdstan(cores = 2)

library(cmdstanr)
library(tidyverse)
library(brms)
library(posterior)     # as_draws_df()
library(bayestestR)    # hdi()

# Set paths ---------------------------------------------------------------
proc.dir <- "Data-processed"
fig.dir  <- "Figures-misc"
mod.dir  <- "Models"
mu.file  <- "02_chlamy_mus.csv"
ind.dir  <- file.path(mod.dir, "individual")     # per-microbe fits go here

pooled.file  <- file.path(mod.dir, "01_TPC_brms")      # brms appends .rds itself
traits.file  <- file.path(proc.dir, "04_TPC_traits_bayes.csv")
effects.file <- file.path(proc.dir, "05_TPC_microbe_effects.csv")

options(mc.cores = parallel::detectCores())

# Lactin II functions -----------------------------------------------------

lactin2 <- function(temp, a, b, tmax, d.t)
  exp(a * temp) - exp(a * tmax - (tmax - temp) / d.t) + b

lactin2_deriv <- function(temp, a, b, tmax, d.t)          # dRate/dT (for Topt)
  a * exp(a * temp) - (1 / d.t) * exp(a * tmax - (tmax - temp) / d.t)

lactin2_half <- function(temp, a, b, tmax, d.t, rh) # for FWHM breadth
  lactin2(temp, a, b, tmax, d.t) - rh

# Topt = where the slope is zero (root of the derivative).
calc_Topt <- function(a, b, tmax, d.t)
  tryCatch(uniroot(function(t) lactin2_deriv(t, a, b, tmax, d.t),
                   interval = c(-10, 45))$root, error = function(e) NA_real_)

# Generic safe root-finder for a crossing between `lo` and `hi`.
root_safe <- function(f, lo, hi) {
  flo <- f(lo); fhi <- f(hi)
  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) return(NA_real_)
  tryCatch(uniroot(f, interval = c(lo, hi))$root, error = function(e) NA_real_)
}

# Derive the full trait set from ONE draw's parameters.
derive_one <- function(a, b, tmax, d.t) {
  Topt  <- calc_Topt(a, b, tmax, d.t)
  rmax  <- if (is.na(Topt)) NA_real_ else lactin2(Topt, a, b, tmax, d.t)
  Tmin <- if (is.na(Topt)) NA_real_ else
    root_safe(function(t) lactin2(t, a, b, tmax, d.t), -50, Topt)
  Tmax <- if (is.na(Topt)) NA_real_ else
    root_safe(function(t) lactin2(t, a, b, tmax, d.t), Topt, 50)
  rh <- rmax / 2
  Tlo   <- if (is.na(Topt)) NA_real_ else
    root_safe(function(t) lactin2_half(t, a, b, tmax, d.t, rh), -10, Topt)
  Thi   <- if (is.na(Topt)) NA_real_ else
    root_safe(function(t) lactin2_half(t, a, b, tmax, d.t, rh), Topt, 50)
  c(Topt = Topt, r.max = rmax, Tmin = Tmin, Tmax = Tmax,
    breadth = Thi - Tlo)
}

# Load the data -----------------------------------------------------------

df <- read.csv(file.path(proc.dir, mu.file))

df.t <- df %>%
  filter(salt == 0, nit == 1000,   # thermal gradient = reference N & salt
         block != 2,               # drop carbon-free block
         rep < 4,                  # drop the few rep-4 (pipetting-corrected) wells
         !is.na(mu)) %>%           # drop wells where mu couldn't be estimated
  mutate(mic   = relevel(factor(mic), ref = "none"),
         block = factor(block))

cat("TPC data:", nrow(df.t), "wells |",
    nlevels(df.t$mic), "microbe levels |",
    "blocks", paste(levels(df.t$block), collapse = ","), "\n")

# Modelling: hierarchical Lactin II ---------------------------------------

tpc.formula <- bf(
  mu ~ exp(a * temp) - exp(a * tmax - (tmax - temp) / dt) + b,
  a     ~ 1 + block + (1 | mic),
  b     ~ 1 + block + (1 | mic),
  tmax  ~ 1 + block + (1 | mic),
  dt    ~ 1 + block + (1 | mic),
  sigma ~ poly(temp, 2),          # log-quadratic residual SD over temp (a hump)
  nl = TRUE
)

tpc.priors <- c(
  prior(normal(0,    0.03), nlpar = "a"),                     # controls the rise
  prior(normal(0.07, 0.03), nlpar = "a",    coef = "Intercept"),

  prior(normal(0,    1   ), nlpar = "b"),                     # the vertical offset
  prior(normal(-1,   1.5 ), nlpar = "b",    coef = "Intercept"),
  
  prior(normal(0,    3   ), nlpar = "tmax"),                  # maximum temp   
  prior(normal(44,   4   ), nlpar = "tmax", coef = "Intercept"),
  
  prior(normal(0,    3   ), nlpar = "dt"),                    # controls the rate of decline past Topt
  prior(normal(10,   3   ), nlpar = "dt",   coef = "Intercept"),
  
  prior(normal(0, 1), class = "b", dpar = "sigma"),           # residual varies with temperature
  
  prior(normal(0, 0.05), class = "sd", nlpar = "a"),           # microbial effects on model parameters
  prior(normal(0, 0.5 ), class = "sd", nlpar = "b"),
  prior(normal(0, 1.5 ), class = "sd", nlpar = "tmax"),
  prior(normal(0, 2   ), class = "sd", nlpar = "dt")
)

n.mic <- length(unique(df.t$mic))

init_fun <- function() list(
  b_a    = as.array(c(0.07,  0, 0)),
  b_b    = as.array(c(-0.94, 0, 0)),
  b_tmax = as.array(c(44,    0, 0)),
  b_dt   = as.array(c(10,    0, 0)),
  sd_1 = as.array(0.01), z_1 = matrix(0, 1, n.mic),
  sd_2 = as.array(0.20), z_2 = matrix(0, 1, n.mic),
  sd_3 = as.array(1.00), z_3 = matrix(0, 1, n.mic),
  sd_4 = as.array(1.00), z_4 = matrix(0, 1, n.mic)
)

tpc.fit <- brm(
  tpc.formula, data = df.t, family = gaussian(),
  prior = tpc.priors, init = init_fun,
  chains = 4, iter = 3000, warmup = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 1,
  file = file.path(mod.dir, "01_TPC_brms"),
  file_refit = "on_change", backend = "cmdstanr"
)

## Check model convergence and fit -----------------------------------------

ss <- summarise_draws(tpc.fit, "rhat", "ess_bulk", "ess_tail")
np <- nuts_params(tpc.fit)

# convergence (covers Rhat, ESS, divergences, treedepth in one)
c(max_rhat     = max(ss$rhat, na.rm = TRUE),
  min_ess_bulk = min(ss$ess_bulk, na.rm = TRUE),
  min_ess_tail = min(ss$ess_tail, na.rm = TRUE),
  divergences  = sum(np$Value[np$Parameter == "divergent__"]),
  treedepth    = sum(np$Value[np$Parameter == "treedepth__"] >= 12))

coverage <- function(fit, prob = 0.95) {
  pp <- posterior_predict(fit)                 # draws x observations
  y  <- fit$data$mu
  a  <- (1 - prob) / 2
  qi <- t(apply(pp, 2, quantile, c(a, 1 - a))) # per-obs [lo, hi]
  mean(y >= qi[, 1] & y <= qi[, 2])
} # calculate coverage

# fit / calibration / predictive accuracy
bayes_R2(tpc.fit)
coverage(tpc.fit)
loo_fit <- loo(tpc.fit); print(loo_fit)
mean(loo_fit$diagnostics$pareto_k > 0.7)

# posterior-predictive checks + microbe variance
pp_check(tpc.fit, ndraws = 100)
pp_check(tpc.fit, type = "stat_grouped", stat = "sd", group = "temp")
print(VarCorr(tpc.fit))
 
# Individual microbe curves -----------------------------------------------

pooled <- readRDS(file.path(mod.dir, "01_TPC_brms.rds"))   # if not already loaded as tpc.fit
bf <- fixef(pooled)
round(bf[grepl("block", rownames(bf)), c("Estimate", "Q2.5", "Q97.5")], 3)

f_ind <- bf(
  mu ~ exp(a * temp) - exp(a * tmax - (tmax - temp) / dt) + b,
  a    ~ 1,   
  b    ~ 1 + block,       
  tmax ~ 1,   
  dt   ~ 1, 
  sigma ~ poly(temp, 2),          # log-quadratic residual SD over temp (a hump)
  nl = TRUE
) # model

p_ind <- c(
  prior(normal(0.07, 0.03), nlpar = "a",  lb = 0, ub = 0.25),
  prior(normal(0,    1   ), nlpar = "b"),
  prior(normal(-1,   1.5 ), nlpar = "b", coef = "Intercept"),
  prior(normal(44,   4   ), nlpar = "tmax"),
  prior(normal(10,   3   ), nlpar = "dt", lb = 0),
  prior(normal(0, 1), class = "b", dpar = "sigma")
)

fit_one <- function(m) {
  d  <- droplevels(subset(df.t, mic == m))
  nb <- nlevels(d$block)                         # blocks present (usually 3)
  init_ind <- function() list(
    b_a = as.array(0.07),
    b_b = as.array(c(-0.9, rep(0, nb - 1))),
    b_tmax = as.array(44),
    b_dt = as.array(10)
  )
  
  brm(f_ind, data = d, family = gaussian(), prior = p_ind, init = init_ind,
      chains = 4, iter = 3000, warmup = 1000,
      control = list(adapt_delta = 0.99, max_treedepth = 12), seed = 1,
      file = file.path(ind.dir, paste0("TPC_mic_", m)),
      file_refit = "on_change", backend = "cmdstanr")
} # filter to each microbe and fit

mics     <- sort(unique(df.t$mic))
fits.ind <- setNames(lapply(mics, fit_one), mics)

ind.diag <- purrr::map_dfr(mics, function(m) {
  f  <- fits.ind[[m]]
  np <- brms::nuts_params(f)
  tibble(mic         = m,
         max_rhat    = round(max(brms::rhat(f), na.rm = TRUE), 3),
         divergences = sum(np$Value[np$Parameter == "divergent__"]),
         td_hits     = sum(np$Value[np$Parameter == "treedepth__"] >= 12))
}) # fit all 17 models. 

print(ind.diag, n = Inf)

## Check models' convergence and fits --------------------------------------

diag_one <- function(fit) {
  ss <- summarise_draws(fit, "rhat", "ess_bulk", "ess_tail")
  np <- nuts_params(fit)
  lo <- suppressWarnings(loo(fit))
  tibble(
    max_rhat     = max(ss$rhat, na.rm = TRUE),
    min_ess_bulk = min(ss$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(ss$ess_tail, na.rm = TRUE),
    divergences  = sum(np$Value[np$Parameter == "divergent__"]),
    treedepth    = sum(np$Value[np$Parameter == "treedepth__"] >= 12),
    bayes_R2     = as.numeric(bayes_R2(fit)[1, "Estimate"]),
    coverage     = coverage(fit),
    pareto_bad   = mean(lo$diagnostics$pareto_k > 0.7)
  )
}

ind.full.diag <- purrr::map_dfr(mics, function(m)
  diag_one(fits.ind[[m]]) |> dplyr::mutate(mic = m, .before = 1))

print(dplyr::mutate(ind.full.diag, dplyr::across(where(is.numeric), ~round(.x, 3))),
      n = Inf)

# Derive thermal trait posteriors and compute microbe-v-none constrasts ---------------------------------------------------

mics.chr <- as.character(mics)
set.seed(1); NSUB <- 8000 
traits.set <- c("Topt", "r.max", "Tmin", "Tmax", "breadth")
traits_of  <- function(post) t(mapply(derive_one, post$a, post$b, post$tmax, post$dt))

# block-average a parameter's fixed part: intercept + mean(block coefs) if present
bavg <- function(D, p) {
  ic <- D[, paste0("b_", p, "_Intercept")]
  bl <- grep(paste0("^b_", p, "_block"), colnames(D), value = TRUE)
  if (length(bl)) ic + rowSums(D[, bl, drop = FALSE]) / (length(bl) + 1) else ic
}

# pooled: block-averaged fixed part + microbe random deviation
draws_pooled <- function(D, m) data.frame(
  a    = bavg(D, "a")    + D[, sprintf("r_mic__a[%s,Intercept]",    m)],
  b    = bavg(D, "b")    + D[, sprintf("r_mic__b[%s,Intercept]",    m)],
  tmax = bavg(D, "tmax") + D[, sprintf("r_mic__tmax[%s,Intercept]", m)],
  dt   = bavg(D, "dt")   + D[, sprintf("r_mic__dt[%s,Intercept]",   m)])

# individual: fixed part only (block sits on b, so bavg returns intercepts for a/tmax/dt)
draws_ind <- function(D) data.frame(
  a = bavg(D, "a"), b = bavg(D, "b"), tmax = bavg(D, "tmax"), dt = bavg(D, "dt"))

# per-draw trait matrices
Dp    <- as.matrix(pooled)
idx.c <- if (is.infinite(NSUB)) seq_len(nrow(Dp)) else sample(nrow(Dp), NSUB)

# pooled: the SAME draws for every microbe -> paired (within-joint-posterior) contrasts
tr_pooled <- setNames(lapply(mics.chr, function(m)
  traits_of(draws_pooled(Dp, m)[idx.c, ])), mics.chr)

# individual: independent draws per fit
tr_ind <- setNames(lapply(mics.chr, function(m) {
  post <- draws_ind(as.matrix(fits.ind[[m]]))
  if (nrow(post) > NSUB) post <- post[sample(nrow(post), NSUB), ]
  traits_of(post)
}), mics.chr)

saveRDS(tr_ind, file.path(mod.dir, "02_TPC_tr_ind.rds"))

# trait table: 04_TPC_traits_bayes.csv
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

# microbe v none contrasts: 05_TPC_microbe_effects.csv
contrast_tbl <- function(tr_list) {
  none <- tr_list[["none"]]
  purrr::map_dfr(setdiff(names(tr_list), "none"), function(m) {
    mm <- tr_list[[m]]
    purrr::map_dfr(traits.set, function(v) {
      d <- mm[, v] - none[, v]; ok <- d[is.finite(d)]
      h <- bayestestR::hdi(ok, ci = 0.95)
      tibble(mic = m, trait = v, delta = median(ok),
             lo = h$CI_low, hi = h$CI_high,
             excl0 = h$CI_low > 0 | h$CI_high < 0,   # 95% HDI excludes 0
             pNA = mean(!is.finite(d)))
    })
  })
}

effects <- bind_rows(
  contrast_tbl(tr_pooled) |> mutate(method = "pooled"),
  contrast_tbl(tr_ind)    |> mutate(method = "individual")
)

write.csv(effects, effects.file, row.names = FALSE)

# Which microbes clear 0?
effects |> group_by(trait, method) |> summarise(n = sum(excl0), .groups = "drop") |>
  tidyr::pivot_wider(names_from = method, values_from = n) |> print()

# Fitting models with R2jags to double check outputs ----------------------
# Just to generate a comparison using an approach I previously relied upon, I will fit the individual models using R2jags

library(R2jags)

jags.file <- file.path("lactin_block.txt")

temp.grid <- seq(min(df.t$temp), max(df.t$temp), length.out = 100)
mics.chr  <- as.character(sort(unique(df.t$mic)))

# function to fit R2jags models to each microbe
run_jags_one <- function(m) {
  d <- droplevels(subset(df.t, mic == m))
  d$blk <- as.integer(factor(d$block))                 # 1..N.block integer codes
  jd <- list(N.obs = nrow(d), temp = d$temp, trait = d$mu,
             block = d$blk, N.block = max(d$blk),
             Temp.xs = temp.grid, N.Temp.xs = length(temp.grid))
  inits <- function() list(a = 0.07, tmax = 44, d.t = 10,
                           b = rep(-0.9, jd$N.block), sigma = 0.3)
  jags(data = jd, inits = inits,
       parameters.to.save = c("a","tmax","d.t","b","sigma","r.pred"),
       model.file = jags.file, n.chains = 4,
       n.iter = 22000, n.burnin = 2000, n.thin = 10)
} 

# fit the models
jags.curves <- purrr::map_dfr(mics.chr, function(m) {
  fit <- run_jags_one(m)
  sm  <- fit$BUGSoutput$summary
  rows <- grep("^r\\.pred\\[", rownames(sm))
  idx  <- as.integer(gsub(".*\\[(\\d+)\\].*", "\\1", rownames(sm)[rows]))
  rp   <- sm[rows, c("mean","2.5%","97.5%")][order(idx), ]
  tibble(mic = m, temp = temp.grid, med = rp[,"mean"], lo = rp[,"2.5%"], hi = rp[,"97.5%"],
         max_rhat = max(sm[, "Rhat"], na.rm = TRUE))
})

cat("JAGS max Rhat across microbes:", round(max(jags.curves$max_rhat), 3), "\n")

# draw the brms posteriors to compile curves with errors. 
brms.curves <- purrr::map_dfr(mics.chr, function(m) {
  post <- draws_ind(as.matrix(fits.ind[[m]]))
  M <- sapply(temp.grid, function(Te) lactin2(Te, post$a, post$b, post$tmax, post$dt))
  tibble(mic = m, temp = temp.grid, med = apply(M, 2, median),
         lo = apply(M, 2, quantile, .025), hi = apply(M, 2, quantile, .975))
})

# plot both models to compare fits
ggplot() +
  geom_point(data = mutate(df.t, mic = as.character(mic)), aes(temp, mu),
             alpha = 0.15, size = 0.5) +
  geom_ribbon(data = brms.curves, aes(temp, ymin = lo, ymax = hi, fill = "brms"), alpha = 0.15) +
  geom_ribbon(data = jags.curves, aes(temp, ymin = lo, ymax = hi, fill = "JAGS"), alpha = 0.15) +
  geom_line(data = brms.curves, aes(temp, med, colour = "brms"), linewidth = 0.6) +
  geom_line(data = jags.curves, aes(temp, med, colour = "JAGS"), linewidth = 0.6) +
  facet_wrap(~ mic, ncol = 5) +
  scale_colour_manual(values = c(brms = "#1b9e77", JAGS = "#d95f02"), name = "engine") +
  scale_fill_manual(  values = c(brms = "#1b9e77", JAGS = "#d95f02"), name = "engine") +
  labs(x = "Temperature (°C)", y = "Growth rate (mu)") +
  theme_bw() + theme(legend.position = "top")

# Figures -----------------------------------------------------------------

temp.grid <- seq(min(df.t$temp), max(df.t$temp), length.out = 120)

# (1) fitted curves over the raw data: pooled vs individual, per microbe
curve_from <- function(post) {
  M  <- sapply(temp.grid, function(Te) lactin2(Te, post$a, post$b, post$tmax, post$dt))
  hp <- apply(M, 2, function(x) { h <- bayestestR::hdi(x, ci = .95); c(h$CI_low, h$CI_high) })
  tibble(temp = temp.grid, med = apply(M, 2, median), lo = hp[1, ], hi = hp[2, ])
}

fit.curves <- bind_rows(
  purrr::map_dfr(mics.chr, ~ curve_from(draws_pooled(Dp, .x)[idx.c, ]) |>
                   mutate(mic = .x, method = "pooled")),
  purrr::map_dfr(mics.chr, ~ curve_from(draws_ind(as.matrix(fits.ind[[.x]]))) |>
                   mutate(mic = .x, method = "individual"))
) |> mutate(mic = factor(mic, levels = c("none", "all", as.character(1:15))))

fig.fits <- ggplot() +
  geom_point(data = mutate(df.t, mic = factor(mic, levels = levels(fit.curves$mic))),
             aes(temp, mu), alpha = 0.18, size = 0.5) +
  geom_ribbon(data = fit.curves, aes(temp, ymin = lo, ymax = hi, fill = method), alpha = 0.22) +
  geom_line(data = fit.curves, aes(temp, med, colour = method), linewidth = 0.6) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  facet_wrap(~ mic, ncol = 5) +
  scale_colour_manual(values = c(pooled = "#1b9e77", individual = "#d95f02")) +
  scale_fill_manual(  values = c(pooled = "#1b9e77", individual = "#d95f02")) +
  labs(x = "Temperature (°C)", y = "Growth rate (mu)",
       title = "Lactin-II fits over raw data: pooled vs individual") +
  theme_bw() + theme(legend.position = "top")

ggsave(file.path(fig.dir, "04_fig_TPC_fits.png"), fig.fits, width = 13, height = 9, dpi = 300)

# (2) microbe - none effect posteriors: forest plot per trait, both methods
fig.eff <- effects |>
  mutate(mic = factor(mic, levels = rev(c("all", as.character(1:15))))) |>
  ggplot(aes(delta, mic, colour = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.6), fatten = 1.3) +
  facet_wrap(~ trait, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c(pooled = "#1b9e77", individual = "#d95f02")) +
  labs(x = "microbe - none (95% HDI)", y = "Microbe",
       title = "Thermal-trait effects of microbial treatment (vs. none)") +
  theme_bw() + theme(legend.position = "top")

ggsave(file.path(fig.dir, "05_fig_TPC_trait_effects.png"), fig.eff, width = 12, height = 5, dpi = 300)

# Description -------------------------------------------------------------
#
# 12_peak_abd.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# October 8, 2026
#
# An exploration of how peak abundance (peak RFU count) varies across
# temperature, nitrogen, and salt gradients.
#
# We then fit Lactin II TPCs,
#
# What this script produces:
#
# Figures-main:
#   07_concept_fig.png            conceptual figure - microbial effects on TPC and interactions with stress gradients
#   
# Load packages -----------------------------------------------------------

library(cmdstanr)
library(tidyverse)
library(brms)
library(posterior)     # as_draws_df()
library(bayestestR)    # hdi()
library(patchwork)
library(ggrepel)

# Set paths ---------------------------------------------------------------

proc.dir <- "Data-processed"
fig.dir  <- "Figures-misc"

mod.dir  <- "Models"
ind.dir  <- file.path(mod.dir, "individual")     # per-microbe fits go here

mu <- read.csv(file.path(proc.dir, "02_chlamy_mus.csv"))

mu <- mu %>% filter(carbon.added == TRUE) # remove block 2

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

# Exploration -------------------------------------------------------------

slice_axes <- function(d) bind_rows(
  d %>% filter(nit == 1000, salt == 0)  %>% transmute(axis = "Temperature (C)", level = temp, peak.RFU, mic),
  d %>% filter(temp == 30,  salt == 0)  %>% transmute(axis = "Nitrogen (uM)",   level = nit,  peak.RFU, mic),
  d %>% filter(temp == 30,  nit == 1000) %>% transmute(axis = "Salt (g/L)",     level = salt, peak.RFU, mic)
)
ax.levels <- c("Temperature (C)", "Nitrogen (uM)", "Salt (g/L)")

none.ab <- slice_axes(filter(mu, mic == "none")) %>% mutate(axis = factor(axis, ax.levels))
all.ab  <- slice_axes(mu)                        %>% mutate(axis = factor(axis, ax.levels))

none.sum <- none.ab %>% group_by(axis, level) %>%
  summarise(mean = mean(peak.RFU), se = sd(peak.RFU) / sqrt(n()), n = n(), .groups = "drop")

f1 <- ggplot(none.ab, aes(level, peak.RFU)) +
  geom_point(alpha = 0.25, colour = "grey45") +
  geom_errorbar(data = none.sum, aes(level, ymin = mean - se, ymax = mean + se),
                inherit.aes = FALSE, width = 0, colour = "#1b9e77") +
  geom_line(data = none.sum, aes(level, mean), colour = "#1b9e77") +
  geom_point(data = none.sum, aes(level, mean), colour = "#1b9e77") +
  facet_wrap(~ axis, scales = "free_x") +
  labs(x = "gradient level", y = "peak RFU (alga alone)",
       title = "Peak algal abundance across gradients (alga alone)") +
  theme_bw()

f1

all.sum <- all.ab %>% group_by(axis, level, mic) %>%
  summarise(mean = mean(peak.RFU), .groups = "drop")
f2 <- ggplot(all.sum, aes(level, mean, group = mic)) +
  geom_line(colour = "grey60", alpha = 0.7) +
  geom_line(data = filter(all.sum, mic == "none"), colour = "black", linewidth = 1) +
  facet_wrap(~ axis, scales = "free_x") +
  labs(x = "gradient level", y = "mean peak RFU",
       title = "Peak abundance across gradients by microbial treatment (control in black)") +
  theme_bw()

f2

# Fit Lactin II TPCs ------------------------------------------------------

df.t <- mu %>%
  filter(salt == 0, nit == 1000, block != 2, rep < 4) %>%     # same subset as 04
  mutate(block = factor(block), y = peak.RFU / 1000)           # k-RFU

f_ind <- bf(
  y ~ exp(a * temp) - exp(a * tmax - (tmax - temp) / dt) + b,   # <- y, not mu
  a    ~ 1,
  b    ~ 1 + block,
  tmax ~ 1,
  dt   ~ 1,
  sigma ~ poly(temp, 2),
  nl = TRUE
) # model

p_ind <- c(
  prior(normal(0.09, 0.04), nlpar = "a",  lb = 0, ub = 0.25),
  prior(normal(0,    1.5 ), nlpar = "b",  coef = "Intercept"),   # k-RFU offset
  prior(normal(0,    1   ), nlpar = "b"),                        # block deviations (k-RFU)
  prior(normal(44,   3   ), nlpar = "tmax"),
  prior(normal(9,    3   ), nlpar = "dt", lb = 0),
  prior(normal(-1,   1   ), class = "Intercept", dpar = "sigma"),# log(resid SD) ~ log(0.3)
  prior(normal(0,    1   ), class = "b",         dpar = "sigma")
)

init_ind <- function() list(
  b_a    = as.array(0.09),
  b_b    = as.array(c(0, rep(0, nb - 1))),
  b_tmax = as.array(44),
  b_dt   = as.array(9)
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
      file = file.path(ind.dir, paste0("TPC_peak_abd_mic_", m)),
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

coverage <- function(fit, prob = 0.95) {
  pp <- posterior_predict(fit)                 # draws x observations
  y  <- fit$data$y
  a  <- (1 - prob) / 2
  qi <- t(apply(pp, 2, quantile, c(a, 1 - a))) # per-obs [lo, hi]
  mean(y >= qi[, 1] & y <= qi[, 2])
} # calculate coverage

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

## Plot the fits -----------------------------------------------------------

bp <- function(D, p) {                                # pop intercept + block-fixed average
  ic <- D[, paste0("b_", p, "_Intercept")]
  bl <- grep(paste0("^b_", p, "_block"), colnames(D), value = TRUE)
  if (length(bl)) ic + rowSums(D[, bl, drop = FALSE]) / (length(bl) + 1) else ic
}

lev <- c("none", as.character(1:15), "all")
tg  <- seq(5, 44, 0.25)

fit_band <- function(m) {
  D  <- as.matrix(fits.ind[[m]])
  a  <- bp(D,"a"); b <- bp(D,"b"); tm <- bp(D,"tmax"); dt <- bp(D,"dt")
  M  <- vapply(tg, function(t) 1000 * lactin2(t, a, b, tm, dt), numeric(length(a)))
  tibble(mic = m, temp = tg,
         med = apply(M, 2, median),
         lo  = apply(M, 2, quantile, 0.025),
         hi  = apply(M, 2, quantile, 0.975))
}
band.df <- purrr::map_dfr(mics, fit_band) |> dplyr::mutate(mic = factor(mic, levels = lev))
dat     <- dplyr::mutate(df.t, mic = factor(as.character(mic), levels = lev))

p1 <- ggplot(band.df, aes(temp, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "red", alpha = 0.15) +
  geom_line(colour = "red3", linewidth = 0.8) +
  geom_point(data = dat, aes(temp, peak.RFU), inherit.aes = FALSE, alpha = 0.4, size = 0.6) +
  facet_wrap(~ mic, scales = "free_y") +
  labs(x = "Temperature (C)", y = "peak RFU",
       title = "Peak-abundance TPCs: pointwise posterior fitted curves (median + 95%)") +
  theme_bw()

ggsave(file.path(fig.dir, "11_fig_TPC_peak_abd.png"),
       p1, width = 9, height = 6, dpi = 300)

# Nitrogen Monod curves ---------------------------------------------------

df.n <- mu %>%
  filter(salt == 0, temp == 30, block != 2, rep < 4, !is.na(peak.RFU)) %>%
  mutate(mic = factor(mic), block = factor(block), y = peak.RFU / 1000)   # k-RFU

f_n <- bf(
  y ~ mumax * nit / (ks + nit),
  mumax ~ 1 + block,
  ks    ~ 1 + block,            # block carried real signal on ks in 05 — keep it
  sigma ~ poly(log1p(nit), 2),
  nl = TRUE)

p_n <- c(
  prior(normal(6.25, 2.5), nlpar = "mumax", coef = "Intercept"),  # nls mumax (k-RFU), wide (non-sat)
  prior(normal(0,    1.0), nlpar = "mumax"),                      # block contrasts (k-RFU)
  prior(normal(897,  400), nlpar = "ks",    coef = "Intercept"),  # nls ks ~ gradient max, very wide
  prior(normal(0,    200), nlpar = "ks"),                         # block contrasts (uM)
  prior(normal(-0.5, 1  ), class = "Intercept", dpar = "sigma"),  # log resid SD ~ log(0.6)
  prior(normal(0,    1  ), class = "b",         dpar = "sigma"))

fit_one_n <- function(m) {
  d  <- droplevels(subset(df.n, mic == m))
  nb <- nlevels(d$block)
  init_n <- function() list(
    b_mumax = as.array(c(6.25, rep(0, nb - 1))),
    b_ks    = as.array(c(897,  rep(0, nb - 1))))
  brm(f_n, data = d, family = gaussian(), prior = p_n, init = init_n,
      chains = 4, iter = 3000, warmup = 1000,
      control = list(adapt_delta = 0.99, max_treedepth = 12), seed = 1,
      file = file.path(ind.dir, paste0("nit_peak_abd_mic_", m)),
      file_refit = "on_change", backend = "cmdstanr")
}

mics.n <- as.character(sort(unique(df.n$mic)))
fits.n <- setNames(lapply(mics.n, fit_one_n), mics.n)

## Diagnostics and ploting -------------------------------------------------------------

ind.diag.n <- purrr::map_dfr(mics.n, function(m)
  diag_one(fits.n[[m]]) |> dplyr::mutate(mic = m, .before = 1))
print(dplyr::mutate(ind.diag.n, dplyr::across(where(is.numeric), ~round(.x, 3))), n = Inf)

nit.grid <- seq(0, max(df.n$nit), length.out = 1000)

draws_n <- function(D) data.frame(mumax = bp(D, "mumax"), ks = bp(D, "ks"))

band_n <- function(m) {
  P <- draws_n(as.matrix(fits.n[[m]]))
  M <- vapply(nit.grid, function(N) 1000 * P$mumax * N / (P$ks + N),   # Monod inline, raw RFU
              numeric(nrow(P)))                                         # draws x grid
  tibble(mic = m, nit = nit.grid, med = apply(M, 2, median),
         lo = apply(M, 2, quantile, 0.025), hi = apply(M, 2, quantile, 0.975))
}

band.n <- purrr::map_dfr(mics.n, band_n) |> dplyr::mutate(mic = factor(mic, levels = lev))
dat.n  <- dplyr::mutate(df.n, mic = factor(as.character(mic), levels = lev))

p2 <- ggplot(band.n, aes(nit, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "red", alpha = 0.15) +
  geom_line(colour = "red3", linewidth = 0.8) +
  geom_point(data = dat.n, aes(nit, peak.RFU), inherit.aes = FALSE, alpha = 0.4, size = 0.6) +
  facet_wrap(~ mic, scales = "free_y") +
  labs(x = "Nitrogen (uM)", y = "peak RFU",
       title = "Peak-abundance Monod fits (pointwise posterior, median + 95%)") +
  theme_bw()

p2

ggsave(file.path(fig.dir, "12_fig_nit_peak_abd.png"),
       p2, width = 9, height = 6, dpi = 300)

# Salt tolerance ----------------------------------------------------------

# Running these as an inverse logistic, given the existence of a plateau-like signature at low salt levels.

df.s <- mu %>%
  filter(nit == 1000, temp == 30, block != 2, rep < 4, !is.na(peak.RFU)) %>%
  mutate(mic = factor(mic), block = droplevels(factor(block)), y = peak.RFU / 1000)

f_s <- bf(
  y ~ L + (U - L) / (1 + exp(b * (salt - m))),
  U ~ 1 + block,          # carbon blocks shift the low-salt abundance level
  L ~ 1,
  b ~ 1 + block,
  m ~ 1 + block,
  sigma ~ poly(salt, 2),  # spread shrinks toward the floor (as in 06)
  nl = TRUE)

p_s <- c(
  prior(normal(4,    1.5), nlpar = "U", coef = "Intercept"),  # low-salt plateau (k-RFU)
  prior(normal(0,    1.0), nlpar = "U"),                      # U block contrasts
  
  prior(normal(0.02, 0.1), nlpar = "L", lb = 0),              # background floor >= 0 (intercept only)
  
  prior(normal(1,    0.7), nlpar = "b", coef = "Intercept"),  # crash steepness (held > 0 by prior)
  prior(normal(0,    0.3), nlpar = "b"),                      # b block contrasts (centered 0)
  
  prior(normal(5,    3  ), nlpar = "m", coef = "Intercept"),  # crash midpoint (g/L)
  prior(normal(0,    1.0), nlpar = "m"),                      # m block contrasts (centered 0)
  
  prior(normal(-0.7, 1  ), class = "Intercept", dpar = "sigma"),
  prior(normal(0,    1  ), class = "b",         dpar = "sigma"))

init_s <- function() list(
  b_U = as.array(c(4,    rep(0, nb - 1))),
  b_L = as.array(0.02),
  b_b = as.array(c(1,    rep(0, nb - 1))),
  b_m = as.array(c(5,    rep(0, nb - 1))))

fit_one_s <- function(m) {
  d  <- droplevels(subset(df.s, mic == m))
  nb <- nlevels(d$block)
  init_s <- function() list(
    b_U = as.array(c(4,    rep(0, nb - 1))),
    b_L = as.array(0.02),
    b_b = as.array(c(1,    rep(0, nb - 1))),   
    b_m = as.array(c(5,    rep(0, nb - 1))))   
  brm(f_s, data = d, family = gaussian(), prior = p_s, init = init_s,
      chains = 4, iter = 3000, warmup = 1000,
      control = list(adapt_delta = 0.99, max_treedepth = 12), seed = 1,
      file = file.path(ind.dir, paste0("salt_peak_abd_mic_", m)),
      file_refit = "on_change", backend = "cmdstanr")
}

mics.s <- as.character(sort(unique(df.s$mic)))
fits.s <- setNames(lapply(mics.s, fit_one_s), mics.s)

# Diagnostics and plotting ------------------------------------------------

ind.diag.s <- purrr::map_dfr(mics.s, function(m)
  diag_one(fits.s[[m]]) |> dplyr::mutate(mic = m, .before = 1))
print(dplyr::mutate(ind.diag.s, dplyr::across(where(is.numeric), ~round(.x, 3))), n = Inf)

salt.grid <- seq(0, max(df.s$salt), length.out = 1000)

draws_s <- function(D) data.frame(U = bp(D, "U"), L = bp(D, "L"), b = bp(D, "b"), m = bp(D, "m"))
band_s  <- function(mic_i) {
  P <- draws_s(as.matrix(fits.s[[mic_i]]))
  M <- vapply(salt.grid,
              function(S) 1000 * (P$L + (P$U - P$L) / (1 + exp(P$b * (S - P$m)))),  # logistic inline, raw RFU
              numeric(nrow(P)))
  tibble(mic = mic_i, salt = salt.grid, med = apply(M, 2, median),
         lo = apply(M, 2, quantile, 0.025), hi = apply(M, 2, quantile, 0.975))
}

band.s <- purrr::map_dfr(mics.s, band_s) |> dplyr::mutate(mic = factor(mic, levels = lev))
dat.s  <- dplyr::mutate(df.s, mic = factor(as.character(mic), levels = lev))

p3 <- ggplot(band.s, aes(salt, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "red", alpha = 0.15) +
  geom_line(colour = "red3", linewidth = 0.8) +
  geom_point(data = dat.s, aes(salt, peak.RFU), inherit.aes = FALSE, alpha = 0.4, size = 0.6) +
  facet_wrap(~ mic, scales = "free_y") +
  labs(x = "Salt (g/L)", y = "peak RFU",
       title = "Peak-abundance salt fits (declining logistic, pointwise posterior)") +
  theme_bw()

p3

ggsave(file.path(fig.dir, "13_fig_salt_peak_abd.png"),
       p3, width = 9, height = 6, dpi = 300)

# Microbial effects across gradients --------------------------------------

lactin2  <- function(temp, a, b, tmax, dt) 1000 * (exp(a*temp) - exp(a*tmax - (tmax-temp)/dt) + b)
monod    <- function(nit, mumax, ks)       1000 * (mumax * nit / (ks + nit))
invlogis <- function(salt, U, L, b, m)     1000 * (L + (U - L) / (1 + exp(b * (salt - m))))

bavg <- function(D, p) {
  ic <- D[, paste0("b_", p, "_Intercept")]
  bl <- grep(paste0("^b_", p, "_block"), colnames(D), value = TRUE)
  if (length(bl)) ic + rowSums(D[, bl, drop = FALSE]) / (length(bl) + 1) else ic
}

curve_draws <- function(fit, grid, pars, fun) {
  D <- as.matrix(fit); P <- lapply(pars, function(p) bavg(D, p))
  sapply(grid, function(x) do.call(fun, c(list(x), P)))
}

interaction_draws <- function(fit_trt, fit_ref, grid, pars, fun) {
  A <- curve_draws(fit_trt, grid, pars, fun); R <- curve_draws(fit_ref, grid, pars, fun)
  n <- min(nrow(A), nrow(R)); A[seq_len(n), ] - R[seq_len(n), ]
}

summarise_int <- function(M, grid) {
  h <- apply(M, 2, function(x) { d <- hdi(x, ci = 0.95); c(d$CI_low, d$CI_high) })
  tibble(x = grid, med = apply(M, 2, median), lo = h[1, ], hi = h[2, ],
         sig = h[1, ] > 0 | h[2, ] < 0)
}                                                   

# fits_T <- setNames(lapply(c("none","all",1:15), \(m) readRDS(file.path("Models/individual",     paste0("TPC_peak_abd_mic_",  m,".rds")))), c("none","all",1:15))
# fits_N <- ... "Models/individual"  / "nit_peak_abd_mic_"
# fits_S <- ... "Models/individual" / "salt_peak_abd_mic_"
fits_T <- fits.ind; fits_N <- fits.n; fits_S <- fits.s

trts     <- c("all", as.character(1:15))
tpc.pars <- c("a","b","tmax","dt"); nit.pars <- c("mumax","ks"); salt.pars <- c("U","L","b","m")

temp.grid <- with(fits_T[["none"]]$data, seq(min(temp), max(temp), length.out = 1000))
nit.grid  <- with(fits_N[["none"]]$data, seq(min(nit),  max(nit),  length.out = 1000))
salt.grid <- with(fits_S[["none"]]$data, seq(min(salt), max(salt), length.out = 1000))

int_of <- function(fits, grid, pars, fun)
  purrr::map_dfr(trts, function(m)
    summarise_int(interaction_draws(fits[[m]], fits[["none"]], grid, pars, fun), grid) |>
      mutate(mic = m))

regions <- function(int, trts, dig) {                     # contiguous sig runs
  step <- diff(sort(unique(int$x)))[1]
  int |> filter(sig) |>
    mutate(dir = ifelse(med > 0, "facilitation", "antagonism")) |>
    arrange(mic, dir, x) |> group_by(mic, dir) |>
    mutate(run = cumsum(c(TRUE, diff(x) > 1.5 * step))) |>
    group_by(mic, dir, run) |>
    summarise(from = round(min(x), dig), to = round(max(x), dig),
              peak = round(med[which.max(abs(med))], 1), .groups = "drop") |>
    arrange(match(mic, trts), from) |> select(mic, dir, from, to, peak)
}

mic_cols <- c("1"="chocolate3","2"="goldenrod2","3"="skyblue","4"="olivedrab4",
              "5"="firebrick3","6"="seagreen4","7"="steelblue","8"="tan3",
              "9"="darkorchid3","10"="hotpink3","11"="brown4","12"="turquoise4",
              "13"="plum3","14"="sienna3","15"="grey25","all"="mediumblue","Others"="grey45")
sp <- c("1"="E. meliloti","2"="E. meliloti","3"="S. cerevisiae","4"="P. psychrotolerans",
        "5"="N. halotolerans","6"="P. macmurdoensis","7"="P. sulfinovorans","8"="B. aerius",
        "9"="C. braakii","10"="R. rosettiformans","11"="P. protegens","12"="S. pituitosa",
        "13"="R. qingshengii","14"="R. phycosphaerae","15"="P. putida")
lab_it <- function(b) parse(text = paste0("italic('", sp[b], "')~'(", b, ")'"))

## Temperature -------------------------------------------------------------

int_t         <- int_of(fits_T, temp.grid, tpc.pars, lactin2)
si_all        <- filter(int_t, mic == "all")
sig_regions_t <- regions(int_t, trts, 1)

cat("\n== temperature: significant interaction regions (C) ==\n")
print(as.data.frame(sig_regions_t), row.names = FALSE)

sig_any_t <- setdiff(unique(sig_regions_t$mic), "all")
cat("microbes with any significant region:", paste(sig_any_t, collapse = ", "), "\n") #all, 1, 11, 13

Fn <- curve_draws(fits_T[["none"]], temp.grid, tpc.pars, lactin2)
topt <- temp.grid[which.max(apply(Fn, 2, median))]

## Nitrogen ----------------------------------------------------------------

int_n         <- int_of(fits_N, nit.grid, nit.pars, monod)
sig_regions_n <- regions(int_n, trts, 1)

cat("\n== nitrogen: significant interaction regions (uM) ==\n")
print(as.data.frame(sig_regions_n), row.names = FALSE)

sig_any_n <- setdiff(unique(sig_regions_n$mic), "all")
cat("microbes with any significant region:", paste(sig_any_n, collapse = ", "), "\n") # all, 1, 2, 4, 7, 11, 13

## Salt --------------------------------------------------------------------

int_s         <- int_of(fits_S, salt.grid, salt.pars, invlogis)
sig_regions_s <- regions(int_s, trts, 2)

cat("\n== salt: significant interaction regions (g/L) ==\n")
print(as.data.frame(sig_regions_s), row.names = FALSE)

sig_any_s <- setdiff(unique(sig_regions_s$mic), "all")
cat("microbes with any significant region:", paste(sig_any_s, collapse = ", "), "\n") # all, 1, 2, 8, 11, 12, 13, 15

sig_union <- as.character(sort(as.integer(Reduce(union, list(sig_any_t, sig_any_n, sig_any_s)))))
brks     <- c("all", sig_union, "Others")
labs_vec <- c(expression("Microbial community"), do.call(c, lapply(sig_union, lab_it)), expression("Others"))

## Build the panels --------------------------------------------------------

ctrl_panel <- function(fits, grid, pars, fun, xcol, xlab, ttl, vx) {
  Fs <- summarise_int(curve_draws(fits[["none"]], grid, pars, fun), grid)
  raw <- fits[["none"]]$data; yv <- as.character(fits[["none"]]$formula$formula[[2]])
  ggplot() +
    geom_vline(xintercept = vx, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_ribbon(data = Fs, aes(x, ymin = lo, ymax = hi), fill = "grey40", alpha = 0.2) +
    geom_point(data = raw, aes(.data[[xcol]], .data[[yv]] * 1000), colour = "grey35", size = 1, alpha = 0.55) +
    geom_line(data = Fs, aes(x, med), colour = "black", linewidth = 1) +
    labs(x = xlab, y = "Peak abundance (RFU)", title = ttl) +
    theme_classic() +
    theme(axis.title = element_text(size = 10), axis.text = element_text(size = 10),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.03), plot.margin = margin(3,3,3,3))
}

eff_panel <- function(int, sig_reg, sig_any, xlab, ttl, vx) {
  si_all <- filter(int, mic == "all")
  bg <- filter(int, !mic %in% c(sig_any, "all")); hl <- filter(int, mic %in% sig_any)
  yr <- range(c(si_all$lo, si_all$hi, int$med)); spn <- diff(yr)
  place <- function(d, top) { d <- droplevels(mutate(d, mic = factor(mic, levels = c("all", sig_any))))
  d$row <- as.integer(d$mic); d$y <- if (top) yr[2] + 0.03*spn*d$row else yr[1] - 0.03*spn*d$row; d }
  fac <- place(filter(sig_reg, dir == "facilitation"), TRUE)
  ant <- place(filter(sig_reg, dir == "antagonism"),  FALSE)
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = vx, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
    geom_ribbon(data = si_all, aes(x, ymin = lo, ymax = hi), fill = "mediumblue", alpha = 0.2) +
    geom_line(data = bg, aes(x, med, group = mic, colour = "Others"), linewidth = 0.35) +
    geom_line(data = hl, aes(x, med, group = mic, colour = mic), linewidth = 0.8) +
    geom_line(data = si_all, aes(x, med, colour = "all"), linewidth = 1) +
    geom_segment(data = filter(fac, mic != "all"), aes(from, y, xend = to, yend = y, colour = mic),
                 linewidth = 1.2, lineend = "butt", show.legend = FALSE) +
    geom_segment(data = filter(ant, mic != "all"), aes(from, y, xend = to, yend = y, colour = mic),
                 linewidth = 1.2, lineend = "butt", show.legend = FALSE) +
    geom_segment(data = filter(fac, mic == "all"), aes(from, y, xend = to, yend = y),
                 colour = "mediumblue", linewidth = 1.2, lineend = "butt") +
    geom_segment(data = filter(ant, mic == "all"), aes(from, y, xend = to, yend = y),
                 colour = "mediumblue", linewidth = 1.2, lineend = "butt") +
    scale_colour_manual(values = mic_cols, breaks = brks, labels = labs_vec,
                        limits = brks, drop = FALSE, name = "Microbial treatment") +
    labs(x = xlab, y = "Effect on peak abundance (RFU)", title = ttl) +
    theme_classic() +
    theme(axis.title = element_text(size = 10), axis.text = element_text(size = 10),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.03), plot.margin = margin(3,3,3,3))
}

pA <- ctrl_panel(fits_T, temp.grid, tpc.pars, lactin2, "temp", "Temperature (\u00b0C)", "A) Temperature", topt)
pB <- ctrl_panel(fits_N, nit.grid,  nit.pars,  monod,   "nit",  "Nitrogen (\u00b5M)",   "B) Nitrogen", 1000)
pC <- ctrl_panel(fits_S, salt.grid, salt.pars, invlogis,"salt", expression("Salt (g L"^-1*")"), "C) Salt", 0)
pD <- eff_panel(int_t, sig_regions_t, sig_any_t, "Temperature (\u00b0C)", "D) Temperature", topt)
pE <- eff_panel(int_n, sig_regions_n, sig_any_n, "Nitrogen (\u00b5M)",   "E) Nitrogen", 1000)
pF <- eff_panel(int_s, sig_regions_s, sig_any_s, expression("Salt (g L"^-1*")"), "F) Salt", 0)

pF <- pF + geom_line(
  data = tibble(mic = rep(sig_union, each = 2), x = Inf, y = Inf),
  aes(x, y, colour = mic, group = mic),
  linewidth = 0.8, na.rm = TRUE, show.legend = TRUE)

design <- "ABC#\nDEFG"
fig <- pA + pB + pC + (pD + guides(colour = "none")) + (pE + guides(colour = "none")) + pF + guide_area() +
  plot_layout(design = design, guides = "collect", widths = c(1,1,1,0.5)) &
  theme(legend.key.size = unit(0.7,"lines"), legend.key.width = unit(1.1,"lines"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8),
        legend.background = element_blank(), legend.key = element_blank())
fig

ggsave(file.path(fig.dir, "14_fig_peak_abd_gradient_effects.png"),
       fig, width = 12, height = 7, dpi = 600)

## Summary tables ----------------------------------------------------------

sig_tbl <- function(int, gname, dig) {
  step <- diff(sort(unique(int$x)))[1]
  int |> filter(sig) |>
    mutate(dir = ifelse(med > 0, "raises", "lowers")) |>
    arrange(mic, dir, x) |> group_by(mic, dir) |>
    mutate(run = cumsum(c(TRUE, diff(x) > 1.5 * step))) |>
    group_by(mic, dir, run) |>
    summarise(from = round(min(x), dig), to = round(max(x), dig),
              x_peak   = round(x[which.max(abs(med))], dig),
              peak_RFU = round(med[which.max(abs(med))], 0),
              frac_sig = round(n() * step / diff(range(int$x)), 2),
              .groups = "drop") |>
    transmute(gradient = gname, mic, dir, from, to, x_peak, peak_RFU, frac_sig)
}

sig_summary <- bind_rows(
  sig_tbl(int_t, "temperature", 1),
  sig_tbl(int_n, "nitrogen",    0),
  sig_tbl(int_s, "salt",        2)) |>
  mutate(treatment = ifelse(mic == "all", "community", sp[mic]),
         effect    = ifelse(dir == "raises", "raises peak abundance", "lowers peak abundance")) |>
  select(gradient, mic, treatment, effect, from, to, x_peak, peak_RFU, frac_sig) |>
  arrange(match(gradient, c("temperature","nitrogen","salt")), match(mic, trts), from)

write.csv(sig_summary, file.path(proc.dir, "15_peak_abd_significant_microbial_effects.csv"), row.names = FALSE)
cat("\n== peak abundance: significant microbial effects (95% HDI) ==\n")
print(as.data.frame(sig_summary), row.names = FALSE)

# Microbial effects on traits ---------------------------------------------

derive_T <- function(a,b,tmax,dt) {
  f <- function(x) lactin2(x,a,b,tmax,dt)
  o <- tryCatch(optimize(f, c(5,45), maximum=TRUE), error=function(e) NULL)
  if (is.null(o) || !is.finite(o$objective) || o$objective <= 0) return(c(Topt=NA,peak=NA,Tbr=NA))
  To <- o$maximum; pk <- o$objective; half <- pk/2
  lo <- tryCatch(uniroot(function(x) f(x)-half, c(0,To))$root,  error=function(e) NA_real_)
  hi <- tryCatch(uniroot(function(x) f(x)-half, c(To,50))$root, error=function(e) NA_real_)
  c(Topt=To, peak=pk, Tbr=hi-lo)
}

derive_N <- function(mumax,ks) c(mumax=mumax, ks=ks, affinity=mumax/ks)

derive_S <- function(U,L,b,m) {
  abd0 <- L + (U-L)/(1+exp(-b*m)); t <- abd0/2
  st <- if (is.finite(abd0) && U>t && t>L && b>0) m + log((U-t)/(t-L))/b else NA_real_
  c(abd_max=abd0, salt_tol=st, m=m)
}

trT <- function(D){ M<-t(mapply(derive_T, bp(D,"a"),bp(D,"b"),bp(D,"tmax"),bp(D,"dt")))
M[,"peak"] <- M[,"peak"]*1000; M }

trN <- function(D){ M<-t(mapply(derive_N, bp(D,"mumax"),bp(D,"ks")))
M[,c("mumax","affinity")] <- M[,c("mumax","affinity")]*1000; M }

trS <- function(D){ M<-t(mapply(derive_S, bp(D,"U"),bp(D,"L"),bp(D,"b"),bp(D,"m")))
M[,"abd_max"] <- M[,"abd_max"]*1000; M }

mk_tr <- function(fits, ex) setNames(lapply(names(fits), \(m) ex(as.matrix(fits[[m]]))), names(fits))
tr_ind.T <- mk_tr(fits.ind, trT)
tr_ind.N <- mk_tr(fits.n,   trN)
tr_ind.S <- mk_tr(fits.s,   trS)

saveRDS(tr_ind.T, file.path(mod.dir, "07_temp_peak_abd_tr_ind.rds"))
saveRDS(tr_ind.N, file.path(mod.dir, "08_nit_peak_abd_tr_ind.rds"))
saveRDS(tr_ind.S, file.path(mod.dir, "09_salt_peak_abd_tr_ind.rds"))

summ_traits <- function(tr) purrr::map_dfr(colnames(tr), function(v) {
  x <- tr[, v]; ok <- x[is.finite(x)]; h <- bayestestR::hdi(ok, ci = 0.95)
  tibble(trait = v, median = median(ok), lo = h$CI_low, hi = h$CI_high, pNA = mean(!is.finite(x)))
})

contrast_tbl <- function(tr_list, trait.set) {
  none <- tr_list[["none"]]
  purrr::map_dfr(setdiff(names(tr_list), "none"), function(m) {
    mm <- tr_list[[m]]
    purrr::map_dfr(trait.set, function(v) {
      d <- mm[, v] - none[, v]; ok <- d[is.finite(d)]; h <- bayestestR::hdi(ok, ci = 0.95)
      tibble(mic = m, trait = v, delta = median(ok), lo = h$CI_low, hi = h$CI_high,
             excl0 = h$CI_low > 0 | h$CI_high < 0, pNA = mean(!is.finite(d)))
    })
  })
}

abd_traits <- bind_rows(
  purrr::map_dfr(names(tr_ind.T), ~ summ_traits(tr_ind.T[[.x]]) |> mutate(mic=.x, axis="temperature")),
  purrr::map_dfr(names(tr_ind.N), ~ summ_traits(tr_ind.N[[.x]]) |> mutate(mic=.x, axis="nitrogen")),
  purrr::map_dfr(names(tr_ind.S), ~ summ_traits(tr_ind.S[[.x]]) |> mutate(mic=.x, axis="salt")))
write.csv(abd_traits, file.path(proc.dir, "16_peak_abd_traits.csv"), row.names=FALSE)

abd_eff <- bind_rows(
  contrast_tbl(tr_ind.T, c("Topt","peak","Tbr"))      |> mutate(axis="temperature"),
  contrast_tbl(tr_ind.N, c("mumax","ks","affinity"))  |> mutate(axis="nitrogen"),
  contrast_tbl(tr_ind.S, c("abd_max","salt_tol","m")) |> mutate(axis="salt"))
write.csv(abd_eff, file.path(proc.dir, "17_abd_microbe_effects.csv"), row.names=FALSE)

cat("\n== abundance: treatments clearing 0 vs control, per trait ==\n")
abd_eff |> group_by(axis, trait) |> summarise(n_sig = sum(excl0), .groups="drop") |> print(n=Inf)

## Figure ------------------------------------------------------------------

tr <- list(temperature = tr_ind.T, nitrogen = tr_ind.N, salt = tr_ind.S)

emergent_axis <- function(trA) {
  mics <- setdiff(names(trA), c("none","all")); traits <- colnames(trA[["all"]])
  n <- min(vapply(trA, nrow, integer(1)))
  none <- trA[["none"]][seq_len(n), , drop=FALSE]; allc <- trA[["all"]][seq_len(n), , drop=FALSE]
  purrr::map_dfr(traits, function(v) {
    mic <- sapply(mics, function(m) trA[[m]][seq_len(n), v])          # n x 15
    E.mean <- allc[, v] - rowMeans(mic, na.rm=TRUE)
    E.sum  <- (allc[, v] - none[, v]) - rowSums(mic - none[, v], na.rm=TRUE)
    hm <- bayestestR::hdi(E.mean[is.finite(E.mean)], ci=0.95)
    hs <- bayestestR::hdi(E.sum[is.finite(E.sum)],  ci=0.95)
    tibble(trait=v, d_mean=median(E.mean,na.rm=TRUE), lo_mean=hm$CI_low, hi_mean=hm$CI_high,
           sig_mean = hm$CI_low>0 | hm$CI_high<0,
           d_sum=median(E.sum,na.rm=TRUE), lo_sum=hs$CI_low, hi_sum=hs$CI_high,
           sig_sum = hs$CI_low>0 | hs$CI_high<0)
  })
}

emergent <- purrr::imap_dfr(tr, ~ emergent_axis(.x) |> mutate(axis=.y, .before=1))
cat("\n== emergent community effects on peak abundance (vs mean & sum) ==\n")
print(as.data.frame(emergent |> mutate(across(where(is.numeric), ~round(.x,1)))), row.names=FALSE)

A <- abd_eff |> transmute(axis, trait,
                          type = ifelse(mic=="all","community_vs_none","microbe_vs_none"),
                          contrast = ifelse(mic=="all","all", paste0("mic_",mic)),
                          median = delta, lo, hi, sig = excl0)

B <- emergent |> rename(median_mean=d_mean, median_sum=d_sum) |>
  pivot_longer(-c(axis,trait), names_to=c(".value","null"), names_pattern="(median|lo|hi|sig)_(mean|sum)") |>
  transmute(axis, trait, type=paste0("emergent_",null),
            contrast=ifelse(null=="mean","all - mean(15)","all - sum(15)"), median, lo, hi, sig)

lev <- c("microbe_vs_none","community_vs_none","emergent_mean","emergent_sum")
stat.tbl <- bind_rows(A,B) |> mutate(across(c(median,lo,hi), ~signif(.x,3)),
                                     type=factor(type,levels=lev)) |> arrange(axis,trait,type,contrast)

write.csv(stat.tbl, file.path(proc.dir,"22_abd_effects_stat_sum.csv"), row.names=FALSE)
cat("\n== significant (individual model) ==\n"); print(as.data.frame(filter(stat.tbl, sig)), row.names=FALSE)

sp <- c("1"="Ensifer meliloti Em1021","2"="Ensifer meliloti Em1022","3"="Saccharomyces cerevisiae",
        "4"="Pseudarthrobacter psychrotolerans","5"="Nesterenkonia halotolerans",
        "6"="Paenisporosarcina macmurdoensis","7"="Pseudarthrobacter sulfinovorans","8"="Bacillus aerius",
        "9"="Citrobacter braakii","10"="Rhizobium rosettiformans","11"="Pseudomonas protogens",
        "12"="Sphingomonas pituitosa","13"="Rhodococcus qingshengii","14"="Rhizorhabdus phycosphaerae",
        "15"="Pseudomonas putida")
mic_labs <- c(none="Control (none)", all="Microbial community", sp)
lev <- rev(c("none", as.character(1:15), "all"))          # none top, community bottom
cols <- c(none="red2", microbe="black", community="mediumblue")

tr.sum  <- abd_traits
sig.tab <- abd_eff |> select(axis, mic, trait, excl0)

fdat <- function(ax, tr_name) {
  m <- tr.sum  |> filter(axis==ax, trait==tr_name) |> select(mic, x=median, lo, hi)
  s <- sig.tab |> filter(axis==ax, trait==tr_name) |> select(mic, sig=excl0)
  left_join(m, s, by="mic") |>
    mutate(sig = coalesce(sig, FALSE),
           type = case_when(mic=="none"~"none", mic=="all"~"community", TRUE~"microbe"),
           alpha = ifelse(type=="none" | sig, 1, 0.3),
           mic = factor(mic, levels=lev))
}

expect <- function(ax, tr_name) {
  trA <- tr[[ax]]; mics <- setdiff(names(trA), c("none","all")); n <- min(vapply(trA, nrow, integer(1)))
  median(rowMeans(sapply(mics, function(m) trA[[m]][seq_len(n), tr_name]), na.rm=TRUE), na.rm=TRUE)
}

emg_sig <- function(ax, tr_name) isTRUE(emergent$sig_mean[emergent$axis==ax & emergent$trait==tr_name])

forest_panel <- function(ax, tr_name, xlab, title, show_y=FALSE) {
  d <- fdat(ax, tr_name); none_x <- d$x[d$mic=="none"]; allrow <- filter(d, mic=="all")
  seg <- data.frame(x=expect(ax,tr_name), xend=allrow$x, mic=allrow$mic); es <- emg_sig(ax,tr_name)
  ggplot(d, aes(x, mic, colour=type)) +
    geom_vline(xintercept=none_x, linetype="dashed", linewidth=0.4, colour="grey20") +
    geom_segment(data=seg, aes(x=x, xend=xend, y=mic, yend=mic), colour="mediumblue",
                 linewidth=if(es) 0.7 else 0.4, linetype=if(es) "solid" else "dotted", inherit.aes=FALSE) +
    geom_point(data=seg, aes(x=x, y=mic), shape=21, fill=NA, colour="mediumblue", size=2.6, stroke=0.6, inherit.aes=FALSE) +
    geom_errorbarh(aes(xmin=lo, xmax=hi, alpha=alpha), height=0, linewidth=0.7) +
    geom_point(aes(alpha=alpha), size=2) + scale_alpha_identity() +
    scale_colour_manual(values=cols) + scale_y_discrete(labels=mic_labs) +
    labs(x=xlab, y=NULL, title=title) + theme_classic() +
    theme(axis.title=element_text(size=10), axis.text.x=element_text(size=9),
          axis.text.y = if(show_y) element_text(size=8, hjust=1, face="italic") else element_blank(),
          axis.ticks.y = if(show_y) element_line() else element_blank(),
          plot.title=element_text(size=10, face="bold", hjust=0.03), legend.position="none")
}

xlab_abd <- expression("Peak abundance (RFU)")
pA <- forest_panel("temperature","peak",    xlab_abd, "A) Temperature", show_y=TRUE)
pB <- forest_panel("nitrogen",   "mumax",   xlab_abd, "B) Nitrogen")
pC <- forest_panel("salt",       "abd_max", xlab_abd, "C) Salt")
pD <- forest_panel("temperature","Tbr", expression("Thermal breadth" ~ italic("T")[italic(br)] ~ "(°C)"), "D) Temperature", show_y=TRUE)
pE <- forest_panel("nitrogen","affinity", expression("Nitrogen affinity" ~ italic("\u03bc")[italic(max)]/italic("K")[italic(S)] ~ "(RFU " * mu * "M"^-1 * ")"), "E) Nitrogen")
pF <- forest_panel("salt","salt_tol", expression("Salt tolerance" ~ italic("S")[50] ~ (g~L^-1)), "F) Salt")

fig2 <- (pA | pB | pC) / (pD | pE | pF)
ggsave(file.path(fig.dir, "15_peak_abd_fig_trait_effects.png"), fig2, width=12, height=9, dpi=300)

# Correlational plots -----------------------------------------------------

mumax_trait <- c(temperature="r.max",  nitrogen="mumax", salt="mu_max")   # growth
abd_trait   <- c(temperature="peak",   nitrogen="mumax", salt="abd_max")  # yield
eff.files   <- c(temperature="Data-processed/05_TPC_microbe_effects.csv",
                 nitrogen   ="Data-processed/07_nit_monod_microbe_effects.csv",
                 salt       ="Data-processed/09_salt_microbe_effects.csv")

growth <- purrr::imap_dfr(eff.files, ~ read.csv(.x) |>
                            dplyr::filter(method == "individual", trait == mumax_trait[[.y]]) |>
                            transmute(axis = .y, mic = as.character(mic), g_delta = delta))

abund <- abd_eff |>
  dplyr::filter((axis=="temperature" & trait=="peak") |
                  (axis=="nitrogen"    & trait=="mumax") |
                  (axis=="salt"        & trait=="abd_max")) |>
  transmute(axis, mic = as.character(mic), a_delta = delta)

sp <- c("1"="E. meliloti 1","2"="E. meliloti 2","3"="S. cerevisiae","4"="P. psychrotolerans",
        "5"="N. halotolerans","6"="P. macmurdoensis","7"="P. sulfinovorans","8"="B. aerius",
        "9"="C. braakii","10"="R. rosettiformans","11"="P. protegens","12"="S. pituitosa",
        "13"="R. qingshengii","14"="R. phycosphaerae","15"="P. putida")

ty <- left_join(growth, abund, by = c("axis","mic")) |>
  group_by(axis) |>
  mutate(gz = as.numeric(scale(g_delta)), az = as.numeric(scale(a_delta))) |>
  ungroup() |>
  mutate(kind = ifelse(mic == "all", "Community", "microbe"),
         axis = factor(axis, levels = c("temperature","nitrogen","salt")),
         lab  = ifelse(mic == "all", "Community",
                       ifelse(pmax(abs(gz), abs(az)) > 1.3, sp[mic], NA)))

ty <- ty |>
  mutate(lab = ifelse(mic == "all", "'Community'",
                      ifelse(pmax(abs(gz), abs(az)) > 1.3, paste0("italic('", sp[mic], "')"), NA)))

corr <- ty |> group_by(axis) |> group_modify(~{
  ca <- cor.test(.x$gz, .x$az)                         # all 16
  s  <- .x[.x$mic != "all", ]; cn <- cor.test(s$gz, s$az)   # drop community
  tibble(r_all = ca$estimate, p_all = ca$p.value, lo_all = ca$conf.int[1], hi_all = ca$conf.int[2],
         r_no  = cn$estimate, p_no  = cn$p.value, lo_no  = cn$conf.int[1], hi_no  = cn$conf.int[2])
}) |> ungroup() |> mutate(across(where(is.numeric), ~round(.x, 3)))

print(as.data.frame(corr), row.names = FALSE)
# corrs and P-values: T—-0.165 (0.541), N— -0.300 (0.258), S—0.144 (0.596)

facet_lab <- as_labeller(c(temperature="A) Temperature", nitrogen="B) Nitrogen", salt="C) Salt"))

p.cor <- ggplot(ty, aes(az, gz)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey55", linewidth = 0.5) +
  geom_point(aes(colour = kind, size = kind)) +
  geom_text_repel(aes(label = lab), parse = TRUE, size = 2.8,
                  min.segment.length = 0, max.overlaps = Inf, segment.colour = "grey70") +
  geom_text(data = corr, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.5, size = 3) +
  facet_wrap(~ axis, nrow = 1, labeller = facet_lab) +
  scale_colour_manual(values = c(microbe = "grey30", Community = "mediumblue")) +
  scale_size_manual(values = c(microbe = 2, Community = 3.4)) +
  labs(x = "Effect on peak abundance (standardized)",
       y = "Effect on growth rate (standardized)") +
  theme_bw() + theme(legend.position = "none",
                     strip.text = element_text(face = "bold", hjust = 0.02))

p.cor

ggsave(file.path(fig.dir, "16_grt_peak_abd_t-offs.png"), p.cor, width=12, height=7, dpi=300)

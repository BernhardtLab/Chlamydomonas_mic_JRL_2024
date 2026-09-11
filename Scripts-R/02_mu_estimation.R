# =============================================================================
# 02_mu_estimation.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
#
# Estimate the per-well maximum exponential growth rate (mu) of the alga from
# its RFU trajectories, for every Chlamy-containing well. mu is the rate trait
# that feeds the downstream niche-trait curves (TPC / Monod / salt tolerance).
#
# ---- What this script produces --------------------------------------------
#   Data-processed/02_chlamy_mus.csv   (one row per Chlamy well)
#
# ---- Method ----------------------------------------------------------------
# For each well: sort by time; if the 2nd read is below the 1st (settling /
# condensation), drop the 1st and treat the 2nd as N0. An expanding window
# anchored at the start grows one read at a time; the endpoint of the
# exponential phase is the window with the steepest cumulative slope of
# log(RFU) vs time (max-slope is self-protecting against lag). mu is then the
# rate r from a 1-parameter nonlinear fit RFU ~ N0 * exp(r * days) over that
# window, with N0 fixed at the first (post-trim) read.
#
# ---- Known limitation (disclosed, not corrected) --------------------------
# With daily reads, ~39% of wells reach their peak in <=3 reads, so mu for the
# fastest growers rests on few points. The per-well `n.points` column records
# how many reads entered each fit so this can be filtered/reported downstream.
# mu can be <= 0 for wells that never grew (stress extremes) — as expected.
# =============================================================================

library(tidyverse)
library(nls.multstart)

# ---- Configuration ---------------------------------------------------------
proc.dir <- "Data-processed"
in.file  <- "01_timeseries_data.csv"    # output of 01_data_upload.R
out.file <- "02_chlamy_mus.csv"

# ---- Load ------------------------------------------------------------------
df <- read.csv(file.path(proc.dir, in.file))
df <- df %>% filter(Chlamy.y.n != "BLANK")   # blanks are not analysed here

# ---- Known data corrections (pipetting errors caught in lab notes) ---------
# Wells that received the wrong microbe are relabelled; wells that got an
# unknown inoculum (or Chlamy by accident) are dropped. See lab notebook and 
# 'notes' column in the design maps for each block
drop.wells <- c("b1.t25.p1.w19", "b1.t33.p1.w2", "b1.t33.p1.w12",
                "b1.t39.p2.w105",                       # block 1
                "b4.t30.p19.w1085", "b4.t30.p20.w1152") # block 4

b1.mic5 <- c("b1.t30.p27.w1615", "b1.t30.p28.w1641", "b1.t30.p28.w1649",
             "b1.t30.p28.w1663", "b1.t30.p28.w1671", "b1.t30.p29.w1685",
             "b1.t30.p29.w1699", "b1.t30.p29.w1716", "b1.t30.p29.w1717",
             "b1.t30.p29.w1730", "b1.t30.p29.w1735")   # actually received microbe 5

df <- df %>%
  mutate(Replicate = suppressWarnings(as.numeric(Replicate))) %>%
  filter(!unique.id %in% drop.wells) %>%
  mutate(
    Microbe = case_when(
      unique.id %in% b1.mic5          ~ "5",
      unique.id == "b3.t30.p7.w365"   ~ "8",
      unique.id == "b3.t30.p26.w1542" ~ "4",
      unique.id == "b4.t30.p30.w1764" ~ "10",
      TRUE ~ Microbe),
    Replicate = case_when(
      unique.id == "b1.t30.p29.w1735" ~ 5,                      # rep 5
      unique.id %in% b1.mic5          ~ 4,                      # the other 10 -> rep 4
      unique.id %in% c("b3.t30.p7.w365", "b3.t30.p26.w1542") ~ 4,
      unique.id == "b4.t30.p30.w1764" ~ 4,
      TRUE ~ Replicate))

# ---- Per-well mu estimator -------------------------------------------------
# di: all reads for one well, columns include RFU, days, log.RFU.
estimate_mu_well <- function(di) {
  
  di <- di[order(di$days), ]
  if (nrow(di) < 2) return(list(mu = NA_real_, peak.RFU = max(di$RFU, na.rm = TRUE),
                                time = NA_real_, n.points = nrow(di)))
  
  # Trim a settling first read: if read 2 < read 1, treat read 2 as N0.
  if (!is.na(di$RFU[2]) && di$RFU[2] < di$RFU[1]) di <- di[-1, ]
  di$N0 <- di$RFU[1]
  
  # Expanding-window cumulative slopes of log(RFU) ~ days; pick the steepest.
  tser <- unique(di$days)[-1]                      # endpoints (>= 2-point windows)
  slopes <- vapply(tser, function(z) {
    unname(coef(lm(log.RFU ~ days, data = di[di$days <= z, ]))["days"])
  }, numeric(1))
  s <- max(2L, which.max(slopes))                  # >= 3 points in the fit window
  di.th <- di[di$days <= tser[s], ]
  
  # 1-parameter nonlinear fit for the exponential rate r (N0 fixed as data).
  if (length(unique(na.omit(di.th$RFU))) == 1) {
    mu <- 0                                         # flat well -> no growth
  } else {
    fit <- tryCatch(
      nls_multstart(RFU ~ N0 * exp(r * days), data = di.th,
                    start_lower = c(r = -4.5), start_upper = c(r = 4.5),
                    iter = 500, supp_errors = "Y",
                    control = nls.control(maxiter = 200)),
      error = function(e) NULL)
    mu <- if (is.null(fit)) NA_real_ else coef(fit)[["r"]]
  }
  
  list(mu       = mu,
       peak.RFU = max(di$RFU, na.rm = TRUE),        # yield trait (was k = max OD600)
       time     = max(di.th$days),                  # window endpoint (days)
       n.points = nrow(di.th))                      # reads entering the fit
}

# ---- Run over every Chlamy well --------------------------------------------
df.chlamy <- df %>%
  filter(Chlamy.y.n == "y") %>%
  mutate(log.RFU = log(RFU + 0.001))

wells <- split(df.chlamy, df.chlamy$unique.id)

df.mu <- bind_rows(lapply(wells, function(di) {
  est <- estimate_mu_well(di)
  data.frame(
    unique.id    = di$unique.id[1],
    temp         = di$Temperature.C[1],
    nit          = di$Nitrogen.conc.uM[1],
    salt         = di$Salt.conc.g.l[1],
    block        = di$Block[1],
    carbon.added = di$carbon.added[1],
    mic          = di$Microbe[1],
    plate        = di$Plate.at.T[1],
    well         = di$Well.at.T[1],
    rep          = di$Replicate[1],
    mu           = est$mu,
    peak.RFU     = est$peak.RFU,
    time         = est$time,
    n.points     = est$n.points,
    stringsAsFactors = FALSE)
}))

# ---- Write -----------------------------------------------------------------
write.csv(df.mu, file.path(proc.dir, out.file), row.names = FALSE)
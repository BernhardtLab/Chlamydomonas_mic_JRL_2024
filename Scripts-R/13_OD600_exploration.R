# Description -------------------------------------------------------------
#
# 13_OD600_exploration.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# October 9, 2026
#
# An exploration of what variation in optical density data looks like,
# to determine whether we can fit gradient-specific models to the data
# for microbe-only wells AND wells with microbes and Chlamydomonas.
#
# What this script produces:
#
# Figures-main:
#   07_concept_fig.png            conceptual figure - microbial effects on TPC and interactions with stress gradients
#   
# Load packages -----------------------------------------------------------

library(tidyverse)
library(nls.multstart)

# Set paths ---------------------------------------------------------------

proc.dir <- "Data-processed"
fig.dir  <- "Figures-misc"
raw.file <- "01_timeseries_data.csv" 

raw <- read.csv(file.path(proc.dir, raw.file))
cat("Chlamy.y.n levels:", paste(sort(unique(raw$Chlamy.y.n)), collapse = ", "), "\n")

mo <- raw %>% filter(Chlamy.y.n == "n", Microbe != "none")
cat("microbe-only wells:", length(unique(mo$unique.id)),
    "| microbes:", paste(sort(unique(mo$Microbe)), collapse = ", "), "\n")

# Estimate exponential growth rate ---------------------------------------------

est_mu_od <- function(di) {
  di <- di[order(di$days), ]
  di <- di[is.finite(di$OD600) & is.finite(di$days), ]
  if (nrow(di) < 3) return(c(mu = NA_real_, n.points = nrow(di)))
  if (di$OD600[2] < di$OD600[1]) di <- di[-1, ]          # trim settling first read
  if (nrow(di) < 3) return(c(mu = NA_real_, n.points = nrow(di)))
  N0 <- di$OD600[1]; di$N0 <- N0
  di$logOD <- log(pmax(di$OD600, 1e-4))                  # floor keeps log finite
  tser   <- unique(di$days)[-1]                          # window endpoints (>= 2 pts)
  slopes <- vapply(tser, function(z) unname(coef(lm(logOD ~ days, data = di[di$days <= z, ]))["days"]),
                   numeric(1))
  s     <- max(2L, which.max(replace(slopes, !is.finite(slopes), -Inf)))  # steepest window (>= 3 pts)
  di.th <- di[di$days <= tser[s], ]
  if (length(unique(na.omit(di.th$OD600))) == 1) return(c(mu = 0, n.points = nrow(di.th)))
  fit <- tryCatch(
    nls_multstart(OD600 ~ N0 * exp(r * days), data = di.th,
                  start_lower = c(r = -4.5), start_upper = c(r = 4.5),
                  iter = 500, supp_errors = "Y", control = nls.control(maxiter = 200)),
    error = function(e) NULL)
  c(mu = if (is.null(fit)) NA_real_ else coef(fit)[["r"]], n.points = nrow(di.th))
}

mu.od <- purrr::map_dfr(split(mo, mo$unique.id), function(di) {
  e <- est_mu_od(di)
  tibble(unique.id = di$unique.id[1], mic = di$Microbe[1],
         temp = di$Temperature.C[1], nit = di$Nitrogen.conc.uM[1],
         salt = di$Salt.conc.g.l[1], block = di$Block[1],
         mu = e[["mu"]], n.points = e[["n.points"]])
}) %>% mutate(mic = factor(mic, levels = c(as.character(1:15), "all")))

cat(sprintf("\nmicrobe-only wells fit: %d | NA mu: %d | mu<=0: %d (%.0f%%)\n",
            nrow(mu.od), sum(is.na(mu.od$mu)),
            sum(mu.od$mu <= 0, na.rm = TRUE), 100*mean(mu.od$mu <= 0, na.rm = TRUE)))

# Plotting ----------------------------------------------------------------

mu_panel <- function(d, xvar, xlab, file) {
  ggplot(d, aes(.data[[xvar]], mu)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_point(alpha = 0.4, size = 0.6) +
    facet_wrap(~ mic, scales = "free_y") +
    labs(x = xlab, y = expression("microbial "*mu*" from OD600 (/day)"),
         title = paste("Microbe-only OD600 growth rate vs", xlab)) +
    theme_bw()
  ggsave(file.path(fig.dir, file), width = 11, height = 8, dpi = 300)
}

mu_panel(filter(mu.od, salt == 0,  nit == 1000), "temp", "Temperature (C)", "17_od_mo_mu_temp.png")
mu_panel(filter(mu.od, temp == 30, salt == 0),   "nit",  "Nitrogen (uM)",   "18_od_mo_mu_nit.png")
mu_panel(filter(mu.od, temp == 30, nit == 1000), "salt", "Salt (g/L)",      "19_od_mo_mu_salt.png")

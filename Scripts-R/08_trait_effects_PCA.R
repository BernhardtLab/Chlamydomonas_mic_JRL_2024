# Description -------------------------------------------------------------
#
# 08_trait_effects_PCA.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# September 23, 2026
#
# Multivariate analysis of mutualistic microbial signatures on C. 
# reinhardtii niche-determining traits. 
#
# This script will use Principal Components Analysis to assess
# whether some microbes boost algal performance across multiple 
# axes simultaneously, and identify emergent diversity in mutualistic
# mechanism and function.
#
# What this script produces:
#
# Figures-main:
#   02_fig2_PCA_microbial_effects.png          multivariate microbial effects on 9 traits
#
# Data-processed:
#   11_PCA_variance.csv                        PCA eigenvalues, % variance etc
#   12_PCA_loadings.csv                        PCA trait loadings (9 traits)
#
# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggrepel)

# Set paths ---------------------------------------------------------------
proc.dir     <- "Data-processed"
fig.main.dir <- "Figures-main"

eff.files <- c(
  temperature = file.path(proc.dir, "05_TPC_microbe_effects.csv"),
  nitrogen    = file.path(proc.dir, "07_nit_monod_microbe_effects.csv"),
  salt        = file.path(proc.dir, "09_salt_microbe_effects.csv"))

effects <- imap_dfr(eff.files, ~ read.csv(.x) |> mutate(axis = .y)) |>
  mutate(kind = ifelse(mic == "all", "community", "microbe"))

# Set the 9 traits (axis + data-name + display label) -------------------------
pca.spec <- tribble(
  ~axis,         ~trait,     ~label,
  "temperature", "Tmin",     "Tmin",
  "temperature", "Tmax",     "Tmax",
  "temperature", "Topt",     "Topt",
  "temperature", "breadth",  "Tbr",
  "temperature", "r.max",    "mumax(T)",
  "nitrogen",    "mumax",    "mumax(N)",
  "nitrogen",    "affinity", "affinity",
  "salt",        "mu_max",   "mumax(S)",
  "salt",        "s_lim",    "Scrit")

# PCA ---------------------------------------------------------------------

# microbe x trait matrix of effect medians (delta = microbe - control)
E <- effects |> filter(method == "individual", kind == "microbe") |>
  inner_join(pca.spec, by = c("axis", "trait")) |>
  select(mic, label, delta) |>
  pivot_wider(names_from = label, values_from = delta)
M <- as.matrix(E[, -1]); rownames(M) <- E$mic
stopifnot(sum(is.na(M)) == 0)

pc <- prcomp(M, center = TRUE, scale. = TRUE)
ve <- round(100 * summary(pc)$importance[2, 1:2], 1)

# project the community's effect vector into the same space -
allv <- effects |> filter(method == "individual", kind == "community") |>
  inner_join(pca.spec, by = c("axis", "trait")) |>
  select(label, delta) |> pivot_wider(names_from = label, values_from = delta)
allsc <- as.data.frame(predict(pc, as.data.frame(allv)[, colnames(M), drop = FALSE]))

scores <- as.data.frame(pc$x[, 1:2]); scores$mic <- rownames(pc$x)

# italic, subscripted trait labels
lab_map <- c(
  "Topt"     = "italic(T)[italic(opt)]",
  "Tmin"     = "italic(T)[italic(min)]",
  "Tmax"     = "italic(T)[italic(max)]",
  "Tbr"      = "italic(T)[italic(br)]",
  "Scrit"    = "italic(S)[italic(crit)]",
  "affinity" = "italic('affinity')",
  "mumax(T)" = "italic('\u03bc')[italic(max)]*' (T)'",
  "mumax(N)" = "italic('\u03bc')[italic(max)]*' (N)'",
  "mumax(S)" = "italic('\u03bc')[italic(max)]*' (S)'")

loads <- as.data.frame(pc$rotation[, 1:2]); loads$trait <- rownames(pc$rotation)
sf <- 0.85 * min(max(abs(scores$PC1)) / max(abs(loads$PC1)),
                 max(abs(scores$PC2)) / max(abs(loads$PC2)))
gap <- 0.12
loads <- loads |>
  mutate(expr = lab_map[trait],
         ex = PC1 * sf, ey = PC2 * sf,
         vert = abs(ey) > abs(ex),
         hj = ifelse(vert, 0.5, ifelse(ex > 0, 0, 1)),
         vj = ifelse(vert, ifelse(ey > 0, 0, 1), 0.5),
         lx = ex + ifelse(vert, 0, gap * sign(ex)),
         ly = ey + ifelse(vert, gap * sign(ey), 0))

# manual fixes
loads <- loads |> mutate(
  lx = case_when(trait == "Tmax" ~ lx - 0.10, TRUE ~ lx),
  ly = case_when(trait == "Tmax" ~ ly + 0.35, TRUE ~ ly))

# points as one frame so colour maps to a legend
pts <- bind_rows(
  transmute(scores, PC1, PC2, lab = mic,        type = "Single microbe"),
  transmute(allsc,  PC1, PC2, lab = "Community", type = "Microbial community"),
  tibble(PC1 = 0, PC2 = 0, lab = "Control",      type = "Control (none)")) |>
  mutate(type = factor(type, levels = c("Control (none)", "Single microbe",
                                        "Microbial community")))

cols <- c("Control (none)" = "red2", "Single microbe" = "black",
          "Microbial community" = "mediumblue")

# Figure ------------------------------------------------------------------

fig_pca <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_segment(data = loads, aes(0, 0, xend = ex, yend = ey),
               arrow = arrow(length = unit(0.15, "cm")), colour = "grey45", linewidth = 0.4) +
  geom_text(data = loads, aes(lx, ly, label = expr, hjust = hj, vjust = vj),
            parse = TRUE, colour = "grey15", size = 3.8) +
  geom_point(data = pts, aes(PC1, PC2, colour = type), size = 2.4) +
  geom_text_repel(data = filter(pts, type == "Single microbe"),
                  aes(PC1, PC2, label = lab), size = 3, seed = 1, max.overlaps = Inf) +
  scale_colour_manual(values = cols, name = "Microbial treatment") +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
  coord_equal() + theme_classic() +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 10),
        legend.position = c(0.85, 0.16),
        legend.background = element_blank(), legend.key = element_blank())

ggsave(file.path(fig.main.dir, "02_fig2_PCA_microbial_effects.png"),
       fig_pca, width = 7, height = 6, dpi = 300)

# Tables ------------------------------------------------------------------

pca_var <- data.frame(
  PC         = paste0("PC", seq_along(pc$sdev)),
  sd         = round(pc$sdev, 3),
  eigenvalue = round(pc$sdev^2, 3),
  prop_var   = round(100 * pc$sdev^2 / sum(pc$sdev^2), 1),
  cum_var    = round(100 * cumsum(pc$sdev^2) / sum(pc$sdev^2), 1))

write.csv(pca_var, file.path(proc.dir, "11_PCA_variance.csv"), row.names = FALSE)

pca_load <- data.frame(trait = rownames(pc$rotation),
                       round(pc$rotation, 3), check.names = FALSE)

write.csv(pca_load, file.path(proc.dir, "12_PCA_loadings.csv"), row.names = FALSE)

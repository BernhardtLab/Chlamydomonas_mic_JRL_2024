# Description -------------------------------------------------------------
#
# 11_concept_fig.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
# September 29, 2026
#
# A conceptual figure to establish key predictions and basic approach. 
#
# What this script produces:
#
# Figures-main:
#   07_concept_fig.png            conceptual figure - microbial effects on TPC and interactions with stress gradients
#   
# Load packages -----------------------------------------------------------

library(ggplot2)
library(patchwork)

# Create the panels -------------------------------------------------------

blue <- "blue3" 

# Lactin II function
lactin2 <- function(temp, a, b, tmax, dt) exp(a*temp) - exp(a*tmax - (tmax-temp)/dt) + b

bavg <- function(D, p) {
  ic <- D[, paste0("b_", p, "_Intercept")]
  bl <- grep(paste0("^b_", p, "_block"), colnames(D), value = TRUE)
  if (length(bl)) ic + rowSums(D[, bl, drop = FALSE]) / (length(bl) + 1) else ic
}

ind.t <- file.path("Models", "individual")
D  <- as.matrix(readRDS(file.path(ind.t, "TPC_mic_none.rds")))
p0 <- sapply(c("a","b","tmax","dt"), function(p) median(bavg(D, p)))
ctrl_fun <- function(T) lactin2(T, p0["a"], p0["b"], p0["tmax"], p0["dt"])

Tfine <- seq(0, p0["tmax"] + 1, length.out = 4000)
Topt  <- Tfine[which.max(ctrl_fun(Tfine))]
mumax <- max(ctrl_fun(Tfine))

s  <- 1.30                                        # breadth stretch (the one knob)
Tg <- seq(0, Topt + s*(p0["tmax"] - Topt) + 2, length.out = 1000)   # 0
ctrl <- ctrl_fun(Tg)
part <- ctrl_fun(Topt + (Tg - Topt)/s)
df   <- data.frame(T = Tg, ctrl = ctrl, part = part)

hm <- mumax/2
cross <- function(y, h) { it <- which.max(y)
c(Tg[which.min(abs(y[seq_len(it)] - h))],
  Tg[it - 1 + which.min(abs(y[it:length(y)] - h))]) }
cc <- cross(ctrl, hm); pc <- cross(part, hm)
upz <- function(y) { it <- which.max(y); Tg[it - 1 + which.min(abs(y[it:length(y)]))] }
tmc <- upz(ctrl); tmp <- upz(part)

lowz  <- function(y) { it <- which.max(y); Tg[which.min(abs(y[seq_len(it)]))] }
tminc <- lowz(ctrl)     # control Tmin (lower zero-crossing)

# Panel A: microbe expands thermal breadth and shifts Tmax, Tmin.

pA <- ggplot(df, aes(T)) +
  geom_line(aes(y = ctrl), colour = "black", linewidth = 1) +
  geom_line(aes(y = part), colour = "red3",    linewidth = 1.1) +
  annotate("segment", x = cc[1], xend = cc[2], y = hm, yend = hm,
           linetype = "dashed", colour = "black", linewidth = 0.9) +
  annotate("text", x = mean(cc), y = hm + 0.15, label = "italic(T)[italic(br)]", parse = TRUE, size = 3.3) +
  annotate("segment", x = cc[1], xend = pc[1] + 0.45, y = hm, yend = hm, colour = "red3", linewidth = 0.9,
           arrow = arrow(length = unit(2, "mm"))) +
  annotate("segment", x = cc[2], xend = pc[2] - 0.45, y = hm, yend = hm, colour = "red3", linewidth = 0.9,
           arrow = arrow(length = unit(2, "mm"))) +
  annotate("segment", x = tmc, xend = tmp - 0.15, y = 0, yend = 0, colour = "red3", linewidth = 0.9,
           arrow = arrow(length = unit(2, "mm"))) +
  annotate("text", x = mean(c(tmc, tmp)) - 6, y = 0.1, label = "italic(T)[italic(max)]",
           parse = TRUE, size = 3.3) +
  annotate("segment", x = tminc, xend = 0.15, y = 0, yend = 0, colour = "red3", linewidth = 0.9,
           arrow = arrow(length = unit(2, "mm"))) +
  annotate("text", x = mean(c(tminc, 0)) + 9, y = 0.12, label = "italic(T)[italic(min)]",
           parse = TRUE, size = 3.3) +
  labs(x = "Temperature (\u00b0C)", y = expression("Growth rate (" * day^-1 * ")"), title = "A) Expansion of thermal niche breadth") +
  coord_cartesian(ylim = c(-0.25, mumax * 1.12)) +      # clip the deep post-crash dive
  theme_classic(base_size = 10) +
  xlim(0,48) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  theme(plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(size = 10))

pA

# Panel B

mfac  <- 1.25                         # peak-growth boost factor (the one knob)
partB <- mfac * ctrl                  # vertical scale about zero: peak up, edges/Topt fixed
dfB   <- data.frame(T = Tg, ctrl = ctrl, part = partB)

pB <- ggplot(dfB, aes(T)) +
  geom_line(aes(y = ctrl), colour = "black", linewidth = 1) +
  geom_line(aes(y = part), colour = blue,    linewidth = 1.1) +
  # mu_max increase: blue up-arrow at Topt between the two peaks
  annotate("segment", x = Topt, xend = Topt, y = mumax, yend = mfac * mumax -0.02,
           colour = blue, linewidth = 0.9, arrow = arrow(length = unit(2, "mm"))) +
  annotate("text", x = Topt, y = (mumax + mfac * mumax) / 2 - 0.4,
           label = "italic('\u03bc')[max]", parse = TRUE, size = 3.3) +
  labs(x = "Temperature (\u00b0C)", y = expression("Growth rate (" * day^-1 * ")"),
       title = "B) Improved optimal performance") +
  coord_cartesian(ylim = c(-0.25, mfac * mumax * 1.12)) +
  theme_classic(base_size = 10) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  xlim(0,48) +
  theme(plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(size = 10))

pB

# Panel C

Tmin_c <- lowz(ctrl); Tmax_c <- upz(ctrl)      # control's thermal limits
Topt   <- Tg[which.max(ctrl)]
rng <- Tg >= Tmin_c & Tg <= Tmax_c             # plot only where the control TPC is valid

dfC <- data.frame(T = Tg[rng], eff = part[rng]  - ctrl[rng])   # A: breadth widener
dfD <- data.frame(T = Tg[rng], eff = partB[rng] - ctrl[rng])   # B: mu_max booster

pC <- ggplot(dfC, aes(T, eff)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = Topt, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(colour = "red3", linewidth = 1.1) +
  labs(x = "Temperature (\u00b0C)", y = expression("Microbial growth effect (" * day^-1 * ")"),
       title = "C) Facilitation increases with stress") +
  theme_classic(base_size = 10) +
  xlim(0,48) +
  theme(plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(size = 10))

pC

# Panel D

pD <- ggplot(dfD, aes(T, eff)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = Topt, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(colour = blue, linewidth = 1.1) +
  labs(x = "Temperature (\u00b0C)", y = expression("Microbial growth effect (" * day^-1 * ")"),
       title = "D) Stress limits facilitation") +
  theme_classic(base_size = 10) +
  xlim(0,48) +
  theme(plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(size = 10))

pD

fig_concept <- (pA | pB) / (pC | pD)
ggsave("Figures-main/07_concept_fig.png", fig_concept,
       width = 6, height = 5, dpi = 1000)
fig_concept

# Another version where each panel shows both red and blue ----------------

red <- "red3"; blue <- "blue3"          # red = tolerance/breadth, blue = mu_max
s <- 1.30; mfac <- 1.25

Tfine <- seq(0, p0["tmax"] + 1, length.out = 4000)
Topt  <- Tfine[which.max(ctrl_fun(Tfine))]
mumax <- max(ctrl_fun(Tfine))
Tg    <- seq(0, Topt + s*(p0["tmax"] - Topt) + 2, length.out = 700)
ctrl      <- ctrl_fun(Tg)
part_red  <- ctrl_fun(Topt + (Tg - Topt)/s)     # breadth widener (red)
part_blue <- mfac * ctrl                          # mu_max booster (blue)

hm <- mumax/2
cross <- function(y,h){it<-which.max(y);c(Tg[which.min(abs(y[seq_len(it)]-h))],Tg[it-1+which.min(abs(y[it:length(y)]-h))])}
upz  <- function(y){it<-which.max(y);Tg[it-1+which.min(abs(y[it:length(y)]))]}
lowz <- function(y){it<-which.max(y);Tg[which.min(abs(y[seq_len(it)]))]}
cc <- cross(ctrl,hm); pc <- cross(part_red,hm)
tmc <- upz(ctrl); tmp <- upz(part_red); tminc <- lowz(ctrl)
Tmin_c <- tminc; Tmax_c <- tmc

pA2 <- ggplot(data.frame(T=Tg), aes(T)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50", linewidth=0.4) +
  geom_line(aes(y=ctrl),      colour="black", linewidth=1) +
  geom_line(aes(y=part_red),  colour=red,  linewidth=1.1) +
  geom_line(aes(y=part_blue), colour=blue, linewidth=1.1) +
  # Tbr: black dashed across control FWHM + red expansion arrows
  annotate("segment", x=cc[1], xend=cc[2], y=hm, yend=hm, linetype="dashed", colour="black", linewidth=0.9) +
  annotate("text", x=mean(cc), y=hm+0.15, label="italic(T)[italic(br)]", colour = red, parse=TRUE, size=3.3) +
  annotate("segment", x=cc[1], xend=pc[1]+0.45, y=hm, yend=hm, colour=red, linewidth=0.9, arrow=arrow(length=unit(2,"mm"))) +
  annotate("segment", x=cc[2], xend=pc[2]-0.45, y=hm, yend=hm, colour=red, linewidth=0.9, arrow=arrow(length=unit(2,"mm"))) +
  # Tmin: red arrow to ~0
  annotate("segment", x=tminc, xend=0.3, y=0, yend=0, colour=red, linewidth=0.9, arrow=arrow(length=unit(2,"mm"))) +
  annotate("text", x=mean(c(tminc,0)) + 9, y=0.12, label="italic(T)[italic(min)]", colour = red, parse=TRUE, size=3.3) +
  # Tmax: red arrow control -> red
  annotate("segment", x=tmc, xend=tmp, y=0, yend=0, colour=red, linewidth=0.9, arrow=arrow(length=unit(2,"mm"))) +
  annotate("text", x=mean(c(tmc,tmp))-6, y=0.1, label="italic(T)[italic(max)]", colour = red, parse=TRUE, size=3.3) +
  # mu_max: blue up-arrow at Topt
  annotate("segment", x=Topt, xend=Topt, y=mumax, yend=mfac*mumax, colour=blue, linewidth=0.9, arrow=arrow(length=unit(2,"mm"))) +
  annotate("text", x=Topt, y=(mumax+mfac*mumax)/2 + 0.4, label="italic('\u03bc')[italic(max)]", parse=TRUE, size=3.3, colour=blue) +
  labs(x="Temperature (\u00b0C)", y=expression("Growth rate ("*day^-1*")"),
       title="A) Partners reshape the thermal niche") +
  coord_cartesian(ylim=c(-0.25, mfac*mumax*1.12)) +
  xlim(0, 48) +
  theme_classic(base_size=10) +
  theme(plot.title=element_text(size=10,face="bold"), axis.title=element_text(size=10))

pA2

rng <- Tg >= Tmin_c & Tg <= Tmax_c
dfC <- data.frame(T=Tg[rng], red=part_red[rng]-ctrl[rng], blue=part_blue[rng]-ctrl[rng])

pC2 <- ggplot(dfC, aes(T)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50", linewidth=0.4) +
  geom_vline(xintercept=Topt, linetype="dashed", colour="grey50", linewidth=0.4) +
  geom_line(aes(y=red),  colour=red,  linewidth=1.1) +
  geom_line(aes(y=blue), colour=blue, linewidth=1.1) +
  labs(x="Temperature (\u00b0C)", y=expression("Microbial growth effect ("*day^-1*")"),
       title="C) Effects across a thermal gradient") +
  theme_classic(base_size=10) +
  xlim(0, 48) +
  theme(plot.title=element_text(size=10,face="bold"), axis.title=element_text(size=10))

pC2

# Salt

expdecay <- function(salt, mmin, A, k) mmin + A * exp(-k * salt)
ind.s <- file.path("Models", "individual")           # salt fits (adjust if elsewhere)
Ds <- as.matrix(readRDS(file.path(ind.s, "salt_mic_none.rds")))
s0 <- sapply(c("mmin","A","k"), function(p) median(bavg(Ds, p)))
sctrl_fun <- function(S) expdecay(S, s0["mmin"], s0["A"], s0["k"])

Sg    <- seq(0, 15, length.out = 1000)
sctrl <- sctrl_fun(Sg)
sblue <- expdecay(Sg, s0["mmin"], s0["A"] * 1.30, s0["k"] * 1.70)   # more growth, crashes sooner
sred  <- expdecay(Sg, s0["mmin"], s0["A"] * 0.78, s0["k"] * 0.60)   # less growth, far more tolerant

scrit <- function(y) { i <- which(y <= 0)[1]; if (is.na(i)) max(Sg) else Sg[i] }
sc_c <- scrit(sctrl); sc_b <- scrit(sblue); sc_r <- scrit(sred)
mu_c <- sctrl[1]; mu_b <- sblue[1]; mu_r <- sred[1]

pB2 <- ggplot(data.frame(S = Sg), aes(S)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(aes(y = sctrl), colour = "black", linewidth = 1) +
  geom_line(aes(y = sred),  colour = red,  linewidth = 1.1) +
  geom_line(aes(y = sblue), colour = blue, linewidth = 1.1) +
  # mu_max at salt 0: blue up, red down from control
  annotate("segment", x = 0, xend = 0, y = mu_c, yend = mu_b, colour = blue,
           linewidth = 0.9, arrow = arrow(length = unit(2, "mm"))) +
  annotate("segment", x = 0, xend = 0, y = mu_c, yend = mu_r, colour = red,
           linewidth = 0.9, arrow = arrow(length = unit(2, "mm"))) +
  annotate("text", x = 1.7, y = 1.6, label = "italic('\u03bc')[italic(max)]", parse = TRUE, size = 3.3) +
  annotate("segment", x = sc_c, xend = sc_b, y = 0, yend = 0, colour = blue,
           linewidth = 0.9, arrow = arrow(length = unit(2, "mm"))) +             # decreased (left)
  annotate("segment", x = sc_c, xend = sc_r, y = 0, yend = 0, colour = red,
           linewidth = 0.9, arrow = arrow(length = unit(2, "mm"))) +             # increased (right)
  annotate("text", x = sc_c + 2, y = 0.27, label = "italic(S)[italic(crit)]", parse = TRUE, size = 3.3) +
  labs(x = expression("Salt (g L"^-1*")"), y = expression("Growth rate ("*day^-1*")"),
       title = "B) Partner-induced trade-offs") +
  coord_cartesian(ylim = c(-0.2, mu_b * 1.12)) +
  theme_classic(base_size = 10) +
  xlim(0,12) +
  theme(plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(size = 10))

pB2

rng_s <- Sg >= 0 & Sg <= sc_c                         
dfD <- data.frame(S = Sg, red = sred - sctrl, blue = sblue - sctrl)   # full range, incl. ctrl < 0


pD2 <- ggplot(dfD, aes(S)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(aes(y = red),  colour = red,  linewidth = 1.1) +
  geom_line(aes(y = blue), colour = blue, linewidth = 1.1) +
  labs(x = expression("Salt (g L"^-1*")"), y = expression("Microbial growth effect ("*day^-1*")"),
       title = "D) Effects across a salt gradient") +
  geom_vline(xintercept=0, linetype="dashed", colour="grey50", linewidth=0.4) +
  theme_classic(base_size = 10) +
  xlim(0,12) +
  theme(plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(size = 10))

pD2

fig_concept2 <- (pA2 | pB2) / (pC2 | pD2)
fig_concept2

ggsave("Figures-main/07b_concept_fig.png", fig_concept2,
       width = 8, height = 7.5, dpi = 1000)


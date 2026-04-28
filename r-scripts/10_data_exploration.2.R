# Jason R Laurich
# April 21st, 2026

# We're going to start doing some data exploration

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(cowplot)
library(emmeans)
library(lme4)
library(car)
library(ggridges)
library(vegan)
library(scales)
library(metafor)

# Load & examine the data -------------------------------------------------

df.n <- read_csv("processed-data/100a_chlamy_Monod_nit_bayes.csv") # This file has the nitrogen monod data for each replicate.
head(df.n)

df.s <- read_csv("processed-data/12a_chlamy_salt_tol_bayes.csv") # This file has the salt tolerance data for each replicate.
head(df.s)

df.t <- read_csv("processed-data/101a_chlamy_TPCs_bayes.csv") # This file has the TPC data for each replicate.
head(df.t)

df.mu <- read_csv("processed-data/02_chlamy_µs.csv")
head(df.mu)

df.summ <- df.t %>%  # combine the data!
  transmute(block = block,
            mic = mic,
            rep = rep,
            id = id,
            T.max = T.max,
            T.max.lwr = T.max.min,
            T.max.upr = T.max.max,
            T.min = T.min,
            T.min.lwr = T.min.min,
            T.min.upr = T.min.max,
            T.µ.max = r.max,
            T.µ.max.lwr = r.max.min,
            T.µ.max.upr = r.max.max,
            T.br = T.br.max - T.br.min,
            T.br.lwr = T.br.lwr,
            T.br.upr = T.br.upr,
            T.opt = T.opt) %>% 
  
  left_join(df.n %>% 
              transmute(N.µ.max = r.max.post,
                        N.µ.max.lwr = r.max.min,
                        N.µ.max.upr = r.max.max,
                        N.aff = aff.post,
                        N.aff.lwr = aff.min,
                        N.aff.upr = aff.max,
                        id = id),
            by = c("id" = "id")) %>% 
  
  left_join(df.s %>% 
              transmute(S.µ.max = r.max.post,
                        S.µ.max.lwr = r.max.min,
                        S.µ.max.upr = r.max.max,
                        S.c = c.post,
                        S.c.lwr = c.post.min,
                        S.c.upr = c.post.max,
                        id = id),
            by = c("id" = "id")) 

head(df.summ)

df.summ <- df.summ %>% 
  filter(block != 2) %>% 
  mutate(block = factor(block),
         mic =factor (mic))

df.summ$mic <- factor( # reorder the microbial treatments as a factor
  df.summ$mic,
  levels = c("none", as.character(1:15), "all")
)

# Effects of microbes on chlamy traits ------------------------------------

###### Nitrogen: µ max ######

df.summ <- df.summ %>%
  mutate(
    N.µ.max.se = (N.µ.max.upr - N.µ.max.lwr) / (2 * 1.96),
    N.µ.max.var = N.µ.max.se^2
  )

mod.n.mu.meta <- rma.mv(
  yi = N.µ.max,
  V = N.µ.max.var,
  mods = ~ mic,
  random = ~ 1 | block,
  data = df.summ
)

df.sum.n.mu <- summary(mod.n.mu.meta)
df.anova.n.mu <- as.data.frame(anova(mod.n.mu.meta))

df.anova.n.mu$term <- rownames(df.anova.n.mu)
rownames(df.anova.n.mu) <- NULL
df.anova.n.mu$trait <- "N.µ.max"

df.sum.n.mu <- data.frame(
  term = rownames(df.sum.n.mu$beta),
  estimate = as.numeric(df.sum.n.mu$beta),
  se = df.sum.n.mu$se,
  zval = df.sum.n.mu$zval,
  pval = df.sum.n.mu$pval,
  ci.lb = df.sum.n.mu$ci.lb,
  ci.ub = df.sum.n.mu$ci.ub
)

# Plotting

ctrl.mean <- as.data.frame(df.sum.n.mu)[1,2]

plot.key <- as.data.frame(df.sum.n.mu) %>%
  mutate(
    mic = sub("mic", "", term),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = pval < 0.1
  ) %>%
  select(term, mic, estimate, se, ci.lb, ci.ub, pval, effect.dir, sig)

plot.key <- plot.key %>%
  mutate(
    mic = if_else(term == "intrcpt", "none", mic),
    estimate = if_else(term == "intrcpt", 0, estimate),
    pval = if_else(term == "intrcpt", NA_real_, pval),
    effect.dir = if_else(term == "intrcpt", "control", effect.dir),
    sig = if_else(term == "intrcpt", FALSE, sig)
  )

df.plot <- df.summ %>%
  select(mic, N.µ.max) %>%
  left_join(plot.key, by = "mic") %>%
  mutate(
    mic = factor(mic),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

df.plot <- df.plot %>%
  mutate(mic = factor(
    mic,
    levels = rev(c("none", "all", as.character(1:15)))
  ))

plot.est <- plot.key %>%
  mutate(
    est.plot = ctrl.mean + estimate,
    ci.lb.plot = if_else(mic == "none", ci.lb, ctrl.mean + ci.lb),
    ci.ub.plot = if_else(mic == "none", ci.ub, ctrl.mean + ci.ub),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

p.n.mu <- ggplot(df.plot, aes(x = N.µ.max, y = mic)) +
  
  labs(x = expression("Maximum growth rate (" * italic(mu)[max] * ")"),  
       y = "Microbe", 
       title = "A — Growth rate") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  geom_segment(data = plot.est,
               aes(x = ci.lb.plot, xend = ci.ub.plot,
                   y = mic, yend = mic,
                   colour = effect.dir,
                   alpha = alpha.group),
               inherit.aes = FALSE,
               linewidth = 1) +
  
  geom_point(data = plot.est,
             aes(x = est.plot, y = mic,
                 fill = effect.dir,
                 alpha = alpha.group),
             inherit.aes = FALSE,
             shape = 21,
             colour = "black",
             size = 3) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey70"
  )) +
  
  scale_colour_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey40"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.9,
    nonsig = 0.25
  )) +
  
  xlim(1.2, 1.8) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.n.mu

###### Nitrogen: affinity ######

df.summ <- df.summ %>%
  mutate(
    N.aff.se = (N.aff.upr - N.aff.lwr) / (2 * 1.96),
    N.aff.var = N.aff.se^2
  )

mod.n.aff.meta <- rma.mv(
  yi = N.aff,
  V = N.aff.var,
  mods = ~ mic,
  random = ~ 1 | block,
  data = df.summ
)

df.sum.n.aff <- summary(mod.n.aff.meta)
df.anova.n.aff <- as.data.frame(anova(mod.n.aff.meta))

df.anova.n.aff$term <- rownames(df.anova.n.aff)
rownames(df.anova.n.aff) <- NULL
df.anova.n.aff$trait <- "N.aff"

df.sum.n.aff <- data.frame(
  term = rownames(df.sum.n.aff$beta),
  estimate = as.numeric(df.sum.n.aff$beta),
  se = df.sum.n.aff$se,
  zval = df.sum.n.aff$zval,
  pval = df.sum.n.aff$pval,
  ci.lb = df.sum.n.aff$ci.lb,
  ci.ub = df.sum.n.aff$ci.ub
)

# Plotting

ctrl.mean <- as.data.frame(df.sum.n.aff)[1,2]

plot.key <- as.data.frame(df.sum.n.aff) %>%
  mutate(
    mic = sub("mic", "", term),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = pval < 0.1
  ) %>%
  select(term, mic, estimate, se, ci.lb, ci.ub, pval, effect.dir, sig)

plot.key <- plot.key %>%
  mutate(
    mic = if_else(term == "intrcpt", "none", mic),
    estimate = if_else(term == "intrcpt", 0, estimate),
    pval = if_else(term == "intrcpt", NA_real_, pval),
    effect.dir = if_else(term == "intrcpt", "control", effect.dir),
    sig = if_else(term == "intrcpt", FALSE, sig)
  )

df.plot <- df.summ %>%
  select(mic, N.aff) %>%
  left_join(plot.key, by = "mic") %>%
  mutate(
    mic = factor(mic),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

df.plot <- df.plot %>%
  mutate(mic = factor(
    mic,
    levels = rev(c("none", "all", as.character(1:15)))
  ))

plot.est <- plot.key %>%
  mutate(
    est.plot = ctrl.mean + estimate,
    ci.lb.plot = if_else(mic == "none", ci.lb, ctrl.mean + ci.lb),
    ci.ub.plot = if_else(mic == "none", ci.ub, ctrl.mean + ci.ub),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

p.n.aff <- ggplot(df.plot, aes(x = N.aff, y = mic)) +
  
  labs(x = "Nitrogen affinity (1/Ks)",  
       y = "Microbe", 
       title = "B — Affinity") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  geom_segment(data = plot.est,
               aes(x = ci.lb.plot, xend = ci.ub.plot,
                   y = mic, yend = mic,
                   colour = effect.dir,
                   alpha = alpha.group),
               inherit.aes = FALSE,
               linewidth = 1) +
  
  geom_point(data = plot.est,
             aes(x = est.plot, y = mic,
                 fill = effect.dir,
                 alpha = alpha.group),
             inherit.aes = FALSE,
             shape = 21,
             colour = "black",
             size = 3) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey70"
  )) +
  
  scale_colour_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey40"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.9,
    nonsig = 0.25
  )) +
  
  xlim(0.01, 0.04) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.n.aff

plot_grid(p.n.mu, p.n.aff)

# Effects on salinity traits ----------------------------------------------

###### Salt: µ max ######

df.summ <- df.summ %>%
  mutate(
    S.µ.max.se = (S.µ.max.upr - S.µ.max.lwr) / (2 * 1.96),
    S.µ.max.var = S.µ.max.se^2
  )

mod.s.mu.meta <- rma.mv(
  yi = S.µ.max,
  V = S.µ.max.var,
  mods = ~ mic,
  random = ~ 1 | block,
  data = df.summ
)

df.sum.s.mu <- summary(mod.s.mu.meta)
df.anova.s.mu <- as.data.frame(anova(mod.s.mu.meta))

df.anova.s.mu$term <- rownames(df.anova.s.mu)
rownames(df.anova.s.mu) <- NULL
df.anova.s.mu$trait <- "S.µ.max"

df.sum.s.mu <- data.frame(
  term = rownames(df.sum.s.mu$beta),
  estimate = as.numeric(df.sum.s.mu$beta),
  se = df.sum.s.mu$se,
  zval = df.sum.s.mu$zval,
  pval = df.sum.s.mu$pval,
  ci.lb = df.sum.s.mu$ci.lb,
  ci.ub = df.sum.s.mu$ci.ub
)

# Plotting

ctrl.mean <- as.data.frame(df.sum.s.mu)[1,2]

plot.key <- as.data.frame(df.sum.s.mu) %>%
  mutate(
    mic = sub("mic", "", term),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = pval < 0.1
  ) %>%
  select(term, mic, estimate, se, ci.lb, ci.ub, pval, effect.dir, sig)

plot.key <- plot.key %>%
  mutate(
    mic = if_else(term == "intrcpt", "none", mic),
    estimate = if_else(term == "intrcpt", 0, estimate),
    pval = if_else(term == "intrcpt", NA_real_, pval),
    effect.dir = if_else(term == "intrcpt", "control", effect.dir),
    sig = if_else(term == "intrcpt", FALSE, sig)
  )

df.plot <- df.summ %>%
  select(mic, S.µ.max) %>%
  left_join(plot.key, by = "mic") %>%
  mutate(
    mic = factor(mic),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

df.plot <- df.plot %>%
  mutate(mic = factor(
    mic,
    levels = rev(c("none", "all", as.character(1:15)))
  ))

plot.est <- plot.key %>%
  mutate(
    est.plot = ctrl.mean + estimate,
    ci.lb.plot = if_else(mic == "none", ci.lb, ctrl.mean + ci.lb),
    ci.ub.plot = if_else(mic == "none", ci.ub, ctrl.mean + ci.ub),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

p.s.mu <- ggplot(df.plot, aes(x = S.µ.max, y = mic)) +
  
  labs(x = expression("Maximum growth rate (" * italic(mu)[max] * ")"),  
       y = "Microbe", 
       title = "A — Growth rate") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  geom_segment(data = plot.est,
               aes(x = ci.lb.plot, xend = ci.ub.plot,
                   y = mic, yend = mic,
                   colour = effect.dir,
                   alpha = alpha.group),
               inherit.aes = FALSE,
               linewidth = 1) +
  
  geom_point(data = plot.est,
             aes(x = est.plot, y = mic,
                 fill = effect.dir,
                 alpha = alpha.group),
             inherit.aes = FALSE,
             shape = 21,
             colour = "black",
             size = 3) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey70"
  )) +
  
  scale_colour_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey40"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.9,
    nonsig = 0.25
  )) +
  
  # xlim(1.2, 1.8) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.s.mu

###### Salt: µ max ######

df.summ <- df.summ %>%
  mutate(
    S.c.se = (S.c.upr - S.c.lwr) / (2 * 1.96),
    S.c.var = S.c.se^2
  )

mod.s.tol.meta <- rma.mv(
  yi = S.c,
  V = S.c.var,
  mods = ~ mic,
  random = ~ 1 | block,
  data = df.summ
)

df.sum.s.tol <- summary(mod.s.tol.meta)
df.anova.s.tol <- as.data.frame(anova(mod.s.tol.meta))

df.anova.s.tol$term <- rownames(df.anova.s.tol)
rownames(df.anova.s.tol) <- NULL
df.anova.s.tol$trait <- "S.tol"

df.sum.s.tol <- data.frame(
  term = rownames(df.sum.s.tol$beta),
  estimate = as.numeric(df.sum.s.tol$beta),
  se = df.sum.s.tol$se,
  zval = df.sum.s.tol$zval,
  pval = df.sum.s.tol$pval,
  ci.lb = df.sum.s.tol$ci.lb,
  ci.ub = df.sum.s.tol$ci.ub
)

# Plotting

ctrl.mean <- as.data.frame(df.sum.s.tol)[1,2]

plot.key <- as.data.frame(df.sum.s.tol) %>%
  mutate(
    mic = sub("mic", "", term),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = pval < 0.1
  ) %>%
  select(term, mic, estimate, se, ci.lb, ci.ub, pval, effect.dir, sig)

plot.key <- plot.key %>%
  mutate(
    mic = if_else(term == "intrcpt", "none", mic),
    estimate = if_else(term == "intrcpt", 0, estimate),
    pval = if_else(term == "intrcpt", NA_real_, pval),
    effect.dir = if_else(term == "intrcpt", "control", effect.dir),
    sig = if_else(term == "intrcpt", FALSE, sig)
  )

df.plot <- df.summ %>%
  select(mic, S.c) %>%
  left_join(plot.key, by = "mic") %>%
  mutate(
    mic = factor(mic),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

df.plot <- df.plot %>%
  mutate(mic = factor(
    mic,
    levels = rev(c("none", "all", as.character(1:15)))
  ))

plot.est <- plot.key %>%
  mutate(
    est.plot = ctrl.mean + estimate,
    ci.lb.plot = if_else(mic == "none", ci.lb, ctrl.mean + ci.lb),
    ci.ub.plot = if_else(mic == "none", ci.ub, ctrl.mean + ci.ub),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

p.s.tol <- ggplot(df.plot, aes(x = S.c, y = mic)) +
  
  labs(x = "Salt tolerance (g/L)",  
       y = "Microbe", 
       title = "B — Salt tolerance") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  geom_segment(data = plot.est,
               aes(x = ci.lb.plot, xend = ci.ub.plot,
                   y = mic, yend = mic,
                   colour = effect.dir,
                   alpha = alpha.group),
               inherit.aes = FALSE,
               linewidth = 1) +
  
  geom_point(data = plot.est,
             aes(x = est.plot, y = mic,
                 fill = effect.dir,
                 alpha = alpha.group),
             inherit.aes = FALSE,
             shape = 21,
             colour = "black",
             size = 3) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey70"
  )) +
  
  scale_colour_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey40"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.9,
    nonsig = 0.25
  )) +
  
  # xlim(1.2, 1.8) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.s.tol

# Effects on thermal traits -----------------------------------------------

###### Temperature: µ max ######

df.summ <- df.summ %>%
  mutate(
    T.µ.max.se = (T.µ.max.upr - T.µ.max.lwr) / (2 * 1.96),
    T.µ.max.var = T.µ.max.se^2
  )

mod.t.mu.meta <- rma.mv(
  yi = T.µ.max,
  V = T.µ.max.var,
  mods = ~ mic,
  random = ~ 1 | block,
  data = df.summ
)

df.sum.t.mu <- summary(mod.t.mu.meta)
df.anova.t.mu <- as.data.frame(anova(mod.t.mu.meta))

df.anova.t.mu$term <- rownames(df.anova.t.mu)
rownames(df.anova.t.mu) <- NULL
df.anova.t.mu$trait <- "T.µ.max"

df.sum.t.mu <- data.frame(
  term = rownames(df.sum.t.mu$beta),
  estimate = as.numeric(df.sum.t.mu$beta),
  se = df.sum.t.mu$se,
  zval = df.sum.t.mu$zval,
  pval = df.sum.t.mu$pval,
  ci.lb = df.sum.t.mu$ci.lb,
  ci.ub = df.sum.t.mu$ci.ub
)

# Plotting

ctrl.mean <- as.data.frame(df.sum.t.mu)[1,2]

plot.key <- as.data.frame(df.sum.t.mu) %>%
  mutate(
    mic = sub("mic", "", term),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = pval < 0.1
  ) %>%
  select(term, mic, estimate, se, ci.lb, ci.ub, pval, effect.dir, sig)

plot.key <- plot.key %>%
  mutate(
    mic = if_else(term == "intrcpt", "none", mic),
    estimate = if_else(term == "intrcpt", 0, estimate),
    pval = if_else(term == "intrcpt", NA_real_, pval),
    effect.dir = if_else(term == "intrcpt", "control", effect.dir),
    sig = if_else(term == "intrcpt", FALSE, sig)
  )

df.plot <- df.summ %>%
  select(mic, T.µ.max) %>%
  left_join(plot.key, by = "mic") %>%
  mutate(
    mic = factor(mic),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

df.plot <- df.plot %>%
  mutate(mic = factor(
    mic,
    levels = rev(c("none", "all", as.character(1:15)))
  ))

plot.est <- plot.key %>%
  mutate(
    est.plot = ctrl.mean + estimate,
    ci.lb.plot = if_else(mic == "none", ci.lb, ctrl.mean + ci.lb),
    ci.ub.plot = if_else(mic == "none", ci.ub, ctrl.mean + ci.ub),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

p.t.mu <- ggplot(df.plot, aes(x = T.µ.max, y = mic)) +
  
  labs(x = expression("Maximum growth rate (" * italic(mu)[max] * ")"),  
       y = "Microbe", 
       title = "A — Growth rate") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  geom_segment(data = plot.est,
               aes(x = ci.lb.plot, xend = ci.ub.plot,
                   y = mic, yend = mic,
                   colour = effect.dir,
                   alpha = alpha.group),
               inherit.aes = FALSE,
               linewidth = 1) +
  
  geom_point(data = plot.est,
             aes(x = est.plot, y = mic,
                 fill = effect.dir,
                 alpha = alpha.group),
             inherit.aes = FALSE,
             shape = 21,
             colour = "black",
             size = 3) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey70"
  )) +
  
  scale_colour_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey40"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.9,
    nonsig = 0.25
  )) +
  
  # xlim(1.2, 1.8) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.t.mu

###### Temperature: thermal breadth ######

df.summ <- df.summ %>%
  mutate(
    T.br.se = (T.br.upr - T.br.lwr) / (2 * 1.96),
    T.br.var = T.br.se^2
  )

mod.t.br.meta <- rma.mv(
  yi = T.br,
  V = T.br.var,
  mods = ~ mic,
  random = ~ 1 | block,
  data = df.summ
)

df.sum.t.br <- summary(mod.t.br.meta)
df.anova.t.br <- as.data.frame(anova(mod.t.br.meta))

df.anova.t.br$term <- rownames(df.anova.t.br)
rownames(df.anova.t.br) <- NULL
df.anova.t.br$trait <- "T.br"

df.sum.t.br <- data.frame(
  term = rownames(df.sum.t.br$beta),
  estimate = as.numeric(df.sum.t.br$beta),
  se = df.sum.t.br$se,
  zval = df.sum.t.br$zval,
  pval = df.sum.t.br$pval,
  ci.lb = df.sum.t.br$ci.lb,
  ci.ub = df.sum.t.br$ci.ub
)

# Plotting

ctrl.mean <- as.data.frame(df.sum.t.br)[1,2]

plot.key <- as.data.frame(df.sum.t.br) %>%
  mutate(
    mic = sub("mic", "", term),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = pval < 0.1
  ) %>%
  select(term, mic, estimate, se, ci.lb, ci.ub, pval, effect.dir, sig)

plot.key <- plot.key %>%
  mutate(
    mic = if_else(term == "intrcpt", "none", mic),
    estimate = if_else(term == "intrcpt", 0, estimate),
    pval = if_else(term == "intrcpt", NA_real_, pval),
    effect.dir = if_else(term == "intrcpt", "control", effect.dir),
    sig = if_else(term == "intrcpt", FALSE, sig)
  )

df.plot <- df.summ %>%
  select(mic, T.br) %>%
  left_join(plot.key, by = "mic") %>%
  mutate(
    mic = factor(mic),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

df.plot <- df.plot %>%
  mutate(mic = factor(
    mic,
    levels = rev(c("none", "all", as.character(1:15)))
  ))

plot.est <- plot.key %>%
  mutate(
    est.plot = ctrl.mean + estimate,
    ci.lb.plot = if_else(mic == "none", ci.lb, ctrl.mean + ci.lb),
    ci.ub.plot = if_else(mic == "none", ci.ub, ctrl.mean + ci.ub),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

p.t.br <- ggplot(df.plot, aes(x = T.br, y = mic)) +
  
  labs(x = "Thermal breadth (°C)",  
       y = "Microbe", 
       title = "B — Thermal breadth") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  geom_segment(data = plot.est,
               aes(x = ci.lb.plot, xend = ci.ub.plot,
                   y = mic, yend = mic,
                   colour = effect.dir,
                   alpha = alpha.group),
               inherit.aes = FALSE,
               linewidth = 1) +
  
  geom_point(data = plot.est,
             aes(x = est.plot, y = mic,
                 fill = effect.dir,
                 alpha = alpha.group),
             inherit.aes = FALSE,
             shape = 21,
             colour = "black",
             size = 3) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey70"
  )) +
  
  scale_colour_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey40"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.9,
    nonsig = 0.25
  )) +
  
  # xlim(1.2, 1.8) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.t.br

###### Temperature: thermal minimum ######

df.summ <- df.summ %>%
  mutate(
    T.min.se = (T.min.upr - T.min.lwr) / (2 * 1.96),
    T.min.var = T.min.se^2
  )

mod.t.min.meta <- rma.mv(
  yi = T.min,
  V = T.min.var,
  mods = ~ mic,
  random = ~ 1 | block,
  data = df.summ
)

df.sum.t.min <- summary(mod.t.min.meta)
df.anova.t.min <- as.data.frame(anova(mod.t.min.meta))

df.anova.t.min$term <- rownames(df.anova.t.min)
rownames(df.anova.t.min) <- NULL
df.anova.t.min$trait <- "T.min"

df.sum.t.min <- data.frame(
  term = rownames(df.sum.t.min$beta),
  estimate = as.numeric(df.sum.t.min$beta),
  se = df.sum.t.min$se,
  zval = df.sum.t.min$zval,
  pval = df.sum.t.min$pval,
  ci.lb = df.sum.t.min$ci.lb,
  ci.ub = df.sum.t.min$ci.ub
)

# Plotting

ctrl.mean <- as.data.frame(df.sum.t.min)[1,2]

plot.key <- as.data.frame(df.sum.t.min) %>%
  mutate(
    mic = sub("mic", "", term),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = pval < 0.1
  ) %>%
  select(term, mic, estimate, se, ci.lb, ci.ub, pval, effect.dir, sig)

plot.key <- plot.key %>%
  mutate(
    mic = if_else(term == "intrcpt", "none", mic),
    estimate = if_else(term == "intrcpt", 0, estimate),
    pval = if_else(term == "intrcpt", NA_real_, pval),
    effect.dir = if_else(term == "intrcpt", "control", effect.dir),
    sig = if_else(term == "intrcpt", FALSE, sig)
  )

df.plot <- df.summ %>%
  select(mic, T.min) %>%
  left_join(plot.key, by = "mic") %>%
  mutate(
    mic = factor(mic),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

df.plot <- df.plot %>%
  mutate(mic = factor(
    mic,
    levels = rev(c("none", "all", as.character(1:15)))
  ))

plot.est <- plot.key %>%
  mutate(
    est.plot = ctrl.mean + estimate,
    ci.lb.plot = if_else(mic == "none", ci.lb, ctrl.mean + ci.lb),
    ci.ub.plot = if_else(mic == "none", ci.ub, ctrl.mean + ci.ub),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

p.t.min <- ggplot(df.plot, aes(x = T.min, y = mic)) +
  
  labs(x = "Thermal minimum (°C)",  
       y = "Microbe", 
       title = "C — Thermal minimum") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  geom_segment(data = plot.est,
               aes(x = ci.lb.plot, xend = ci.ub.plot,
                   y = mic, yend = mic,
                   colour = effect.dir,
                   alpha = alpha.group),
               inherit.aes = FALSE,
               linewidth = 1) +
  
  geom_point(data = plot.est,
             aes(x = est.plot, y = mic,
                 fill = effect.dir,
                 alpha = alpha.group),
             inherit.aes = FALSE,
             shape = 21,
             colour = "black",
             size = 3) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey70"
  )) +
  
  scale_colour_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey40"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.9,
    nonsig = 0.25
  )) +
  
  # xlim(1.2, 1.8) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.t.min

###### Temperature: thermal maximum ######

df.summ <- df.summ %>%
  mutate(
    T.max.se = (T.max.upr - T.max.lwr) / (2 * 1.96),
    T.max.var = T.max.se^2
  )

mod.t.max.meta <- rma.mv(
  yi = T.max,
  V = T.max.var,
  mods = ~ mic,
  random = ~ 1 | block,
  data = df.summ
)

df.sum.t.max <- summary(mod.t.max.meta)
df.anova.t.max <- as.data.frame(anova(mod.t.max.meta))

df.anova.t.max$term <- rownames(df.anova.t.max)
rownames(df.anova.t.max) <- NULL
df.anova.t.max$trait <- "T.max"

df.sum.t.max <- data.frame(
  term = rownames(df.sum.t.max$beta),
  estimate = as.numeric(df.sum.t.max$beta),
  se = df.sum.t.max$se,
  zval = df.sum.t.max$zval,
  pval = df.sum.t.max$pval,
  ci.lb = df.sum.t.max$ci.lb,
  ci.ub = df.sum.t.max$ci.ub
)

# Plotting

ctrl.mean <- as.data.frame(df.sum.t.max)[1,2]

plot.key <- as.data.frame(df.sum.t.max) %>%
  mutate(
    mic = sub("mic", "", term),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = pval < 0.1
  ) %>%
  select(term, mic, estimate, se, ci.lb, ci.ub, pval, effect.dir, sig)

plot.key <- plot.key %>%
  mutate(
    mic = if_else(term == "intrcpt", "none", mic),
    estimate = if_else(term == "intrcpt", 0, estimate),
    pval = if_else(term == "intrcpt", NA_real_, pval),
    effect.dir = if_else(term == "intrcpt", "control", effect.dir),
    sig = if_else(term == "intrcpt", FALSE, sig)
  )

df.plot <- df.summ %>%
  select(mic, T.max) %>%
  left_join(plot.key, by = "mic") %>%
  mutate(
    mic = factor(mic),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

df.plot <- df.plot %>%
  mutate(mic = factor(
    mic,
    levels = rev(c("none", "all", as.character(1:15)))
  ))

plot.est <- plot.key %>%
  mutate(
    est.plot = ctrl.mean + estimate,
    ci.lb.plot = if_else(mic == "none", ci.lb, ctrl.mean + ci.lb),
    ci.ub.plot = if_else(mic == "none", ci.ub, ctrl.mean + ci.ub),
    alpha.group = if_else(sig, "sig", "nonsig")
  )

p.t.max <- ggplot(df.plot, aes(x = T.max, y = mic)) +
  
  labs(x = "Thermal maximum (°C)",  
       y = "Microbe", 
       title = "D — Thermal maximum") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  geom_segment(data = plot.est,
               aes(x = ci.lb.plot, xend = ci.ub.plot,
                   y = mic, yend = mic,
                   colour = effect.dir,
                   alpha = alpha.group),
               inherit.aes = FALSE,
               linewidth = 1) +
  
  geom_point(data = plot.est,
             aes(x = est.plot, y = mic,
                 fill = effect.dir,
                 alpha = alpha.group),
             inherit.aes = FALSE,
             shape = 21,
             colour = "black",
             size = 3) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey70"
  )) +
  
  scale_colour_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control = "grey40"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.9,
    nonsig = 0.25
  )) +
  
  # xlim(1.2, 1.8) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.t.max

# Assemble plots ----------------------------------------------------------

p.n.effects <- plot_grid(p.n.mu, p.n.aff,
                                ncol = 2,
                                align = "hv",
                                axis = "tblr")

ggsave("figures/01_nitrogen_mic_effects.jpeg", p.n.effects, width = 8, height = 5)

p.s.effects <- plot_grid(p.s.mu, p.s.tol,
                         ncol = 2,
                         align = "hv",
                         axis = "tblr")

ggsave("figures/02_salt_mic_effects.jpeg", p.s.effects, width = 8, height = 5)

p.t.effects <- plot_grid(p.t.mu, p.t.br,
                         p.t.min, p.t.max,
                         ncol = 2,
                         align = "hv",
                         axis = "tblr")

ggsave("figures/03_temp_mic_effects.jpeg", p.t.effects, width = 8, height = 8)

# Trade-off plots ---------------------------------------------------------

df.ctrl <- df.summ %>%
  filter(mic == "none")

df.all <- df.summ %>%
  filter(mic == "all")

df.mic <- df.summ %>%
  filter(!mic %in% c("none", "all"))

p.n.trade <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = N.µ.max, y = N.aff, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = N.µ.max, y = N.aff, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = N.µ.max, y = N.aff),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = N.µ.max, y = N.aff),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  geom_point(data = df.all,
             aes(x = N.µ.max, y = N.aff),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = N.µ.max, y = N.aff),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.5) +
  
  xlim(1.2, 1.75) +
  ylim(0.015, 0.065) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (N)")),
       y = "Affinity (1/Ks)")

p.n.trade

p.n.trade2 <- ggplot() +
  
  geom_point(data = df.ctrl,
             aes(x = N.µ.max, y = N.aff),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = N.µ.max, y = N.aff),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  geom_point(data = df.all,
             aes(x = N.µ.max, y = N.aff),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = N.µ.max, y = N.aff),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.5) +
  
  xlim(1.2, 1.75) +
  ylim(0.015, 0.065) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (N)")),
       y = "Affinity (1/Ks)")

p.n.trade2

p.n.trade1 <- ggplot() +
  
  geom_point(data = df.ctrl,
             aes(x = N.µ.max, y = N.aff),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = N.µ.max, y = N.aff),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  xlim(1.2, 1.75) +
  ylim(0.015, 0.065) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (N)")),
       y = "Affinity (1/Ks)")

p.n.trade1

p.s.trade <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = S.µ.max, y = S.c, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = S.µ.max, y = S.c, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = S.µ.max, y = S.c),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = S.µ.max, y = S.c),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  geom_point(data = df.all,
             aes(x = S.µ.max, y = S.c),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = S.µ.max, y = S.c),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.5) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (S)")),
       y = "Salt tolerance (g/L)")

p.s.trade

p.t.trade <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = T.µ.max, y = T.br, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = T.µ.max, y = T.br, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = T.µ.max, y = T.br),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = T.µ.max, y = T.br),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  geom_point(data = df.all,
             aes(x = T.µ.max, y = T.br),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = T.µ.max, y = T.br),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.5) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (T)")),
       y = "Thermal breadth (°C)")

p.t.trade

p.nt.trade <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = T.µ.max, y = N.µ.max, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = T.µ.max, y = N.µ.max, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = T.µ.max, y = N.µ.max),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = T.µ.max, y = N.µ.max),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  geom_point(data = df.all,
             aes(x = T.µ.max, y = N.µ.max),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = T.µ.max, y = N.µ.max),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.5) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (T)")),
       y = expression(paste("Maximum growth rate, ", mu[max], " (N)")))

p.nt.trade


p.nt2.trade <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = T.br, y = N.aff, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = T.br, y = N.aff, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = T.br, y = N.aff),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = T.br, y = N.aff),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  geom_point(data = df.all,
             aes(x = T.br, y = N.aff),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = T.br, y = N.aff),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.5) +
  
  theme_classic() +
  labs(x = "Affinity (1/Ks)",
       y = "Thermal breadth (°C)")

p.nt2.trade

p.ns.trade <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = S.µ.max, y = N.µ.max, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = S.µ.max, y = N.µ.max, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = S.µ.max, y = N.µ.max),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = S.µ.max, y = N.µ.max),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  geom_point(data = df.all,
             aes(x = S.µ.max, y = N.µ.max),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = S.µ.max, y = N.µ.max),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.5) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (S)")),
       y = expression(paste("Maximum growth rate, ", mu[max], " (N)")))

p.ns.trade

p.ns.trade2 <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = S.c, y = N.aff, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = S.c, y = N.aff, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = S.c, y = N.aff),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = S.c, y = N.aff),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  geom_point(data = df.all,
             aes(x = S.c, y = N.aff),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = S.c, y = N.aff),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.5) +
  
  theme_classic() +
  labs(x = "Salt tolerance (g/L)",
       y = "Affinity (1/Ks)")

p.ns.trade2

p.ts.trade <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = S.µ.max, y = T.µ.max, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = S.µ.max, y = T.µ.max, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = S.µ.max, y = T.µ.max),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = S.µ.max, y = T.µ.max),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  geom_point(data = df.all,
             aes(x = S.µ.max, y = T.µ.max),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = S.µ.max, y = T.µ.max),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.5) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (S)")),
       y = expression(paste("Maximum growth rate, ", mu[max], " (T)")))

p.ts.trade

p.ts.trade2 <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = S.c, y = T.br, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = S.c, y = T.br, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = S.c, y = T.br, group = mic),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = S.c, y = T.br, group = mic),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 1.5) +
  
  geom_point(data = df.all,
             aes(x = S.c, y = T.br, group = mic),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = S.c, y = T.br, group = mic),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.5) +
  
  theme_classic() +
  labs(x = "Salt tolerance (g/L)",
       y = "Thermal breadth (°C)")

p.ts.trade2

# Assemble trade-off plots ------------------------------------------------

p.n.example <- plot_grid(p.n.trade1, p.n.trade2, p.n.trade, nrow = 1)

ggsave("figures/04_nit_trade-offs_example.jpeg", p.n.example, width = 10, height = 3.5)

grid.toffs <- plot_grid(p.n.trade, p.ns.trade, p.nt.trade,
                                p.ns.trade2, p.s.trade, p.ts.trade, 
                                p.nt2.trade, p.ts.trade2, p.t.trade,
                                ncol = 3,
                                align = "hv",
                                axis = "tblr")

ggsave("figures/05_toffs.grid.jpeg", grid.toffs, width = 12, height = 12)

# SGH conceptual figure ---------------------------------------------------

pred_lact <- function(temp, a, b, delta_t, tmax){
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / delta_t)) + b
}

pred_mon <- function(res, r.max, k.s){
  r.max * res / (k.s + res)
}

pred_salt <- function(salt, a, b, c){
  a / (1 + exp(b * (salt - c)))
}

###### Temperature ######

x.t <- tibble(temp = seq(0, 50, length.out = 300))

df.t <- x.t %>%
  mutate(
    control = pred_lact(temp, .08, -0.6, 4, 40),
    facil   = pred_lact(temp, .075, -0.10, 5, 42),
    collapse= pred_lact(temp, .15, -2, 6, 39)
  ) %>%
  pivot_longer(-temp)

p1 <- ggplot(df.t, aes(temp, value, colour = name)) +
  geom_line(linewidth = 1.4) +
  labs(title = "A — Temperature",
       x = "Temperature (°C)",
       y = expression(mu)) +
  
  ylim(-1,15) +
  
  scale_colour_manual(values = c(
    control = "forestgreen",
    facil = "black",
    collapse = "firebrick3"
  )) +
  theme_classic() +
  theme(legend.position = "none")

p1

p1a <- ggplot(df.t[df.t$name == "control",], aes(temp, value, colour = name)) +
  geom_line(linewidth = 1.4) +
  labs(title = "A — Temperature",
       x = "Temperature (°C)",
       y = expression(mu)) +
  
  ylim(-1,15) +
  
  scale_colour_manual(values = c(
    control = "forestgreen",
    facil = "black",
    collapse = "firebrick3"
  )) +
  theme_classic() +
  theme(legend.position = "none")

p1a

p1b <- ggplot(df.t[df.t$name != "collapse",], aes(temp, value, colour = name)) +
  geom_line(linewidth = 1.4) +
  labs(title = "A — Temperature",
       x = "Temperature (°C)",
       y = expression(mu)) +
  
  ylim(-1,15) +
  
  scale_colour_manual(values = c(
    control = "forestgreen",
    facil = "black",
    collapse = "firebrick3"
  )) +
  theme_classic() +
  theme(legend.position = "none")

p1b

###### Nitrogen ######

x.n <- tibble(N = seq(0, 1000, length.out = 300))

df.n <- x.n %>%
  mutate(
    control = pred_mon(N, 1.0, 120),
    facil   = pred_mon(N, 0.8, 40),
    collapse= pred_mon(N, 1.2, 250)
  ) %>%
  pivot_longer(-N)

p2 <- ggplot(df.n, aes(N, value, colour = name)) +
  geom_line(linewidth = 1.4) +
  labs(title = "B — Nitrogen",
       x = expression("Nitrogen ("*mu*"M)"),
       y = expression(mu)) +
  scale_colour_manual(values = c(
    control = "forestgreen",
    facil = "black",
    collapse = "firebrick3"
  )) +
  theme_classic() +
  theme(legend.position = "none")

p2

###### Salt ######

x.s <- tibble(salt = seq(0, 12, length.out = 300))

df.s <- x.s %>%
  mutate(
    control = pred_salt(salt, 1.0, 1.2, 5),
    facil   = pred_salt(salt, 0.8, 2, 7),
    collapse= pred_salt(salt, 1.5, 0.65, 2)
  ) %>%
  pivot_longer(-salt)

p3 <- ggplot(df.s, aes(salt, value, colour = name)) +
  geom_line(linewidth = 1.4) +
  labs(title = "C — Salt",
       x = "Salt (g/L)",
       y = expression(mu)) +
  scale_colour_manual(values = c(
    control = "forestgreen",
    facil = "black",
    collapse = "firebrick3"
  )) +
  theme_classic() +
  theme(legend.position = "none")

p3

###### Predictions plot ######

df.ant <- tibble(
  dist = seq(0, 1, length.out = 300)
) %>%
  mutate(
    SGH = -0.25 + 1.15 * dist + 0.05 * sin(3 * pi * dist),   # increasing
    destabil = 0.45 * exp(-3 * dist) - 0.35 * dist - 0.05    # decreasing
  ) %>%
  pivot_longer(-dist, names_to = "hypothesis", values_to = "benefit")

p.ant <- ggplot(df.ant, aes(x = dist, y = benefit, colour = hypothesis)) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.6) +
  
  geom_line(linewidth = 1.5) +
  
  scale_colour_manual(values = c(
    SGH = "black",
    destabil = "firebrick3"
  )) +
  
  labs(
    x = "Distance from optimal conditions",
    y = "Relative benefit of microbial partners"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  
  annotate("text", x = 0.6, y = 0.75,
           label = "Stress-gradient \nhypothesis",
           colour = "black", size = 5, fontface='bold') +
  annotate("text", x = 0.57, y = -0.1,
           label = "Stress destabilization",
           colour = "firebrick3", hjust = 0, size = 5, fontface='bold')

p.ant

# Assemble the plots ------------------------------------------------------

p.temp.sgh <- plot_grid(p1a, p1b, p1, nrow = 1)

ggsave("figures/06_conceptual_fig_sgh.temp.jpeg", p.temp.sgh, width = 12, height = 4)

p.sgh <- plot_grid(p1, p2, p3, nrow = 1)

ggsave("figures/07_conceptual_fig_sgh.jpeg", p.sgh, width = 12, height = 4)


# Stress gradient hypotheses ----------------------------------------------

###### Temperature ######

df.sgh.t <- df.mu %>% 
  filter(salt == 0,
         nit == 1000)

ctrl.topt <- df.summ %>%
  filter(mic == "none") %>%
  summarise(ctrl.topt = mean(T.opt, na.rm = TRUE)) %>%
  pull(ctrl.topt)

df.ctrl <- df.sgh.t %>%
  filter(mic == "none") %>%
  group_by(temp) %>%
  summarise(mu.ctrl = mean(µ, na.rm = TRUE),
            .groups = "drop")
  
df.sgh.t <- df.sgh.t %>%
  left_join(df.ctrl, by = "temp")

df.sgh.t <- df.sgh.t %>%
  mutate(rel.mu = µ / mu.ctrl)

df.sgh.t <- df.sgh.t %>%
  mutate(dist.opt = abs(temp - ctrl.topt))

df.all <- df.sgh.t %>%
  filter(mic == "all")

df.mic <- df.sgh.t %>%
  filter(!mic %in% c("none", "all"))

p.sgh.t1 <- ggplot(df.all, aes(x = dist.opt, y = rel.mu, group = mic)) +
  
  geom_jitter(alpha = 0.5) +
  
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  theme_classic() +
  
  geom_smooth(se = FALSE, colour = "navyblue", linewidth = 1.2) +
  
  labs(
    x = "Distance from thermal optimum (°C)",
    y = "Microbial effects on algal growth rate"
  )

p.sgh.t1

p.sgh.t <- ggplot() +
  
  geom_smooth(data = df.mic,
              aes(x = dist.opt, y = rel.mu, group = mic),
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.6,
              span = 0.9) +
  
  geom_smooth(data = df.all,
              aes(x = dist.opt, y = rel.mu),
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.2,
              span = 0.9) +
  
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  labs(
    x = "Distance from thermal optimum (°C)",
    y = "Microbial effects on algal growth rate"
  ) +
  
  theme_classic()

p.sgh.t

sgh.t <- plot_grid(p.sgh.t1, p.sgh.t)

ggsave("figures/08_temp.sgh.jpeg", sgh.t, width = 8, height = 4)

###### Nitrogen ######

df.sgh.n <- df.mu %>% 
  filter(salt == 0,
         temp == 30)

df.ctrl <- df.sgh.n %>%
  filter(mic == "none") %>%
  group_by(nit) %>%
  summarise(mu.ctrl = mean(µ, na.rm = TRUE),
            .groups = "drop")

df.sgh.n <- df.sgh.n %>%
  left_join(df.ctrl, by = "nit")

df.sgh.n <- df.sgh.n %>%
  mutate(rel.mu = µ / mu.ctrl)

df.sgh.n <- df.sgh.n %>%
  mutate(dist.opt = abs(nit - 1000))

df.all <- df.sgh.n %>%
  filter(mic == "all")

df.mic <- df.sgh.n %>%
  filter(!mic %in% c("none", "all"))

p.sgh.n <- ggplot() +
  
  geom_smooth(data = df.mic,
              aes(x = dist.opt, y = rel.mu, group = mic),
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.6,
              span = 0.9) +
  
  geom_smooth(data = df.all,
              aes(x = dist.opt, y = rel.mu),
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.2,
              span = 0.9) +
  
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  labs(
    x = "Distance from optimal conditions (µM N)",
    y = "Microbial effects on algal growth rate"
  ) +
  
  theme_classic()

p.sgh.n

###### Salt ######

df.sgh.s <- df.mu %>% 
  filter(nit == 1000,
         temp == 30)

df.ctrl <- df.sgh.s %>%
  filter(mic == "none") %>%
  group_by(salt) %>%
  summarise(mu.ctrl = mean(µ, na.rm = TRUE),
            .groups = "drop")

df.sgh.s <- df.sgh.s %>%
  left_join(df.ctrl, by = "salt")

df.sgh.s <- df.sgh.s %>%
  mutate(rel.mu = µ / mu.ctrl)

df.sgh.s <- df.sgh.s %>%
  mutate(dist.opt = abs(salt - 0))

df.all <- df.sgh.s %>%
  filter(mic == "all")

df.mic <- df.sgh.s %>%
  filter(!mic %in% c("none", "all"))

p.sgh.s <- ggplot() +
  
  geom_smooth(data = df.mic,
              aes(x = dist.opt, y = rel.mu, group = mic),
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.6,
              span = 0.9) +
  
  geom_smooth(data = df.all,
              aes(x = dist.opt, y = rel.mu),
              se = FALSE,
              colour = "navyblue",
              linewidth = 1.2,
              span = 0.9) +
  
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  labs(
    x = "Distance from optimal conditions (g/L salt)",
    y = "Microbial effects on algal growth rate"
  ) +
  
  theme_classic()

p.sgh.s

# Assemble the plot -------------------------------------------------------

sgh.plots <- plot_grid(p.sgh.t, p.sgh.n, p.sgh.s, nrow = 1)

ggsave("figures/09_conceptual_fig_sgh.jpeg", sgh.plots, width = 12, height = 4)

# Miscellaneous exploration -----------------------------------------------

###### raw control growth? ######

df.c <- df.mu %>% 
  filter(mic == "none")

c1 <- ggplot(df.c %>% filter(salt == 0, nit == 1000), aes(x = temp, y = µ, colour = block)) + 
  
  geom_point() + 
  
  theme(legend.position = "none")

c2 <- ggplot(df.c %>% filter(salt == 0, temp == 30), aes(x = nit, y = µ, colour = block)) + 
  
  geom_point() + 
  
  theme(legend.position = "none")

c3 <- ggplot(df.c %>% filter(nit == 1000, temp == 30), aes(x = salt, y = µ, colour = block)) + 
  
  geom_point()

plot_grid(c1, c2, c3, nrow = 1)

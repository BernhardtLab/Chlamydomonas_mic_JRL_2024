# Jason R Laurich
# April 7th, 2026

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

# Load & examine the data -------------------------------------------------

df.n <- read_csv("processed-data/06a_chlamy_Monod_nit.csv") # This file has the nitrogen monod data for each replicate.
head(df.n)

df.s <- read_csv("processed-data/07a_chlamy_salt_tol.csv") # This file has the salt tolerance data for each replicate.
head(df.s)

df.t <- read_csv("processed-data/05a_chlamy_TPCs.csv") # This file has the TPC data for each replicate.
head(df.t)

df.summ <- df.t %>%  # combine the data!
  transmute(block = block,
            mic = mic,
            rep = rep,
            id = id,
            T.max = T.max,
            T.min = T.min,
            T.µ.max = r.max,
            T.br = T.br.max - T.br.min,
            T.br.min = T.br.min,
            T.br.max = T.br.max,
            T.opt = T.opt,
            T.a = a,
            T.b = b,
            T.tmax = tmax,
            T.d.t = d.t) %>% 
  
  left_join(df.n %>% 
              transmute(N.µ.max = r.max.post,
                        N = R.post,
                        N.comp = 1/N,
                        N.K.s = K.s.post,
                        id = id),
            by = c("id" = "id")) %>% 
  
  left_join(df.s %>% 
              transmute(S.µ.max = r.max.post,
                        S.b = b.post,
                        S.c = c.post,
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

mod.n.mu <- lmer(N.µ.max ~ mic + (1 | block), data = df.summ, REML = TRUE)

summary(mod.n.mu)
df.anova.n.mu <- as.data.frame(Anova(mod.n.mu, type = 3))

df.anova.n.mu$term <- rownames(df.anova.n.mu)
rownames(df.anova.n.mu) <- NULL
df.anova.n.mu$trait <- "N.µ.max"

emm.n.mu <- emmeans(mod.n.mu, ~ mic)
con.n.mu <- contrast(emm.n.mu, method = "trt.vs.ctrl", ref = "none")

df.mod.n.mu <- as.data.frame(con.n.mu)
df.mod.n.mu$trait <- "N.µ.max"

# Plotting

ctrl.mean <- df.summ %>% # calculate the mean of the control group
  filter(mic == "none") %>%
  summarise(ctrl = mean(N.µ.max, na.rm = TRUE)) %>%
  pull(ctrl)

plot.key <- as.data.frame(con.n.mu) %>%
  mutate(
    mic = sub(" - none", "", contrast),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = p.value < 0.05
  ) %>%
  select(mic, estimate, p.value, effect.dir, sig)

plot.key <- bind_rows(
  tibble(
    mic = "none",
    estimate = 0,
    p.value = NA_real_,
    effect.dir = "control",
    sig = FALSE
  ),
  plot.key
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

p.n.mu <- ggplot(df.plot, aes(x = N.µ.max, y = mic)) +
  
  labs(x = expression("Maximum growth rate (" * italic(mu)[max] * ")"),  
       y = "Microbe", 
       title = "A — Growth rate") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  ggridges::geom_density_ridges(
    aes(fill = effect.dir,
        alpha = alpha.group),
    scale = 0.9,
    rel_min_height = 0.01,
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control  = "grey70"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.81,
    nonsig = 0.15
  )) +
  
  xlim(0.5, 3) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.n.mu

###### Nitrogen: 1/N* ######

mod.n.comp <- lmer(N.comp ~ mic + (1 | block), data = df.summ, REML = TRUE)

summary(mod.n.comp)
df.anova.n.comp <- as.data.frame(Anova(mod.n.comp, type = 3))

df.anova.n.comp$term <- rownames(df.anova.n.comp)
rownames(df.anova.n.comp) <- NULL
df.anova.n.comp$trait <- "N.comp"

emm.n.comp <- emmeans(mod.n.comp, ~ mic)
con.n.comp <- contrast(emm.n.comp, method = "trt.vs.ctrl", ref = "none")

df.mod.n.comp <- as.data.frame(con.n.comp)
df.mod.n.comp$trait <- "N.comp"

# Plotting

ctrl.mean <- df.summ %>% # calculate the mean of the control group
  filter(mic == "none") %>%
  summarise(ctrl = mean(N.comp, na.rm = TRUE)) %>%
  pull(ctrl)

plot.key <- as.data.frame(con.n.comp) %>%
  mutate(
    mic = sub(" - none", "", contrast),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = p.value < 0.05
  ) %>%
  select(mic, estimate, p.value, effect.dir, sig)

plot.key <- bind_rows(
  tibble(
    mic = "none",
    estimate = 0,
    p.value = NA_real_,
    effect.dir = "control",
    sig = FALSE
  ),
  plot.key
)

df.plot <- df.summ %>%
  select(mic, N.comp) %>%
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

p.n.comp <- ggplot(df.plot, aes(x = N.comp, y = mic)) +
  
  labs(x = "Competitive ability (1/N*)",  
       y = "Microbe", 
       title = "B — 1/N*") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  ggridges::geom_density_ridges(
    aes(fill = effect.dir,
        alpha = alpha.group),
    scale = 0.9,
    rel_min_height = 0.01,
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control  = "grey70"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.81,
    nonsig = 0.15
  )) +
  
  #xlim(0.5, 3) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.n.comp

###### Salt: µ max ######

mod.s.mu <- lmer(S.µ.max ~ mic + (1 | block), data = df.summ, REML = TRUE)

summary(mod.s.mu)
df.anova.s.mu <- as.data.frame(Anova(mod.s.mu, type = 3))

df.anova.s.mu$term <- rownames(df.anova.s.mu)
rownames(df.anova.s.mu) <- NULL
df.anova.s.mu$trait <- "S.µ.max"

emm.s.mu <- emmeans(mod.s.mu, ~ mic)
con.s.mu <- contrast(emm.s.mu, method = "trt.vs.ctrl", ref = "none")

df.mod.s.mu <- as.data.frame(con.s.mu)
df.mod.s.mu$trait <- "S.µ.max"

# Plotting

ctrl.mean <- df.summ %>% # calculate the mean of the control group
  filter(mic == "none") %>%
  summarise(ctrl = mean(S.µ.max, na.rm = TRUE)) %>%
  pull(ctrl)

plot.key <- as.data.frame(con.s.mu) %>%
  mutate(
    mic = sub(" - none", "", contrast),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = p.value < 0.05
  ) %>%
  select(mic, estimate, p.value, effect.dir, sig)

plot.key <- bind_rows(
  tibble(
    mic = "none",
    estimate = 0,
    p.value = NA_real_,
    effect.dir = "control",
    sig = FALSE
  ),
  plot.key
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

p.s.mu <- ggplot(df.plot, aes(x = S.µ.max, y = mic)) +
  
  labs(x = expression("Maximum growth rate (" * italic(mu)[max] * ")"),  
       y = "Microbe", 
       title = "A — Growth rate") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  ggridges::geom_density_ridges(
    aes(fill = effect.dir,
        alpha = alpha.group),
    scale = 0.9,
    rel_min_height = 0.01,
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control  = "grey70"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.81,
    nonsig = 0.15
  )) +
  
  xlim(0, 5) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.s.mu

###### Salt: tolerance ######

mod.s.c <- lmer(S.c ~ mic + (1 | block), data = df.summ, REML = TRUE)

summary(mod.s.c)
df.anova.s.c <- as.data.frame(Anova(mod.s.c, type = 3))

df.anova.s.c$term <- rownames(df.anova.s.c)
rownames(df.anova.s.c) <- NULL
df.anova.s.c$trait <- "S.tolerance"

emm.s.c <- emmeans(mod.s.c, ~ mic)
con.s.c <- contrast(emm.s.c, method = "trt.vs.ctrl", ref = "none")

df.mod.s.c <- as.data.frame(con.s.c)
df.mod.s.c$trait <- "S.tolerance"

# Plotting

ctrl.mean <- df.summ %>% # calculate the mean of the control group
  filter(mic == "none") %>%
  summarise(ctrl = mean(S.c, na.rm = TRUE)) %>%
  pull(ctrl)

plot.key <- as.data.frame(con.s.c) %>%
  mutate(
    mic = sub(" - none", "", contrast),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = p.value < 0.05
  ) %>%
  select(mic, estimate, p.value, effect.dir, sig)

plot.key <- bind_rows(
  tibble(
    mic = "none",
    estimate = 0,
    p.value = NA_real_,
    effect.dir = "control",
    sig = FALSE
  ),
  plot.key
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

p.s.c <- ggplot(df.plot, aes(x = S.c, y = mic)) +
  
  labs(x = "Salt tolerance",  
       y = "Microbe", 
       title = "B — Salt tolerance") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  ggridges::geom_density_ridges(
    aes(fill = effect.dir,
        alpha = alpha.group),
    scale = 0.9,
    rel_min_height = 0.01,
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control  = "grey70"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.81,
    nonsig = 0.15
  )) +
  
  xlim(-0.5, 5) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.s.c

###### Temperature: µ max ######

mod.t.mu <- lmer(T.µ.max ~ mic + (1 | block), data = df.summ, REML = TRUE)

summary(mod.t.mu)
df.anova.t.mu <- as.data.frame(Anova(mod.t.mu, type = 3))

df.anova.t.mu$term <- rownames(df.anova.t.mu)
rownames(df.anova.t.mu) <- NULL
df.anova.t.mu$trait <- "T.µ.max"

emm.t.mu <- emmeans(mod.t.mu, ~ mic)
con.t.mu <- contrast(emm.t.mu, method = "trt.vs.ctrl", ref = "none")

df.mod.t.mu <- as.data.frame(con.t.mu)
df.mod.t.mu$trait <- "T.µ.max"

# Plotting

ctrl.mean <- df.summ %>% # calculate the mean of the control group
  filter(mic == "none") %>%
  summarise(ctrl = mean(T.µ.max, na.rm = TRUE)) %>%
  pull(ctrl)

plot.key <- as.data.frame(con.t.mu) %>%
  mutate(
    mic = sub(" - none", "", contrast),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = p.value < 0.05
  ) %>%
  select(mic, estimate, p.value, effect.dir, sig)

plot.key <- bind_rows(
  tibble(
    mic = "none",
    estimate = 0,
    p.value = NA_real_,
    effect.dir = "control",
    sig = FALSE
  ),
  plot.key
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

p.t.mu <- ggplot(df.plot, aes(x = T.µ.max, y = mic)) +
  
  labs(x = expression("Maximum growth rate (" * italic(mu)[max] * ")"),  
       y = "Microbe", 
       title = "A — Growth rate") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  ggridges::geom_density_ridges(
    aes(fill = effect.dir,
        alpha = alpha.group),
    scale = 0.9,
    rel_min_height = 0.01,
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control  = "grey70"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.81,
    nonsig = 0.15
  )) +
  
  xlim(1.25, 2.5) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.t.mu

###### Temperature: T min ######

mod.t.min <- lmer(T.min ~ mic + (1 | block), data = df.summ, REML = TRUE)

summary(mod.t.min)
df.anova.t.min <- as.data.frame(Anova(mod.t.min, type = 3))

df.anova.t.min$term <- rownames(df.anova.t.min)
rownames(df.anova.t.min) <- NULL
df.anova.t.min$trait <- "T.min"

emm.t.min <- emmeans(mod.t.min, ~ mic)
con.t.min <- contrast(emm.t.min, method = "trt.vs.ctrl", ref = "none")

df.mod.t.min <- as.data.frame(con.t.min)
df.mod.t.min$trait <- "T.µ.max"

# Plotting

ctrl.mean <- df.summ %>% # calculate the mean of the control group
  filter(mic == "none") %>%
  summarise(ctrl = mean(T.min, na.rm = TRUE)) %>%
  pull(ctrl)

plot.key <- as.data.frame(con.t.min) %>%
  mutate(
    mic = sub(" - none", "", contrast),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = p.value < 0.05
  ) %>%
  select(mic, estimate, p.value, effect.dir, sig)

plot.key <- bind_rows(
  tibble(
    mic = "none",
    estimate = 0,
    p.value = NA_real_,
    effect.dir = "control",
    sig = FALSE
  ),
  plot.key
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

p.t.min <- ggplot(df.plot, aes(x = T.min, y = mic)) +
  
  labs(x = "Minimum temperature with positive growth",  
       y = "Microbe", 
       title = "B — T min") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  ggridges::geom_density_ridges(
    aes(fill = effect.dir,
        alpha = alpha.group),
    scale = 0.9,
    rel_min_height = 0.01,
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  scale_fill_manual(values = c(
    positive = "firebrick",
    negative = "dodgerblue3",
    control  = "grey70"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.81,
    nonsig = 0.15
  )) +
  
  #xlim(1.25, 2.5) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.t.min

###### Temperature: T max ######

mod.t.max <- lmer(T.max ~ mic + (1 | block), data = df.summ, REML = TRUE)

summary(mod.t.max)
df.anova.t.max <- as.data.frame(Anova(mod.t.max, type = 3))

df.anova.t.max$term <- rownames(df.anova.t.max)
rownames(df.anova.t.max) <- NULL
df.anova.t.max$trait <- "T.max"

emm.t.max <- emmeans(mod.t.max, ~ mic)
con.t.max <- contrast(emm.t.max, method = "trt.vs.ctrl", ref = "none")

df.mod.t.max <- as.data.frame(con.t.max)
df.mod.t.max$trait <- "T.max"

# Plotting

ctrl.mean <- df.summ %>% # calculate the mean of the control group
  filter(mic == "none") %>%
  summarise(ctrl = mean(T.max, na.rm = TRUE)) %>%
  pull(ctrl)

plot.key <- as.data.frame(con.t.max) %>%
  mutate(
    mic = sub(" - none", "", contrast),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = p.value < 0.05
  ) %>%
  select(mic, estimate, p.value, effect.dir, sig)

plot.key <- bind_rows(
  tibble(
    mic = "none",
    estimate = 0,
    p.value = NA_real_,
    effect.dir = "control",
    sig = FALSE
  ),
  plot.key
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

p.t.max <- ggplot(df.plot, aes(x = T.max, y = mic)) +
  
  labs(x = " Maximum temperature for growth",  
       y = "Microbe", 
       title = "C — T max") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  ggridges::geom_density_ridges(
    aes(fill = effect.dir,
        alpha = alpha.group),
    scale = 0.9,
    rel_min_height = 0.01,
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control  = "grey70"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.81,
    nonsig = 0.15
  )) +
  
  #xlim(1.25, 2.5) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.t.max

###### Temperature: T opt ######

mod.t.opt <- lmer(T.opt ~ mic + (1 | block), data = df.summ, REML = TRUE)

summary(mod.t.opt)
df.anova.t.opt <- as.data.frame(Anova(mod.t.opt, type = 3))

df.anova.t.opt$term <- rownames(df.anova.t.opt)
rownames(df.anova.t.opt) <- NULL
df.anova.t.opt$trait <- "T.opt"

emm.t.opt <- emmeans(mod.t.opt, ~ mic)
con.t.opt <- contrast(emm.t.opt, method = "trt.vs.ctrl", ref = "none")

df.mod.t.opt <- as.data.frame(con.t.opt)
df.mod.t.opt$trait <- "T.opt"

# Plotting

ctrl.mean <- df.summ %>% # calculate the mean of the control group
  filter(mic == "none") %>%
  summarise(ctrl = mean(T.opt, na.rm = TRUE)) %>%
  pull(ctrl)

plot.key <- as.data.frame(con.t.opt) %>%
  mutate(
    mic = sub(" - none", "", contrast),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = p.value < 0.05
  ) %>%
  select(mic, estimate, p.value, effect.dir, sig)

plot.key <- bind_rows(
  tibble(
    mic = "none",
    estimate = 0,
    p.value = NA_real_,
    effect.dir = "control",
    sig = FALSE
  ),
  plot.key
)

df.plot <- df.summ %>%
  select(mic, T.opt) %>%
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

p.t.opt <- ggplot(df.plot, aes(x = T.opt, y = mic)) +
  
  labs(x = "Optimal temperature for growth",  
       y = "Microbe", 
       title = "D — T opt") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  ggridges::geom_density_ridges(
    aes(fill = effect.dir,
        alpha = alpha.group),
    scale = 0.9,
    rel_min_height = 0.01,
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control  = "grey70"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.81,
    nonsig = 0.15
  )) +
  
  #xlim(1.25, 2.5) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.t.opt

###### Temperature: T opt ######

mod.t.br <- lmer(T.br ~ mic + (1 | block), data = df.summ, REML = TRUE)

summary(mod.t.br)
df.anova.t.br <- as.data.frame(Anova(mod.t.br, type = 3))

df.anova.t.br$term <- rownames(df.anova.t.br)
rownames(df.anova.t.br) <- NULL
df.anova.t.br$trait <- "T.br"

emm.t.br <- emmeans(mod.t.br, ~ mic)
con.t.br <- contrast(emm.t.br, method = "trt.vs.ctrl", ref = "none")

df.mod.t.br <- as.data.frame(con.t.br)
df.mod.t.br$trait <- "T.br"

# Plotting

ctrl.mean <- df.summ %>% # calculate the mean of the control group
  filter(mic == "none") %>%
  summarise(ctrl = mean(T.br, na.rm = TRUE)) %>%
  pull(ctrl)

plot.key <- as.data.frame(con.t.br) %>%
  mutate(
    mic = sub(" - none", "", contrast),
    effect.dir = if_else(estimate > 0, "positive", "negative"),
    sig = p.value < 0.05
  ) %>%
  select(mic, estimate, p.value, effect.dir, sig)

plot.key <- bind_rows(
  tibble(
    mic = "none",
    estimate = 0,
    p.value = NA_real_,
    effect.dir = "control",
    sig = FALSE
  ),
  plot.key
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

p.t.br <- ggplot(df.plot, aes(x = T.br, y = mic)) +
  
  labs(x = "Thermal breadth",  
       y = "Microbe", 
       title = "E — T brt") +  # labels
  
  geom_vline(xintercept = ctrl.mean,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.5) +
  
  ggridges::geom_density_ridges(
    aes(fill = effect.dir,
        alpha = alpha.group),
    scale = 0.9,
    rel_min_height = 0.01,
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_point(position = position_jitter(height = 0.06, width = 0),
             colour = "black",
             size = 1.8,
             alpha = 0.7) +
  
  scale_fill_manual(values = c(
    positive = "dodgerblue3",
    negative = "firebrick",
    control  = "grey70"
  )) +
  
  scale_alpha_manual(values = c(
    sig = 0.81,
    nonsig = 0.15
  )) +
  
  #xlim(1.25, 2.5) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

p.t.br

df.anova.sum <- rbind(df.anova.n.mu, df.anova.n.comp,
                      df.anova.s.mu, df.anova.s.c,
                      df.anova.t.mu, df.anova.t.br,
                      df.anova.t.min, df.anova.t.opt, df.anova.t.max)

write.csv(df.anova.sum, "processed-data/08_data_exp_anova_sum.csv")

df.con.sum <- rbind(df.mod.n.mu, df.mod.n.comp,
                    df.mod.s.mu, df.mod.s.c,
                    df.mod.t.mu, df.mod.t.br,
                    df.mod.t.min, df.mod.t.opt, df.mod.t.max)

write.csv(df.con.sum, "processed-data/09_data_exp_effects_sum.csv")

# Trade-off plots ---------------------------------------------------------

df.ctrl <- df.summ %>%
  filter(mic == "none")

df.all <- df.summ %>%
  filter(mic == "all")

df.mic <- df.summ %>%
  filter(!mic %in% c("none", "all"))

p.n.trade <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = N.µ.max, y = N.comp, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = N.µ.max, y = N.comp, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = N.µ.max, y = N.comp),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = N.µ.max, y = N.comp),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 0.9) +
  
  geom_point(data = df.all,
             aes(x = N.µ.max, y = N.comp),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = N.µ.max, y = N.comp),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 0.9) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (N)")),
       y = expression(paste("Competitive ability, ", 1/N^"*")))

p.n.trade

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
              linewidth = 0.9) +
  
  geom_point(data = df.all,
             aes(x = S.µ.max, y = S.c),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = S.µ.max, y = S.c),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 0.9) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (S)")),
       y = "Salt tolerance")

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
              linewidth = 0.9) +
  
  geom_point(data = df.all,
             aes(x = T.µ.max, y = T.br),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = T.µ.max, y = T.br),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 0.9) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (T)")),
       y = "Thermal breadth")

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
              linewidth = 0.9) +
  
  geom_point(data = df.all,
             aes(x = T.µ.max, y = N.µ.max),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = T.µ.max, y = N.µ.max),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 0.9) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (T)")),
       y = expression(paste("Maximum growth rate, ", mu[max], " (N)")))

p.nt.trade


p.nt2.trade <- ggplot() +
  
  geom_point(data = df.mic,
             aes(x = T.br, y = N.comp, group = mic),
             colour = "black",
             alpha = 0.5,
             size = 2) +
  
  geom_smooth(data = df.mic,
              aes(x = T.br, y = N.comp, group = mic),
              method = "lm",
              se = FALSE,
              colour = "black",
              linewidth = 0.6,
              alpha = 0.5) +
  
  geom_point(data = df.ctrl,
             aes(x = T.br, y = N.comp),
             colour = "darkgreen",
             size = 2.5) +
  
  geom_smooth(data = df.ctrl,
              aes(x = T.br, y = N.comp),
              method = "lm",
              se = FALSE,
              colour = "darkgreen",
              linewidth = 0.9) +
  
  geom_point(data = df.all,
             aes(x = T.br, y = N.comp),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = T.br, y = N.comp),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 0.9) +
  
  theme_classic() +
  labs(x = "1/N*",
       y = "Thermal breadth")

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
              linewidth = 0.9) +
  
  geom_point(data = df.all,
             aes(x = S.µ.max, y = N.µ.max),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = S.µ.max, y = N.µ.max),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 0.9) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (S)")),
       y = expression(paste("Maximum growth rate, ", mu[max], " (N)")))

p.ns.trade

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
              linewidth = 0.9) +
  
  geom_point(data = df.all,
             aes(x = S.µ.max, y = T.µ.max),
             colour = "navyblue",
             size = 2.5) +
  
  geom_smooth(data = df.all,
              aes(x = S.µ.max, y = T.µ.max),
              method = "lm",
              se = FALSE,
              colour = "navyblue",
              linewidth = 0.9) +
  
  theme_classic() +
  labs(x = expression(paste("Maximum growth rate, ", mu[max], " (S)")),
       y = expression(paste("Maximum growth rate, ", mu[max], " (T)")))

p.ts.trade

# Now let's do the treatment shift effect plotting ------------------------

ctrl.means <- df.summ %>%
  filter(mic == "none") %>%
  summarise(
    N.mu.ctrl = mean(N.µ.max, na.rm = TRUE),
    N.comp.ctrl = mean(N.comp, na.rm = TRUE),
    S.mu.ctrl = mean(S.µ.max, na.rm = TRUE),
    S.c.ctrl = mean(S.c, na.rm = TRUE),
    T.mu.ctrl = mean(T.µ.max, na.rm = TRUE),
    T.br.ctrl = mean(T.br, na.rm = TRUE),
    T.min.ctrl = mean(T.min, na.rm = TRUE),
    T.opt.ctrl = mean(T.opt, na.rm = TRUE),
    T.max.ctrl = mean(T.max, na.rm = TRUE)
  )

df.delta <- df.summ %>%
  mutate(
    d.N.mu = N.µ.max - ctrl.means$N.mu.ctrl,
    d.N.comp = N.comp - ctrl.means$N.comp.ctrl,
    d.S.mu = S.µ.max - ctrl.means$S.mu.ctrl,
    d.S.c = S.c - ctrl.means$S.c.ctrl,
    d.T.mu = T.µ.max - ctrl.means$T.mu.ctrl,
    d.T.br = T.br - ctrl.means$T.br.ctrl,
    d.T.min = T.min - ctrl.means$T.min.ctrl,
    d.T.opt = T.opt - ctrl.means$T.opt.ctrl,
    d.T.max = T.max - ctrl.means$T.max.ctrl
  )

df.pca.in <- df.delta %>%
  select(d.N.mu, d.N.comp, d.S.mu, d.S.c, d.T.mu, d.T.br, d.T.min, d.T.opt, d.T.max) %>%
  filter(if_all(everything(), ~ !is.na(.)))

df.pca.scaled <- as.data.frame(scale(df.pca.in))

pca <- prcomp(df.pca.scaled, center = FALSE, scale. = FALSE)

summary(pca)
pca$rotation
pca$sdev^2 / sum(pca$sdev^2)

df.pca.meta <- df.delta %>%
  select(block, mic, rep, id,
         d.N.mu, d.N.comp, d.S.mu, d.S.c, d.T.mu, d.T.br, d.T.min, d.T.opt, d.T.max) %>%
  filter(if_all(starts_with("d."), ~ is.finite(.)))

df.pca.plot <- bind_cols(df.pca.meta, as.data.frame(pca$x))

df.pca.meta <- df.delta %>%
  select(block, mic, rep, id,
         d.N.mu, d.N.comp, d.S.mu, d.S.c, d.T.mu, d.T.br, d.T.min, d.T.opt, d.T.max) %>%
  filter(if_all(starts_with("d."), ~ is.finite(.)))

df.pca.plot <- bind_cols(df.pca.meta, as.data.frame(pca$x))

df.pca.plot <- df.pca.plot %>%
  mutate(mic.col = case_when(
    mic == "none" ~ "control",
    mic == "all"  ~ "all",
    TRUE          ~ "single"
  ))

# variance explained for axis labels
var.expl <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

# loadings / arrows
df.load <- as.data.frame(pca$rotation[, 1:2])
df.load$trait <- rownames(df.load)

# scale arrows so they fit nicely on the score plot
arrow.mult <- 3

p.pca <- ggplot(df.pca.plot, aes(x = PC1, y = PC2)) +
  
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  
  geom_point(aes(colour = mic.col),
             size = 2.2,
             alpha = 0.8) +
  
  geom_segment(data = df.load,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow.mult,
                   yend = PC2 * arrow.mult),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 0.5,
               colour = "black") +
  
  geom_text(data = df.load,
            aes(x = PC1 * arrow.mult,
                y = PC2 * arrow.mult,
                label = trait),
            inherit.aes = FALSE,
            size = 3,
            hjust = 0.5,
            vjust = -0.4) +
  
  scale_colour_manual(values = c(
    control = "forestgreen",
    all = "dodgerblue3",
    single = "black"
  )) +
  
  theme_classic() +
  labs(
    x = paste0("PC1 (", var.expl[1], "%)"),
    y = paste0("PC2 (", var.expl[2], "%)"),
    colour = "Treatment"
  )

p.pca

mic.levels <- levels(factor(df.pca.plot$mic))

# assign colors
mic.cols <- setNames(
  c("forestgreen", "dodgerblue3", hue_pal()(length(mic.levels) - 2)),
  c("none", "all", setdiff(mic.levels, c("none", "all")))
)

scale_colour_manual(values = mic.cols)

p.pca2 <- ggplot(df.pca.plot, aes(x = PC1, y = PC2)) +
  
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  
  geom_point(aes(colour = mic),
             size = 2.2,
             alpha = 0.8) +
  
  geom_segment(data = df.load,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow.mult,
                   yend = PC2 * arrow.mult),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 0.5,
               colour = "black") +
  
  geom_text(data = df.load,
            aes(x = PC1 * arrow.mult,
                y = PC2 * arrow.mult,
                label = trait),
            inherit.aes = FALSE,
            size = 3) +
  
  scale_colour_manual(values = mic.cols) +
  
  theme_classic()

p.pca2

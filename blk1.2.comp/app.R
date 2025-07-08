# J R Laurich

# July 7th, 2025

# A Shiny app to explore the effects of adding carbon to COMBO on TPC, Monod curves, and salt tolerance. 

# Load packages -----------------------------------------------------------

library(shiny)
library(tidyverse)
library(DT)

# Load the data --------------------------------------------------------

setwd("C:/Users/jason/OneDrive/Documents/GitHub/Chlamydomonas_mic_JRL_2024/blk1.2.comp")

df.t <- read.csv('data/5a_blocks1.2_TPC_stats.csv')
head(df.t)

df.n <- read.csv('data/5b_blocks1.2_N_monod_stats.csv')
head(df.n)
df.n$Rstar = 1/df.n$R.mth

df.s <- read.csv('data/5c_blocks1.2_salt_stats.csv')
head(df.s)

# Explore the data ----------------------------------------------------

# So what we will first want to do here is show a table summarizing how different bacteria shape key parameters (both with and without C)

# Temperature, nitrogen and salt - overall stats

df.summ <- bind_cols(
  df.t %>% select(Mic, r.max.raw, T.br.raw, T.min.raw, T.max.raw), 
  df.n %>% select(r.max, R.mth, Rstar),
  df.s %>% select(r.max, c.mod)
)

n <- nrow(df.summ) / 2 # Need to split based on block

df.split <- data.frame(Microbe = df.summ$Mic[1:n])

output_names <- c("µ.max (T)", "T.br", "T.min", "T.max", "µ.max (N)", "R", "R*", "µ.max (S)", "c")
columns_to_split <- c("r.max.raw", "T.br.raw", "T.min.raw", "T.max.raw", "r.max...6", "R.mth", "Rstar", "r.max...9", "c.mod")

for(i in seq_along(columns_to_split)) {
  colname <- columns_to_split[i]
  pretty_name <- output_names[i]
  
  block1 <- df.summ[[colname]][1:n]
  block2 <- df.summ[[colname]][(n+1):(2*n)]
  
  df.split[[paste0(pretty_name, " (B1)")]] <- block1
  df.split[[paste0(pretty_name, " (B2)")]] <- block2
}

head(df.split)

# Temperature - effects of microbes

df.t.D <- df.t %>%
  filter(Mic != "none") %>%  # exclude 'none' for delta calculation rows
  left_join(
    df.t %>%
      filter(Mic == "none") %>%
      select(
        Block, 
        T.br.raw.none = T.br.raw,
        T.opt.raw.none = T.opt.raw,
        r.max.raw.none = r.max.raw
      ),
    by = "Block"
  ) %>%
  mutate(
    D.T.br = T.br.raw - T.br.raw.none,
    D.T.opt = T.opt.raw - T.opt.raw.none,
    D.r.max = r.max.raw - r.max.raw.none
  )

# Nitrogen - effects of microbes

df.n.D <- df.n %>%
  filter(Mic != "none") %>%  # exclude 'none' for delta calculation rows
  left_join(
    df.n %>%
      filter(Mic == "none") %>%
      select(
        Block, 
        K.s.none = K.s,
        r.max.none = r.max,
        R.none = R.mth,
        Rstar.none = Rstar,
      ),
    by = "Block"
  ) %>%
  mutate(
    D.K.s = K.s - K.s.none,
    D.r.max = r.max - r.max.none,
    D.R = R.mth - R.none,
    D.Rstar = Rstar - Rstar.none
  )

# Salt - effects of microbes

df.s.D <- df.s %>%
  filter(Mic != "none") %>%  # exclude 'none' for delta calculation rows
  left_join(
    df.s %>%
      filter(Mic == "none") %>%
      select(
        Block, 
        r.max.none = r.max,
        S.c.none = c.mod
      ),
    by = "Block"
  ) %>%
  mutate(
    D.r.max = r.max - r.max.none,
    D.S.c = c.mod - S.c.none
  )

# The app! ----------------------------------------------------------------

pred_lact <- function(temp, a, b, delta_t, tmax) {
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / delta_t)) + b
}

pred_mon <- function(nit, r.max, k.s) {
  r.max * nit / (k.s + nit)
}

pred_salt <- function(salt, a, b, c) {
  a / (1 + exp(b * (salt - c)))
}

# Define UI ---------------------------------------------------------------
ui <- fluidPage(
  
  # App title
  titlePanel("Effects of carbon addition and microbes on Chlamydomonas ecological traits"),
  
  # Preamble text
  p("This Shiny app allows you to explore how adding carbon to COMBO medium influences thermal performance curves (TPCs), Monod curve parameters (nitrogen limitation), and salt tolerance in Chlamydomonas when interacting with different microbial partners."),
  p("Below is a summary table showing key parameter estimates from both experimental blocks. We added sodium acetate to block 1 but not block 2"),
  
  # Data table output
  dataTableOutput("summary_table"),
  
  # Parameter selection menu
  selectInput("param",
              "Select a parameter to plot:",
              choices = c("µ.max (T)", "T.br", "T.min", "T.max", "µ.max (N)", "R", "R*", "µ.max (S)", "c")),
  
  # Plot output
  plotOutput("param_plot"),
  
  # Section for microbe effects
  hr(),
  h3("Effects of microbial inoculation relative to baseline (none)"),
  
  # Data table output for delta effects
  dataTableOutput("effect_summary_table"),
  
  # Parameter selection menu for effects
  selectInput("param_eff",
              "Select a delta parameter to plot:",
              choices = output_names_effect),
  
  # Plot output for delta effects
  plotOutput("param_eff_plot"),
  
  hr(),
  h3("Plot trait performance curves for a selected microbe"),
  
  # Microbe selection input
  selectInput("microbe_choice",
              "Select a microbial treatment (excluding none):",
              choices = c(as.character(1:15), "all")),
  
  # Plot output
  plotOutput("trait_curves_plot")
  
)

# Define server -----------------------------------------------------------
server <- function(input, output, session) {
  
  # Load the data --------------------------------------------------------
  setwd("C:/Users/jason/OneDrive/Documents/GitHub/Chlamydomonas_mic_JRL_2024/blk1.2.comp")
  
  df.t <- read.csv('data/5a_blocks1.2_TPC_stats.csv')
  df.n <- read.csv('data/5b_blocks1.2_N_monod_stats.csv') %>%
    mutate(Rstar = 1/R.mth)
  df.s <- read.csv('data/5c_blocks1.2_salt_stats.csv')
  
  # Combine and process ----------------------------------------------------
  df.summ <- bind_cols(
    df.t %>% select(Mic, r.max.raw, T.br.raw, T.min.raw, T.max.raw), 
    df.n %>% select(r.max, R.mth, Rstar),
    df.s %>% select(r.max, c.mod)
  )
  
  n <- nrow(df.summ) / 2
  
  df.split <- data.frame(Microbe = df.summ$Mic[1:n])
  
  # Define desired microbe order
  microbe_levels <- c("none", as.character(1:15), "all")
  
  # Relevel Microbe as a factor in df.split
  df.split <- df.split %>%
    mutate(Microbe = factor(Microbe, levels = microbe_levels))
  
  output_names <- c("µ.max (T)", "T.br", "T.min", "T.max", "µ.max (N)", "R", "R*", "µ.max (S)", "c")
  columns_to_split <- c("r.max.raw", "T.br.raw", "T.min.raw", "T.max.raw", "r.max...6", "R.mth", "Rstar", "r.max...9", "c.mod")
  
  for(i in seq_along(columns_to_split)) {
    colname <- columns_to_split[i]
    pretty_name <- output_names[i]
    
    block1 <- df.summ[[colname]][1:n]
    block2 <- df.summ[[colname]][(n+1):(2*n)]
    
    df.split[[paste0(pretty_name, " (B1)")]] <- block1
    df.split[[paste0(pretty_name, " (B2)")]] <- block2
  }
  
  # Render the data table -------------------------------------------------
  output$summary_table <- renderDataTable({
    datatable(df.split, options = list(pageLength = 10))
  })
  
  # Render parameter plot -------------------------------------------------
  
  output$param_plot <- renderPlot({
    req(input$param)
    
    param_b1 <- paste0(input$param, " (B1)")
    param_b2 <- paste0(input$param, " (B2)")
    
    plot_df <- df.split %>%
      select(Microbe, !!param_b1, !!param_b2) %>%
      pivot_longer(cols = c(2,3),
                   names_to = "Block",
                   values_to = "Value") %>%
      mutate(Block = recode(Block,
                            !!param_b1 := "Block 1",
                            !!param_b2 := "Block 2"),
             Block = factor(Block, levels = c("Block 1", "Block 2")))
    
    ggplot(plot_df, aes(x = Microbe, y = Value, colour = Block)) +
      geom_point(size = 3) +
      scale_colour_manual(values = c("Block 1" = "goldenrod", "Block 2" = "magenta")) +
      labs(title = paste("Parameter:", input$param),
           y = input$param,
           x = "Microbe treatment") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Combine microbe effect deltas into one summary table
  df.effect.summ <- bind_cols(
    df.t.D %>% select(Mic, D.r.max, D.T.br, D.T.opt),
    df.n.D %>% select(D.K.s, D.r.max, D.R, D.Rstar),
    df.s.D %>% select(D.r.max, D.S.c)
  )
  
  n.eff <- nrow(df.effect.summ) / 2
  
  df.effect.split <- data.frame(Microbe = df.effect.summ$Mic[1:n.eff])
  
  output_names_effect <- c("Δµ.max (T)", "ΔT.br", "ΔT.opt", "ΔK.s", "Δµ.max (N)", "ΔR", "ΔR*", "Δµ.max (S)", "Δc")
  columns_to_split_effect <- c("D.r.max...2", "D.T.br", "D.T.opt", "D.K.s", "D.r.max...6", "D.R", "D.Rstar", "D.r.max...9", "D.S.c")
  
  for(i in seq_along(columns_to_split_effect)) {
    colname <- columns_to_split_effect[i]
    pretty_name <- output_names_effect[i]
    
    block1 <- df.effect.summ[[colname]][1:n.eff]
    block2 <- df.effect.summ[[colname]][(n.eff+1):(2*n.eff)]
    
    df.effect.split[[paste0(pretty_name, " (B1)")]] <- block1
    df.effect.split[[paste0(pretty_name, " (B2)")]] <- block2
  }
  
  # Reorder Microbe factor levels for consistent plotting
  microbe_levels <- c(as.character(1:15), "all")
  df.effect.split <- df.effect.split %>%
    mutate(Microbe = factor(Microbe, levels = microbe_levels))
  
  # Render microbial effect summary table
  output$effect_summary_table <- renderDataTable({
    datatable(df.effect.split, options = list(pageLength = 10))
  })
  
  # Render microbial effect plot
  output$param_eff_plot <- renderPlot({
    req(input$param_eff)
    
    param_b1 <- paste0(input$param_eff, " (B1)")
    param_b2 <- paste0(input$param_eff, " (B2)")
    
    plot_df <- df.effect.split %>%
      select(Microbe, !!param_b1, !!param_b2) %>%
      pivot_longer(cols = c(2,3),
                   names_to = "Block",
                   values_to = "DeltaValue") %>%
      mutate(Block = recode(Block,
                            !!param_b1 := "Block 1",
                            !!param_b2 := "Block 2"),
             Block = factor(Block, levels = c("Block 1", "Block 2")))
    
    ggplot(plot_df, aes(x = Microbe, y = DeltaValue, colour = Block)) +
      geom_point(size = 3) +
      scale_colour_manual(values = c("Block 1" = "goldenrod", "Block 2" = "magenta")) +
      labs(title = paste("Effect of microbial inoculation on", input$param_eff),
           y = paste("Delta", input$param_eff),
           x = "Microbe treatment") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$trait_curves_plot <- renderPlot({
    req(input$microbe_choice)
    
    # Filter parameter rows for the selected microbe and 'none'
    param_rows_t <- df.t %>%
      filter(Mic %in% c("none", input$microbe_choice))
    
    param_rows_n <- df.n %>%
      filter(Mic %in% c("none", input$microbe_choice))
    
    param_rows_s <- df.s %>%
      filter(Mic %in% c("none", input$microbe_choice))
    
    # Generate predicted temperature curves
    temp_seq <- seq(5, 45, length.out = 200)
    
    pred_t <- param_rows_t %>%
      mutate(LineColour = ifelse(Mic == "none", "green", "red"),
             LineType = ifelse(Block == 1, "solid", "dashed")) %>%
      rowwise() %>%
      do({
        tibble(temp = temp_seq,
               rate = pred_lact(temp, .$cf.a, .$cf.b, .$cf.delta_t, .$cf.tmax),
               LineColour = .$LineColour,
               LineType = .$LineType,
               Mic = .$Mic)
      }) %>%
      ungroup()
    
    p.t <- ggplot(pred_t, aes(x = temp, y = rate, colour = LineColour, linetype = LineType)) +
      geom_line(size = 1) +
      labs(title = "Thermal performance curve", x = "Temperature (°C)", y = "Predicted µmax") +
      scale_colour_identity() +
      theme_classic()
    
    # Generate predicted nitrogen curves
    nitrogen_seq <- seq(0, max(df.n$phos, na.rm = TRUE), length.out = 200)
    
    pred_n <- param_rows_n %>%
      mutate(LineColour = ifelse(Mic == "none", "green", "red"),
             LineType = ifelse(Block == 1, "solid", "dashed")) %>%
      rowwise() %>%
      do({
        tibble(phos = nitrogen_seq,
               rate = pred_mon(phos, .$r.max, .$K.s),
               LineColour = .$LineColour,
               LineType = .$LineType,
               Mic = .$Mic)
      }) %>%
      ungroup()
    
    p.n <- ggplot(pred_n, aes(x = phos, y = rate, colour = LineColour, linetype = LineType)) +
      geom_line(size = 1) +
      labs(title = "Nitrogen Monod curve", x = "Nitrogen (µM)", y = "Predicted µmax") +
      scale_colour_identity() +
      theme_classic()
    
    # Generate predicted salt curves
    salt_seq <- seq(0, max(df.s$salt, na.rm = TRUE), length.out = 200)
    
    pred_s <- param_rows_s %>%
      mutate(LineColour = ifelse(Mic == "none", "green", "red"),
             LineType = ifelse(Block == 1, "solid", "dashed")) %>%
      rowwise() %>%
      do({
        tibble(salt = salt_seq,
               rate = pred_salt(salt, .$a, .$b, .$c),
               LineColour = .$LineColour,
               LineType = .$LineType,
               Mic = .$Mic)
      }) %>%
      ungroup()
    
    p.s <- ggplot(pred_s, aes(x = salt, y = rate, colour = LineColour, linetype = LineType)) +
      geom_line(size = 1) +
      labs(title = "Salt tolerance curve", x = "Salt (g/L)", y = "Predicted µmax") +
      scale_colour_identity() +
      theme_classic()
    
    # Combine plots using patchwork
    library(patchwork)
    
    p.t + p.n + p.s + plot_layout(ncol = 1)
  })
  
  
}

# Run the app -------------------------------------------------------------
shinyApp(ui, server)

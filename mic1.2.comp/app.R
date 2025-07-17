# J R Laurich

# July 14th, 2025

# A Shiny app to explore the effects of adding carbon to COMBO on TPC, Monod curves, and salt tolerance. 

# Load packages -----------------------------------------------------------

library(shiny)
library(tidyverse)
library(DT)
library(patchwork)


# The app -----------------------------------------------------------------

# Define prediction functions
pred_lact <- function(temp, a, b, delta_t, tmax) {
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / delta_t)) + b
}

pred_mon <- function(nit, r.max, k.s) {
  r.max * nit / (k.s + nit)
}

pred_salt <- function(salt, a, b, c) {
  a / (1 + exp(b * (salt - c)))
}

# Define UI
ui <- fluidPage(
  titlePanel("Microbial effects on ecological performance"),
  
  selectInput("microbe_choice", "Select a microbial treatment:",
              choices = c(as.character(1:15), "all")),
  
  plotOutput("trait_curves_plot"),
  
  selectInput("param", "Select a parameter to plot:",
              choices = c("µ.max (T)", "T.br", "T.min", "T.max", "µ.max (N)", "R", "R*", "µ.max (S)", "c")),
  
  plotOutput("param_plot")
)

# Define server
server <- function(input, output, session) {
  setwd("C:/Users/jason/OneDrive/Documents/GitHub/Chlamydomonas_mic_JRL_2024/mic1.2.comp")
  
  df.t.fit <- read.csv('data/7b_block1.2_mic_TPC_fits.csv') %>%
    select(Block, Mic, Chlamy, Parameter, mean) %>%
    pivot_wider(names_from = Parameter, values_from = mean)
  
  df.s.fit <- read.csv('data/7e_blocks1.2_salt_fits.csv') %>%
    select(Block, Mic, Chlamy, Parameter, mean) %>%
    pivot_wider(names_from = Parameter, values_from = mean)
  
  df.n <- read.csv('data/7c_blocks1.2_N_monod_stats.csv')
  
  df.t.mu <- read.csv("data/6a_blk1.2_mic_temp_mus.csv") 
  df.n.mu <- read.csv("data/6b_blk1.2_mic_nit_mus.csv")
  df.s.mu <- read.csv("data/6c_blk1.2_mic_salt_mus.csv") 
  
  output$trait_curves_plot <- renderPlot({
    req(input$microbe_choice)
    
    temp_seq <- seq(5, 45, length.out = 200)
    nitrogen_seq <- seq(0, 1000, length.out = 200)
    salt_seq <- seq(0, 12, length.out = 200)
    
    # Filter all datasets to selected microbe
    param_rows_t <- df.t.fit %>% filter(Mic == input$microbe_choice)
    param_rows_n <- df.n %>% filter(Mic == input$microbe_choice)
    param_rows_s <- df.s.fit %>% filter(Mic == input$microbe_choice)
    
    obs_t <- df.t.mu %>% filter(mic == input$microbe_choice)
    obs_n <- df.n.mu %>% filter(mic == input$microbe_choice)
    obs_s <- df.s.mu %>% filter(mic == input$microbe_choice)
    
    # Temperature plot
    pred_t <- param_rows_t %>%
      mutate(Carbon = ifelse(Block == 1, "Sodium acetate", "No carbon")) %>%
      rowwise() %>%
      do(tibble(temp = temp_seq,
                rate = pred_lact(temp, .$cf.a, .$cf.b, .$cf.delta_t, .$cf.tmax),
                Carbon = .$Carbon,
                Chlamy = .$Chlamy)) %>%
      ungroup()
    
    obs_t <- obs_t %>%
      mutate(Carbon = ifelse(block == 1, "Sodium acetate", "No carbon"))
    
    p.t <- ggplot(pred_t, aes(x = temp, y = rate, colour = Chlamy, linetype = Carbon)) +
      geom_line(size = 1) +
      geom_point(data = obs_t, aes(x = temp, y = r.exp, shape = Carbon, colour = Chlamy),
                 position = position_jitter(width = 0.5), size = 2) +
      labs(title = "Thermal performance curve", x = "Temperature (°C)", y = "µmax") +
      theme_classic()
    
    # Nitrogen plot
    pred_n <- param_rows_n %>%
      mutate(Carbon = ifelse(Block == 1, "Sodium acetate", "No carbon")) %>%
      rowwise() %>%
      do(tibble(nit = nitrogen_seq,
                rate = pred_mon(nit, .$r.max, .$K.s),
                Carbon = .$Carbon,
                Chlamy = .$Chlamy)) %>%
      ungroup()
    
    obs_n <- obs_n %>%
      mutate(Carbon = ifelse(block == 1, "Sodium acetate", "No carbon"))
    
    p.n <- ggplot(pred_n, aes(x = nit, y = rate, colour = Chlamy, linetype = Carbon)) +
      geom_line(size = 1) +
      geom_point(data = obs_n, aes(x = nit, y = r.exp, shape = Carbon, colour = Chlamy),
                 position = position_jitter(width = 8), size = 2) +
      labs(title = "Nitrogen Monod curve", x = "Nitrogen (µM)", y = "µmax") +
      theme_classic()
    
    # Salt plot
    pred_s <- param_rows_s %>%
      mutate(Carbon = ifelse(Block == 1, "Sodium acetate", "No carbon")) %>%
      rowwise() %>%
      do(tibble(salt = salt_seq,
                rate = pred_salt(salt, .$a, .$b, .$c),
                Carbon = .$Carbon,
                Chlamy = .$Chlamy)) %>%
      ungroup()
    
    obs_s <- obs_s %>%
      mutate(Carbon = ifelse(block == 1, "Sodium acetate", "No carbon"))
    
    p.s <- ggplot(pred_s, aes(x = salt, y = rate, colour = Chlamy, linetype = Carbon)) +
      geom_line(size = 1) +
      geom_point(data = obs_s, aes(x = salt, y = r.exp, shape = Carbon, colour = Chlamy),
                 position = position_jitter(width = 0.15), size = 2) +
      labs(title = "Salt tolerance curve", x = "Salt (g/L)", y = "µmax") +
      theme_classic()
    
    p.t + p.n + p.s + plot_layout(nrow = 1)
  })
  
  # Summary parameter plot
  
  df.t <- read.csv('data/7a_blocks1.2_mic_TPC_stats.csv')
  df.n <- read.csv('data/7c_blocks1.2_N_monod_stats.csv') %>%
    mutate(Rstar = 1/R.mth)
  df.s <- read.csv('data/7d_blocks1.2_salt_stats.csv')
  
  df.summ <- df.t %>%
    select(Block, Mic, Chlamy, r.max = r.max.raw, T.br = T.br.raw, T.min = T.min.raw, T.max = T.max.raw) %>%
    full_join(df.n %>% select(Block, Mic, Chlamy, r.max.N = r.max, R = R.mth, Rstar),
              by = c("Block", "Mic", "Chlamy")) %>%
    full_join(df.s %>% select(Block, Mic, Chlamy, r.max.S = r.max, c = c.mod),
              by = c("Block", "Mic", "Chlamy"))

  df.long <- df.summ %>%
    pivot_longer(cols = -c(Block, Mic, Chlamy), names_to = "Parameter", values_to = "Value") %>%
    mutate(Parameter = recode(Parameter,
                              "r.max" = "µ.max (T)", 
                              "T.br" = "T.br", 
                              "T.min" = "T.min", 
                              "T.max" = "T.max",
                              "r.max.N" = "µ.max (N)", 
                              "R" = "R", 
                              "Rstar" = "R*", 
                              "r.max.S" = "µ.max (S)", 
                              "c" = "c"))
  
  
  output$param_plot <- renderPlot({
    req(input$param)
    
    df.long %>%
      filter(Parameter == input$param) %>%
      mutate(ChlamyCol = ifelse(Chlamy == "yes", "blue", "red")) %>%
      ggplot(aes(x = Mic, y = Value)) +
      geom_point(
        aes(shape = factor(Block), 
            fill = ChlamyCol),
        size = 3,
        stroke = 1.5,  # Thickness of outline for open circles
        colour = "black",
        position = position_jitter(width = 0.15)
      ) +
      scale_shape_manual(values = c("1" = 16, "2" = 21)) +
      scale_fill_identity() +  # Use literal colors in `ChlamyCol`
      labs(
        title = paste("Parameter:", input$param),
        x = "Microbe", 
        y = input$param,
        shape = "Block"
      ) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
}

# Run the app
shinyApp(ui, server)
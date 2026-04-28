# Jason R Laurich

# April 19th, 2026

# Let's evaluate the fit of my TPCs (blocked, Bayes)

# Load packages -----------------------------------------------------------

library(shiny)
library(tidyverse)
library(nls.multstart)

# Functions -------------------------------------------------

pred_mon <- function(res, r.max, k.s) {
  r.max * res / (k.s + res)
}

pred_salt <- function(salt, a, b, c) {
  a / (1 + exp(b * (salt - c)))
}

# Load & explore the data --------------------------------------------------------

###### Summary data ######

df.nit <- read.csv('data/11a_chlamy_Monod_nit_bayes.csv') # summary data

head(df.nit)

df.nit <- df.nit %>%
  mutate(mic = factor(mic,
                      levels = c("none", as.character(1:15), "all")),
         block = factor(block))

block.order <- levels(df.nit$block)
mic.order <- levels(df.nit$mic)

df.nit <- df.nit %>%
  arrange(block, mic) # Now this is ordered by block and mic

df.salt <- read.csv('data/12a_chlamy_salt_tol_bayes.csv') # summary data

head(df.salt)

df.salt <- df.salt %>%
  mutate(mic = factor(mic,
                      levels = c("none", as.character(1:15), "all")),
         block = factor(block))


df.salt <- df.salt %>%
  arrange(block, mic) # Now this is ordered by block and mic

###### µ estimates ######

df.mu <- read.csv('data/02_chlamy_µs.csv') # Growth data

head(df.mu)

df.mu <- df.mu %>% # trim the none-temperature data and block 2 (no carbon added)
  filter(temp == 30,
         block != 2)

df.mu <- df.mu %>% # convert block and mic to factor
  mutate(mic = factor(mic,
                      levels = c("none", as.character(1:15), "all")),
         block = factor(block))

df.mu <- df.mu %>%
  arrange(block, mic) # Now this is ordered by block and mic

# UI ----------------------------------------------------------------------

ui <- fluidPage(
  
  titlePanel("Pick a block and microbial treatment"),
  
  fluidRow(
    column(
      3,
      selectInput(
        "block", "Select a block:",
        choices  = block.order,        
        selected = block.order[1],
        width    = "100%"
      ),
      actionButton("prev.block", "Previous block"),
      actionButton("next.block", "Next block")
    )
  ),
  
  fluidRow(
    column(
      3,
      selectInput(
        "mic", "Select a microbial treatment:",
        choices  = mic.order,
        selected = mic.order[1],
        width    = "100%"
      ),
      actionButton("prev.mic", "Previous microbial treatment"),
      actionButton("next.mic", "Next microbial treatment")
    )
  ),
  
  hr(),
  
  fluidRow(
    column(5, plotOutput("p.nit")),
    column(2, ""),
    column(5, plotOutput("p.salt"))
  )
  
)

server <- function(input, output, session) {
  
  ###### Load the information ######
  
  observeEvent(input$next.block, { # Next block
    req(input$block)
    i <- match(input$block, block.order)
    next_i <- ifelse(i == length(block.order), 1, i + 1)
    updateSelectInput(session, "block", selected = block.order[next_i])
  })
  
  observeEvent(input$prev.block, {   # Previous block
    req(input$block)
    i <- match(input$block, block.order)
    prev_i <- ifelse(i == 1, length(block.order), i - 1)
    updateSelectInput(session, "block", selected = block.order[prev_i])
  })
  
  observeEvent(input$next.mic, {   # Next microbial treatment
    req(input$mic)
    j <- match(as.numeric(input$mic), mic.order)
    next_j <- ifelse(j == length(mic.order), 1, j + 1)
    updateSelectInput(session, "mic", selected = mic.order[next_j])
  })
  
  observeEvent(input$prev.mic, {   # Previous microbial treatment
    req(input$mic)
    j <- match(as.numeric(input$mic), mic.order)
    prev_j <- ifelse(j == 1, length(mic.order), j - 1)
    updateSelectInput(session, "mic", selected = mic.order[prev_j])
  })
  
  observeEvent(input$next.rep, { # Next replicate
    req(input$rep)
    i <- match(input$rep, rep.order)
    next_i <- ifelse(i == length(rep.order), 1, i + 1)
    updateSelectInput(session, "rep", selected = rep.order[next_i])
  })
  
  observeEvent(input$prev.rep, {   # Previous replicate
    req(input$rep)
    i <- match(input$rep, rep.order)
    prev_i <- ifelse(i == 1, length(rep.order), i - 1)
    updateSelectInput(session, "rep", selected = rep.order[prev_i])
  })
  
  ###### Generate nit plot ######
  
  output$p.nit <- renderPlot({
    
    mu.n <- df.mu %>% filter(block == input$block,
                             mic == input$mic)
    
    mu.n <- mu.n %>% 
      filter(salt == 0)
    
    df.n <- df.nit %>% filter(block == input$block,
                               mic == input$mic)
    
    curve.n <- tibble::tibble(
      res  = seq(0, 1000, length.out = 1000),
      rate = pred_mon(res, k.s = df.n$K.s.post, r.max = df.n$r.max.post)
    )
    
    ggplot(curve.n, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "royalblue") +
      geom_point(data = mu.n,
                 aes(x = nit, y = µ),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Nitrate concentration (µM)",
        y = "Exponential growth rate",
        title = "Monod curve - nitrogen (blocked, Bayesian)"
      ) +
      theme_classic() +
      
      ylim(-0.5, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_hline(yintercept = df.n$r.max.post, linetype = "dashed", colour = 'forestgreen') +
      geom_vline(xintercept = df.n$K.s.post, linetype = "dashed", colour = 'darkred')
  })
  
  ###### Generate salt plot ######
  
  output$p.salt <- renderPlot({
    
    mu.s <- df.mu %>% filter(block == input$block,
                             mic == input$mic)
    
    mu.s <- mu.s %>% 
      filter(nit == 1000)
    
    df.s <- df.salt %>% filter(block == input$block,
                              mic == input$mic)
    
    curve.s <- tibble::tibble(
      res  = seq(0, 12, length.out = 240),
      rate = pred_salt(salt = res, a = df.s$r.max.post, b = df.s$b.post, c = df.s$c.post)
    )
    
    ggplot(curve.s, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "royalblue") +
      geom_point(data = mu.s,
                 aes(x = salt, y = µ),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Salt concentration (g/L)",
        y = "Exponential growth rate",
        title = "Salt tolerance curve (blocked, Bayesian)"
      ) +
      
      theme_classic() +
      
      ylim(-0.5, 3) +
      
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_hline(yintercept = df.s$r.max.post, linetype = "dashed", colour = 'forestgreen') +
      geom_vline(xintercept = df.s$c.post, linetype = "dashed", colour = 'darkred')
  })
  
}

shinyApp(ui, server)
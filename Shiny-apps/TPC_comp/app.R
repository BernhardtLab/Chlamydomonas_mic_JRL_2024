# Jason R Laurich

# April 20th, 2026

# Let's compare my Bayesian TPC fits (original and with constrained priors)

# Load packages -----------------------------------------------------------

library(shiny)
library(tidyverse)
library(nls.multstart)

# Functions -------------------------------------------------

pred_lact <- function(temp, a, b, d.t, tmax) {
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / d.t)) + b
}

# Load & explore the data --------------------------------------------------------

###### Summary data ######

df.sum1 <- read.csv('data/10a_chlamy_TPCs_bayes.csv') # summary data

head(df.sum1)

df.sum1 <- df.sum1 %>%
  mutate(mic = factor(mic,
                      levels = c("none", as.character(1:15), "all")),
         block = factor(block))

block.order <- levels(df.sum1$block)
mic.order <- levels(df.sum1$mic)

df.sum1 <- df.sum1 %>%
  arrange(block, mic) # Now this is ordered by block and mic

df.sum2 <- read.csv('data/10c_chlamy_TPCs_bayes.csv') # summary data

head(df.sum2)

df.sum2 <- df.sum2 %>%
  mutate(mic = factor(mic,
                      levels = c("none", as.character(1:15), "all")),
         block = factor(block))

block.order <- levels(df.sum2$block)
mic.order <- levels(df.sum2$mic)

df.sum2 <- df.sum2 %>%
  arrange(block, mic) # Now this is ordered by block and mic

###### µ estimates ######

df.mu <- read.csv('data/02_chlamy_µs.csv') # Growth data

head(df.mu)

df.mu <- df.mu %>% # trim the none-temperature data and block 2 (no carbon added)
  filter(salt == 0,
         nit == 1000,
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
    column(5, plotOutput("p.temp1")),
    column(2, ),
    column(5, plotOutput("p.temp2"))
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
  
  ###### Generate TPC plot ######
  
  output$p.temp1 <- renderPlot({
    
    mu.t <- df.mu %>% filter(block == input$block,
                             mic == input$mic)
    
    tpc.t1 <- df.sum1 %>% filter(block == input$block,
                               mic == input$mic)
    
    curve.t1 <- tibble::tibble(
      res  = seq(-5, 45, length.out = 250),
      rate = pred_lact(res, a = tpc.t1$a, b = tpc.t1$b, d.t = tpc.t1$d.t, tmax = tpc.t1$tmax)
    )
    
    ggplot(curve.t1, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.t,
                 aes(x = temp, y = µ),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Lactin II TPC fits (older priors)"
      ) +
      theme_classic() +
      ylim(-0.5, 3) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = tpc.t1$T.min, linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.t1$T.max, linetype = "dashed", colour = 'black') +
      geom_vline(xintercept = tpc.t1$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  output$p.temp2 <- renderPlot({
    
    mu.t <- df.mu %>% filter(block == input$block,
                             mic == input$mic)
    
    tpc.t2 <- df.sum2 %>% filter(block == input$block,
                                 mic == input$mic)
    
    curve.t2 <- tibble::tibble(
      res  = seq(-5, 45, length.out = 250),
      rate = pred_lact(res, a = tpc.t2$a, b = tpc.t2$b, d.t = tpc.t2$d.t, tmax = tpc.t2$tmax)
    )
    
    ggplot(curve.t2, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.t,
                 aes(x = temp, y = µ),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Lactin II TPC fits (newer, narrower priors)"
      ) +
      theme_classic() +
      ylim(-0.5, 3) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = tpc.t2$T.min, linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.t2$T.max, linetype = "dashed", colour = 'black') +
      geom_vline(xintercept = tpc.t2$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
}

shinyApp(ui, server)

# Jason R Laurich

# April 19th, 2026

# Let's evaluate the fit of my TPCs (blocked, Bayes)

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

df.sum <- read.csv('data/10a_chlamy_TPCs_bayes.csv') # summary data

head(df.sum)

df.sum <- df.sum %>%
  mutate(mic = factor(mic,
                      levels = c("none", as.character(1:15), "all")),
         block = factor(block))

block.order <- levels(df.sum$block)
mic.order <- levels(df.sum$mic)

df.sum <- df.sum %>%
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

###### Raw growth rate time series data ######

df.raw <- read.csv('data/01_raw_timeseries_data.csv') # Raw growth data

head(df.raw)

df.raw <- df.raw %>% # rename columns to match other data frames
  rename(block = Block,
         mic = Microbe,
         salt = Salt.conc.g.l,
         nit = Nitrogen.conc.µM,
         temp = Temperature.C,
         rep = Replicate)

df.raw <- df.raw %>% # trim the none-temperature data and block 2 (no carbon added). Also trim the blanks in the plates. 
  filter(salt == 0,
         nit == 1000,
         block != 2,
         mic != "BLANK",
         rep < 4)

df.raw <- df.raw %>% # convert block and mic to factor
  mutate(mic = factor(mic,
                      levels = c("none", as.character(1:15), "all")),
         block = factor(block),
         rep = factor(rep))

df.raw <- df.raw %>%
  arrange(block, mic) # Now this is ordered by block and mic

df.raw <- df.raw %>% 
  filter(Chlamy.y.n == 'y') %>% 
  mutate(log.RFU = log(RFU + 0.001))

rep.order <- levels(df.raw$rep)

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
    column(3, ""),
    column(6, plotOutput("p.temp")),
    column(3, "")
  ),
  
  hr(),
  hr(),
  
  fluidRow(
    column(
      3,
      selectInput(
        "rep", "Select a replicate:",
        choices  = rep.order,
        selected = rep.order[1],
        width    = "100%"
      ),
      actionButton("prev.rep", "Previous replicate"),
      actionButton("next.rep", "Next replicate")
    )
  ),
  
  fluidRow(
    
    fluidRow(
      column(4, plotOutput("p.1")),
      column(4, plotOutput("p.2")),
      column(4, plotOutput("p.3"))
    ),
    
    fluidRow(
      column(4, plotOutput("p.4")),
      column(4, plotOutput("p.5")),
      column(4, plotOutput("p.6"))
    ),
    
    fluidRow(
      column(4, plotOutput("p.7")),
      column(4, plotOutput("p.8")),
      column(4, plotOutput("p.9"))
    )
    
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
  
  output$p.temp <- renderPlot({
    
    mu.t <- df.mu %>% filter(block == input$block,
                             mic == input$mic)
    
    tpc.t <- df.sum %>% filter(block == input$block,
                               mic == input$mic)
    
    curve.t <- tibble::tibble(
      res  = seq(-5, 45, length.out = 250),
      rate = pred_lact(res, a = tpc.t$a, b = tpc.t$b, d.t = tpc.t$d.t, tmax = tpc.t$tmax)
    )
    
    ggplot(curve.t, aes(x = res, y = rate)) +
      geom_line(size = 1.5, colour = "firebrick3") +
      geom_point(data = mu.t,
                 aes(x = temp, y = µ),
                 inherit.aes = FALSE, size = 2) +
      labs(
        x = "Temperature (°C)",
        y = "Exponential growth rate",
        title = "Lactin II TPC fits (blocked, Bayesian)"
      ) +
      theme_classic() +
      ylim(-0.5, 3) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
      
      geom_vline(xintercept = tpc.t$T.min, linetype = "dashed", colour = 'royalblue') +
      geom_vline(xintercept = tpc.t$T.max, linetype = "dashed", colour = 'black') +
      geom_vline(xintercept = tpc.t$T.opt, linetype = "dashed", colour = 'forestgreen')
  })
  
  ###### Generate µ estimation plots ######
  
  output$p.1 <- renderPlot({
    
    df.p1 <- df.raw %>% filter(block == input$block,
                               mic == input$mic,
                               rep == input$rep)
    
    levels <- sort(unique(as.numeric(df.p1$temp)))
    
    df.p1 <- df.p1 %>% filter(temp == levels[1])
    
    df.p1 <- df.p1[order(df.p1$days), ]
    
    if (df.p1$RFU[2] < df.p1$RFU[1]) {
      df.p1 <- df.p1[-1, ]
    } # So if there is a drop in RFUs between the first and second read, we will treat the second data point as N0.
    
    df.p1 <- df.p1 %>% 
      mutate(N0 = RFU[1])
    
    t.series <- unique(df.p1$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series){   # Can't consider the slope just including days 0 
      
      df.p1.sl <- df.p1[df.p1$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(log.RFU~days, data = df.p1.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p1.th <- df.p1[df.p1$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
    
    if (length(unique(na.omit(df.p1.th$RFU))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
      µ.est <- 0
      
    } else {
      
      µ.mod <- nls_multstart(
          RFU ~ N0 * exp(r * days),
          data = df.p1.th,
          start_lower = c(r = -4.5),
          start_upper = c(r = 4.5),
          iter = 500,
          supp_errors = "Y",
          control = nls.control(maxiter = 200)
        )
    
    }
    
    smt.days <- seq(min(df.p1.th$days, na.rm = TRUE), max(df.p1.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(µ.mod, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p1, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Temperature:", df.p1$temp, "°C"), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.2 <- renderPlot({
    
    df.p2 <- df.raw %>% filter(block == input$block,
                               mic == input$mic,
                               rep == input$rep)
    
    levels <- sort(unique(as.numeric(df.p2$temp)))
    
    df.p2 <- df.p2 %>% filter(temp == levels[2])
    
    df.p2 <- df.p2[order(df.p2$days), ]
    
    if (df.p2$RFU[2] < df.p2$RFU[1]) {
      df.p2 <- df.p2[-1, ]
    } # So if there is a drop in RFUs between the first and second read, we will treat the second data point as N0.
    
    df.p2 <- df.p2 %>% 
      mutate(N0 = RFU[1])
    
    t.series <- unique(df.p2$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series){   # Can't consider the slope just including days 0 
      
      df.p2.sl <- df.p2[df.p2$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(log.RFU~days, data = df.p2.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p2.th <- df.p2[df.p2$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
    
    if (length(unique(na.omit(df.p2.th$RFU))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
      µ.est <- 0
      
    } else {
      
      µ.mod <- nls_multstart(
        RFU ~ N0 * exp(r * days),
        data = df.p2.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      )
      
    }
    
    smt.days <- seq(min(df.p2.th$days, na.rm = TRUE), max(df.p2.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(µ.mod, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p2, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Temperature:", df.p2$temp, "°C"), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.3 <- renderPlot({
    
    df.p3 <- df.raw %>% filter(block == input$block,
                               mic == input$mic,
                               rep == input$rep)
    
    levels <- sort(unique(as.numeric(df.p3$temp)))
    
    df.p3 <- df.p3 %>% filter(temp == levels[3])
    
    df.p3 <- df.p3[order(df.p3$days), ]
    
    if (df.p3$RFU[2] < df.p3$RFU[1]) {
      df.p3 <- df.p3[-1, ]
    } # So if there is a drop in RFUs between the first and second read, we will treat the second data point as N0.
    
    df.p3 <- df.p3 %>% 
      mutate(N0 = RFU[1])
    
    t.series <- unique(df.p3$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series){   # Can't consider the slope just including days 0 
      
      df.p3.sl <- df.p3[df.p3$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(log.RFU~days, data = df.p3.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p3.th <- df.p3[df.p3$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
    
    if (length(unique(na.omit(df.p3.th$RFU))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
      µ.est <- 0
      
    } else {
      
      µ.mod <- nls_multstart(
        RFU ~ N0 * exp(r * days),
        data = df.p3.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      )
      
    }
    
    smt.days <- seq(min(df.p3.th$days, na.rm = TRUE), max(df.p3.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(µ.mod, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p3, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Temperature:", df.p3$temp, "°C"), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.4 <- renderPlot({
    
    df.p4 <- df.raw %>% filter(block == input$block,
                               mic == input$mic,
                               rep == input$rep)
    
    levels <- sort(unique(as.numeric(df.p4$temp)))
    
    df.p4 <- df.p4 %>% filter(temp == levels[4])
    
    df.p4 <- df.p4[order(df.p4$days), ]
    
    if (df.p4$RFU[2] < df.p4$RFU[1]) {
      df.p4 <- df.p4[-1, ]
    } # So if there is a drop in RFUs between the first and second read, we will treat the second data point as N0.
    
    df.p4 <- df.p4 %>% 
      mutate(N0 = RFU[1])
    
    t.series <- unique(df.p4$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series){   # Can't consider the slope just including days 0 
      
      df.p4.sl <- df.p4[df.p4$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(log.RFU~days, data = df.p4.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p4.th <- df.p4[df.p4$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
    
    if (length(unique(na.omit(df.p4.th$RFU))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
      µ.est <- 0
      
    } else {
      
      µ.mod <- nls_multstart(
        RFU ~ N0 * exp(r * days),
        data = df.p4.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      )
      
    }
    
    smt.days <- seq(min(df.p4.th$days, na.rm = TRUE), max(df.p4.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(µ.mod, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p4, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Temperature:", df.p4$temp, "°C"), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.5 <- renderPlot({
    
    df.p5 <- df.raw %>% filter(block == input$block,
                               mic == input$mic,
                               rep == input$rep)
    
    levels <- sort(unique(as.numeric(df.p5$temp)))
    
    df.p5 <- df.p5 %>% filter(temp == levels[5])
    
    df.p5 <- df.p5[order(df.p5$days), ]
    
    if (df.p5$RFU[2] < df.p5$RFU[1]) {
      df.p5 <- df.p5[-1, ]
    } # So if there is a drop in RFUs between the first and second read, we will treat the second data point as N0.
    
    df.p5 <- df.p5 %>% 
      mutate(N0 = RFU[1])
    
    t.series <- unique(df.p5$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series){   # Can't consider the slope just including days 0 
      
      df.p5.sl <- df.p5[df.p5$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(log.RFU~days, data = df.p5.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p5.th <- df.p5[df.p5$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
    
    if (length(unique(na.omit(df.p5.th$RFU))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
      µ.est <- 0
      
    } else {
      
      µ.mod <- nls_multstart(
        RFU ~ N0 * exp(r * days),
        data = df.p5.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      )
      
    }
    
    smt.days <- seq(min(df.p5.th$days, na.rm = TRUE), max(df.p5.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(µ.mod, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p5, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Temperature:", df.p5$temp, "°C"), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.6 <- renderPlot({
    
    df.p6 <- df.raw %>% filter(block == input$block,
                               mic == input$mic,
                               rep == input$rep)
    
    levels <- sort(unique(as.numeric(df.p6$temp)))
    
    df.p6 <- df.p6 %>% filter(temp == levels[6])
    
    df.p6 <- df.p6[order(df.p6$days), ]
    
    if (df.p6$RFU[2] < df.p6$RFU[1]) {
      df.p6 <- df.p6[-1, ]
    } # So if there is a drop in RFUs between the first and second read, we will treat the second data point as N0.
    
    df.p6 <- df.p6 %>% 
      mutate(N0 = RFU[1])
    
    t.series <- unique(df.p6$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series){   # Can't consider the slope just including days 0 
      
      df.p6.sl <- df.p6[df.p6$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(log.RFU~days, data = df.p6.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p6.th <- df.p6[df.p6$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
    
    if (length(unique(na.omit(df.p6.th$RFU))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
      µ.est <- 0
      
    } else {
      
      µ.mod <- nls_multstart(
        RFU ~ N0 * exp(r * days),
        data = df.p6.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      )
      
    }
    
    smt.days <- seq(min(df.p6.th$days, na.rm = TRUE), max(df.p6.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(µ.mod, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p6, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Temperature:", df.p6$temp, "°C"), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.7 <- renderPlot({
    
    df.p7 <- df.raw %>% filter(block == input$block,
                               mic == input$mic,
                               rep == input$rep)
    
    levels <- sort(unique(as.numeric(df.p7$temp)))
    
    df.p7 <- df.p7 %>% filter(temp == levels[7])
    
    df.p7 <- df.p7[order(df.p7$days), ]
    
    if (df.p7$RFU[2] < df.p7$RFU[1]) {
      df.p7 <- df.p7[-1, ]
    } # So if there is a drop in RFUs between the first and second read, we will treat the second data point as N0.
    
    df.p7 <- df.p7 %>% 
      mutate(N0 = RFU[1])
    
    t.series <- unique(df.p7$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series){   # Can't consider the slope just including days 0 
      
      df.p7.sl <- df.p7[df.p7$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(log.RFU~days, data = df.p7.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p7.th <- df.p7[df.p7$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
    
    if (length(unique(na.omit(df.p7.th$RFU))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
      µ.est <- 0
      
    } else {
      
      µ.mod <- nls_multstart(
        RFU ~ N0 * exp(r * days),
        data = df.p7.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      )
      
    }
    
    smt.days <- seq(min(df.p7.th$days, na.rm = TRUE), max(df.p7.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(µ.mod, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p7, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Temperature:", df.p7$temp, "°C"), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.8 <- renderPlot({
    
    df.p8 <- df.raw %>% filter(block == input$block,
                               mic == input$mic,
                               rep == input$rep)
    
    levels <- sort(unique(as.numeric(df.p8$temp)))
    
    df.p8 <- df.p8 %>% filter(temp == levels[8])
    
    df.p8 <- df.p8[order(df.p8$days), ]
    
    if (df.p8$RFU[2] < df.p8$RFU[1]) {
      df.p8 <- df.p8[-1, ]
    } # So if there is a drop in RFUs between the first and second read, we will treat the second data point as N0.
    
    df.p8 <- df.p8 %>% 
      mutate(N0 = RFU[1])
    
    t.series <- unique(df.p8$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series){   # Can't consider the slope just including days 0 
      
      df.p8.sl <- df.p8[df.p8$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(log.RFU~days, data = df.p8.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p8.th <- df.p8[df.p8$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
    
    if (length(unique(na.omit(df.p8.th$RFU))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
      µ.est <- 0
      
    } else {
      
      µ.mod <- nls_multstart(
        RFU ~ N0 * exp(r * days),
        data = df.p8.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      )
      
    }
    
    smt.days <- seq(min(df.p8.th$days, na.rm = TRUE), max(df.p8.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(µ.mod, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p8, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Temperature:", df.p8$temp, "°C"), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  output$p.9 <- renderPlot({
    
    df.p9 <- df.raw %>% filter(block == input$block,
                               mic == input$mic,
                               rep == input$rep)
    
    levels <- sort(unique(as.numeric(df.p9$temp)))
    
    df.p9 <- df.p9 %>% filter(temp == levels[9])
    
    df.p9 <- df.p9[order(df.p9$days), ]
    
    if (df.p9$RFU[2] < df.p9$RFU[1]) {
      df.p9 <- df.p9[-1, ]
    } # So if there is a drop in RFUs between the first and second read, we will treat the second data point as N0.
    
    df.p9 <- df.p9 %>% 
      mutate(N0 = RFU[1])
    
    t.series <- unique(df.p9$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
    t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
    
    ln.slopes <- c() # Re-initialize this too!
    
    for (z in t.series){   # Can't consider the slope just including days 0 
      
      df.p9.sl <- df.p9[df.p9$days <= z, ]        # Subset the data to exclude time points above our window
      
      ln_slope <- lm(log.RFU~days, data = df.p9.sl) # calculate the log-linear slope
      
      ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
      
    }
    
    s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
    
    df.p9.th <- df.p9[df.p9$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
    
    if (length(unique(na.omit(df.p9.th$RFU))) == 1) { # If all the numbers in df.i.th are the same. set µ to 0...
      µ.est <- 0
      
    } else {
      
      µ.mod <- nls_multstart(
        RFU ~ N0 * exp(r * days),
        data = df.p9.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      )
      
    }
    
    smt.days <- seq(min(df.p9.th$days, na.rm = TRUE), max(df.p9.th$days, na.rm = TRUE), length.out = 100) # Get a smooth distribution of time points
    
    # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
    fit.RFU.mod <- predict(µ.mod, newdata = data.frame(days = smt.days))
    
    ggplot() + 
      geom_point(data=df.p9, aes(x=days, y=RFU), size = 2) + 
      geom_line(aes(x = smt.days, y = fit.RFU.mod), size = 1) + # Fitted line
      labs(title = paste("Temperature:", df.p9$temp, "°C"), x = "Days", y = "RFU") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
}

shinyApp(ui, server)
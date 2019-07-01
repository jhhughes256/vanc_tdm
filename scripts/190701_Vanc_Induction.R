# Vancomycin PopPK Model Simulation - SA Health Regimen
# ------------------------------------------------------------------------------
# Simulate dosing regimen for vancomycin based on SA Health guidelines
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare workspace
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../vanc_tdm/")

# Load package libraries
  library(dplyr)	# dplyr required for mrgsolve
  library(mrgsolve)	 # Metrum differential equation solver for pharmacometrics
  library(ggplot2)  # graphical visualisation
  # library(MASS)  # mvrnorm

# Source external scripts
  source("scripts/functions_utility.R")  # functions utility
  source("models/vanc_pk.R")  # PopPK model script 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Number of individuals
  nid <- 4*1000  # Number of individuals
  ID <- 1:nid  # Sequence of individual ID's
  
# Source population
  source("scripts/190701_Vanc_Population.R")
  
# Create simulation input dataset
# Define time points
  conc_times <- seq(from = 0, 24*7, by = 0.25)  # 1 day of half hourly data
  freq_bd <- 12*0:14
  freq_d <- c(0, 24*0:6 + 12)
  amt_load <- 2000  # loading dose (mg)
  amt_mnt1 <- 1500  # maintenence dose CrCL = 100 mL/min (mg)
  amt_mnt2 <- 1000  # maintenence dose CrCL = 75 mL/min (mg) 
  amt_mnt3 <- 750  # maintenence dose CrCL = 50 mL/min (mg) 
  amt_mnt4 <- 1000  # maintenence dose CrCL = 25 mL/min (mg)
  
# Create concentration dataset
  input_tb <- tibble::tibble(
      ID = rep(ID, each = length(conc_times)),
      time = rep(conc_times, times = nid),
      amt = 0,
      cmt = 1,
      evid = 0,
      rate = 0
  )
  
# Incorporate dose times (lazy, no tidyverse here)
  dose_cols <- c("amt", "evid", "rate")
  sub_mnt1 <- with(input_tb, ID %in% 1:(nid/4) & time %in% freq_bd)
  input_tb[sub_mnt1, dose_cols] <- data.frame(amt = amt_mnt1, evid = 1, rate = -2)
  
  sub_mnt2 <- with(input_tb, ID %in% (nid/4+1):(nid*2/4) & time %in% freq_bd)
  input_tb[sub_mnt2, dose_cols] <- data.frame(amt = amt_mnt2, evid = 1, rate = -2)
  
  sub_mnt3 <- with(input_tb, ID %in% (nid*2/4+1):(nid*3/4) & time %in% freq_bd)
  input_tb[sub_mnt3, dose_cols] <- data.frame(amt = amt_mnt3, evid = 1, rate = -2)
  
  sub_mnt4 <- with(input_tb, ID %in% (nid*3/4+1):nid & time %in% freq_d)
  input_tb[sub_mnt4, dose_cols] <- data.frame(amt = amt_mnt4, evid = 1, rate = -2)
  
  input_tb[input_tb$time == 0, "amt"] <- amt_load  # loading dose
  
# Create residual unexplained variability dataset
  input_dim <- nid*length(conc_times)
  eps_tb <- mrgsolve::smat(mod, make = TRUE) %>%  # Omega matrix from model
    {MASS::mvrnorm(n = input_dim, mu = rep(0, dim(.)[1]), Sigma = .)} %>%
    tibble::as_tibble() %>%  # matrix used to sample from normal distribution
    dplyr::rename_all(function(x) paste0("EPS", readr::parse_number(x)))
  
# Combine to make final input dataset
  input_tb <- dplyr::bind_cols(input_tb, eps_tb)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate standard dosing regimen
# Pipe dataset to model
  output_tb <- mod %>%
    mrgsolve::data_set(input_tb) %>%  # set input data (observations/events)
    mrgsolve::idata_set(pop_tb) %>%  # set individual data (sets ETAs)
    mrgsolve::carry_out(amt, evid, rate, cmt) %>%  # copy to simulated output
    mrgsolve::mrgsim() %>%  # simulate using mrgsolve
    tibble::as_tibble()
  
# Save data as .RDS for later use
  saveRDS(output_tb, "output/vanc_induction.rds")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot patient data
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Plot individual patient vancomycin concentrations
  p <- NULL
  p <- ggplot(data = output_tb)
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "line", fun.y = median, 
    colour = "red", size = 1)
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "ribbon", 
    fun.ymin = CI90lo,  fun.ymax = CI90hi, fill = "red", size = 1, alpha = 0.25)
  p <- p + labs(x = "Time (hours)", y = "Caffeine Concentration (mg/L)")
  p <- p + coord_cartesian(xlim = c(0, 72), ylim = NULL)
  p <- p + facet_wrap(~CRCL)
  p

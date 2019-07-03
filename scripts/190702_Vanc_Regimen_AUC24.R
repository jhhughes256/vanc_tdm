# Vancomycin Dose-Adjustment Regimen - NCA AUC Proportional Regimen
# -----------------------------------------------------------------------------
# Takes induction dataset and changes the dosage after the time when the first
#   blood sample would be taken. Uses AUC as determined by NCA from the 
#   previous 24 hours to provide proportional dose change to achieve AUC within
#   target range 400 - 700 mg.h/L.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# Load package libraries
  library(dplyr)  # Split and rearrange data - required for mrgsolve
  library(mrgsolve)  # Metrum differential equation solver for pharmacometrics
  library(ggplot2)  # Graphical package

# Source PopPK model script
  source("scripts/functions_utility.R")  # functions utility
  source("models/vanc_pk.R")  # PK model
  
# Read .rds files for population, input and induction datasets
  pop_tb <- readr::read_rds("output/vanc_population.rds")
  input_tb <- readr::read_rds("output/vanc_induction.rds")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  NCAprop_fn <- function(induction_tb) {
  # Change ID2 column name to ID and determine final time
    induction_tb <- dplyr::rename(induction_tb, ID = ID2)
    final_time <- tail(induction_tb, 1) %>% dplyr::pull(time)
  # Determine when to take sample according to exerpt from guidelines:
  # "For patients with normal renal function, the first TDM (trough) should 
  #  occur just prior (within one hour) to the fourth dose or on day 3, 
  #  whichever occurs earlier. In patients with GFR between 20-39 mL/min check 
  #  trough level before (within one hour) of the third dose."
    subj_crcl <- unique(induction_tb$CRCL)
    if (subj_crcl <= 39) {
      if (subj_crcl < 20) {
        stop("Function not designed for patients with CrCL < 20 mL/min.")
      }
      dose_num <- 3
    } else {
      dose_num <- 4
    }
    ncain_tb <- induction_tb
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Begin loop for successive dose adjustment
    repeat {
    # Determine sample time and sample concentration based on the designated
    #   sample number assigned initially, or after dose adjustment.
      sample_time <- dplyr::filter(ncain_tb, evid != 0) %>%
        dplyr::slice(dose_num) %>%
        dplyr::pull(time)
    # Determine previous dose and frequency
      last_amt <- dplyr::filter(ncain_tb, time == sample_time) %>%
        dplyr::pull(amt)
      last_frq <- dplyr::filter(ncain_tb, evid != 0 & time >= sample_time) %>%
        dplyr::pull(time) %>% diff() %>% unique()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Determine what the dose should be for future dose times based on patient 
    #   AUC0-24 prior to the dose adjustment time.
    # First determine the AUC for the previous 24 hours.
      linear_auc24_fn <- function(data) {
      # Prepare function environment
        Ct1 <- head(data$DV, -1)
        Ct2 <- tail(data$DV, -1)
        h <- diff(data$time)
      # Calculate AUC
        auc <- sum((Ct1 + Ct2)*h/2)
      }
      auc24_val <- ncain_tb %>% 
        dplyr::filter(time <= sample_time & time >= (sample_time - 24)) %>%
        linear_auc24_fn()
    # Then determine whether it is within the range. If it, keep the dose the
    #   the same. If it's outside, proportionally change the dose to achieve
    #   500 mg.h/L.
      if (auc24_val >= 400 & auc24_val <= 550) {
        dose_amt <- last_amt
        dose_frq <- last_frq
      } else {
        auc24_target <- 400
        exact_amt <- last_amt*auc24_target/auc24_val
        dose_amt <- ceiling(exact_amt/250)*250  # round to closest 250 mg
    # If dose is not within 0.5 - 1.5 g consider changing dosing frequency
        if (dose_amt >= 500 & dose_amt <= 2000) {
          dose_frq <- last_frq
        } else {
          frq_list <- c(6, 8, 12, 24)
          prev_frq <- last_frq
          repeat{
            frq_num <- which(frq_list == prev_frq) - ifelse(dose_amt < 500, -1, 1)
            if (frq_num %in% 1:length(frq_list)) {
              dose_frq <- frq_list[frq_num]
              exact_amt <- dose_amt*dose_frq/prev_frq
              dose_amt <- ceiling(exact_amt/250)*250  # round to closest 250 mg
            } else {
              dose_frq <- prev_frq
            }
            if (dose_frq == prev_frq | (dose_amt >= 500 & dose_amt <= 2000)) break
            prev_frq <- dose_frq
          }
        }
      }
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Set dose and simulate concentration that result from the dose adjustment.
      if (dose_frq == last_frq & dose_amt == last_amt) {
        mrginput_tb <- ncain_tb %>%
          dplyr::select(ID, time, evid, amt, cmt, rate, EPS1, EPS2)
      } else {
        dose_times <- seq(sample_time, final_time, by = dose_frq)
        adjusted_tb <- dplyr::filter(ncain_tb, time >= sample_time) %>%
          dplyr::mutate(
            evid = dplyr::if_else(time %in% dose_times, 1, 0),
            amt = dplyr::if_else(evid == 1, dose_amt, 0),
            rate = dplyr::if_else(amt != 0, -2, 0)
          )
        mrginput_tb <- dplyr::filter(ncain_tb, time < sample_time) %>%
          dplyr::bind_rows(adjusted_tb) %>%
          dplyr::select(ID, time, evid, amt, cmt, rate, EPS1, EPS2)
      }
    # Simulate new trough concentrations
      ncaout_tb <- mod %>%
        mrgsolve::data_set(mrginput_tb) %>%
        mrgsolve::idata_set(pop_tb) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
        mrgsolve::mrgsim() %>%
        tibble::as_tibble() 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Determine when the next blood sample should be taken. If dose resulted in
    #   an AUC >= 800, then the next sample should in 24 hours. Otherwise, 
    #   next sample should be in 72 hours.
      sample_next <- dplyr::if_else(auc24_val >= 800, 24, 72)
      dose_num <- dose_num + sample_next/dose_frq
      if (dose_num >= dim(dplyr::filter(ncaout_tb, evid != 0))[1]) break
      ncain_tb <- ncaout_tb
    }  # end repeat
    dplyr::select(ncaout_tb, -ID)
  }  # brackets closing "NCAprop_fn"
  
  tictoc::tic()
  output_tb <- input_tb %>%
    # dplyr::filter(ID %in% 1:1000) %>%
  { tibble::add_column(., ID2 = .$ID) } %>%  # so that ID is carried inside of the nest structure
    dplyr::group_by(ID) %>% tidyr::nest() %>%  # create list column for ID data
    dplyr::mutate(bayes = purrr::map(data, NCAprop_fn))  # create new list column using NCA_fn
  tictoc::toc() 
    
  readr::write_rds(output_tb, path = "output/vanc_regimen_auc24.rds")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot patient data
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# Set palette
  cbPalette <- data.frame(
		grey = "#999999",
		orange = "#E69F00",
		skyblue = "#56B4E9",
		green = "#009E73",
		yellow = "#F0E442",
		blue = "#0072B2",
		red = "#D55E00",
		pink = "#CC79A7",
		stringsAsFactors = F
	)
  
# Plot median patient vancomycin concentrations with 90% confidence intervals
  plot_tb <- dplyr::select(output_tb, ID, bayes) %>%
    tidyr::unnest()
  p <- NULL
  p <- ggplot(data = plot_tb)
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "line", fun.y = median,
    colour = "red", size = 1)
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "ribbon",
    fun.ymin = CI90lo,  fun.ymax = CI90hi, fill = "red", size = 1, alpha = 0.25)
  p <- p + labs(x = "Time (hours)", y = "Vancomycin Concentration (mg/L)")
  p <- p + coord_cartesian(xlim = c(0, 168), ylim = NULL)
  p <- p + facet_wrap(~CRCL)
  p
  
# Plot % patients with AUC > 400 mg.h/L
  auc_tb <- dplyr::filter(plot_tb, time %in% (0:7*24)) %>%
    dplyr::group_by(ID) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(data) {
      dplyr::mutate(data, dAUC = c(0, diff(AUC)))
    })) %>% tidyr::unnest() %>%
    dplyr::filter(time != 0) %>%
    dplyr::mutate(time = time/24)
  
  auc_target_bycrcl <- auc_tb %>%
    dplyr::group_by(time, CRCL) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(data) {
      tibble::tibble(
        gt400 = sum(data$dAUC > 400)/dim(data)[1],
        gt700 = sum(data$dAUC > 700)/dim(data)[1]
      )
    })) %>% tidyr::unnest()
  
  p <- NULL
  p <- ggplot(data = auc_target_bycrcl)
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + scale_x_continuous("Time (days)", breaks = 0:7)
  p <- p + scale_y_continuous("IPRE AUC24 > Target AUC24 (%)")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p <- p + facet_wrap(~CRCL)
  p
  
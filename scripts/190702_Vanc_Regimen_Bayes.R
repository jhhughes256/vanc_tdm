# Nivolumab PK Model with Tumour Growth - TDM Step-wise Dosing
# -----------------------------------------------------------------------------
# Simulation of dosing 240 mg every 2 weeks initially, before using TDM with
#   proportional dosage changes.
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
  source("scripts/190320_Nivo_Induction.R")  # induction data
  # sourcing of induction data sources model and population
  # this will be swapped for one script sourcing all scripts for final run
  dose_interval <- 14
  dose_min <- 40
  dose_max <- 800
  dose_opts <- c(40, seq(80, 800, by = 20))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  bayes_fn <- function(induction_tb) {
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
    bayesin_tb <- induction_tb
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Begin loop for successive dose adjustment
    repeat {
      browser()
    # Determine sample time and sample concentration based on the designated
    #   sample number assigned initially, or after dose adjustment.
      sample_time <- dplyr::filter(bayesin_tb, evid != 0) %>%
        dplyr::slice(dose_num) %>%
        dplyr::pull(time)
      sample_conc <- dplyr::filter(bayesin_tb, time == sample_time) %>%
        dplyr::pull(DV)
    # Determine previous dose and frequency
      last_amt <- dplyr::filter(bayesin_tb, time == sample_time) %>%
        dplyr::pull(amt)
      last_frq <- dplyr::filter(bayesin_tb, evid != 0 & time >= sample_time) %>%
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
      auc24_val <- bayesin_tb %>% 
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
        mrginput_tb <- bayesin_tb %>%
          dplyr::select(ID, time, evid, amt, cmt, rate, EPS1, EPS2)
      } else {
        dose_times <- seq(sample_time, final_time, by = dose_frq)
        adjusted_tb <- dplyr::filter(bayesin_tb, time >= sample_time) %>%
          dplyr::mutate(
            evid = dplyr::if_else(time %in% dose_times, 1, 0),
            amt = dplyr::if_else(evid == 1, dose_amt, 0),
            rate = dplyr::if_else(amt != 0, -2, 0)
          )
        mrginput_tb <- dplyr::filter(bayesin_tb, time < sample_time) %>%
          dplyr::bind_rows(adjusted_tb) %>%
          dplyr::select(ID, time, evid, amt, cmt, rate, EPS1, EPS2)
      }
    # Simulate new trough concentrations
      bayesout_tb <- mod %>%
        mrgsolve::data_set(mrginput_tb) %>%
        mrgsolve::idata_set(pop_tb) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
        mrgsolve::mrgsim() %>%
        tibble::as_tibble() 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Determine when the next blood sample should be taken. If dose was held
    #   next blood sample should be in 24 hours. Otherwise, next sample should
    #   be in 72 hours.
      sample_next <- dplyr::if_else(auc24_val >= 800, 24, 72)
      dose_num <- dose_num + sample_next/dose_frq
      if (dose_num >= dim(dplyr::filter(bayesout_tb, evid != 0))[1]) break
      bayesin_tb <- bayesout_tb
    }  # end repeat
    dplyr::select(bayesout_tb, -ID)
  }  # brackets closing "bayes_fn"
  
  tictoc::tic()
  output_tdmstep_df <- trough_flat_df %>%
    dplyr::filter(ID %in% 1:100) %>%
  { tibble::add_column(., ID2 = .$ID) } %>%  # so that ID is carried inside of the nest structure
    dplyr::group_by(ID) %>% tidyr::nest() %>%  # create list column for ID data
    dplyr::mutate(bayes = purrr::map(data, TDMstep_fn))  # create new list column using bayes_fn
  tictoc::toc() 
    
  readr::write_rds(output_tdmstep_df, path = "output/stepwise_tdm.rds")
  
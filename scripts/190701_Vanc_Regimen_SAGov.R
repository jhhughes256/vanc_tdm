# Vancomycin Dose-Adjustment Regimen - SA Health Clinical Guidelines
# -----------------------------------------------------------------------------
# Takes induction dataset and changes the dosage after the time when the first
#   blood sample would be taken.
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
  SAHealth_fn <- function(induction_tb) {
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
    sahin_tb <- induction_tb
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Begin loop for successive dose adjustment
    repeat {
    # Determine sample time and sample concentration based on the designated
    #   sample number assigned initially, or after dose adjustment.
      sample_time <- dplyr::filter(sahin_tb, evid != 0) %>%
        dplyr::slice(dose_num) %>%
        dplyr::pull(time)
      sample_conc <- dplyr::filter(sahin_tb, time == sample_time) %>%
        dplyr::pull(DV)
    # Determine previous dose and frequency
      last_amt <- dplyr::filter(sahin_tb, time == sample_time) %>%
        dplyr::pull(amt)
      last_frq <- dplyr::filter(sahin_tb, evid != 0 & time >= sample_time) %>%
        dplyr::pull(time) %>% diff() %>% unique()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Determine what the dose should be for future dose times based on patient 
    #   GFR and trough concentration.
    # If trough concentration more than 30 mg/L hold dose regardless of GFR
      if (sample_conc > 30) {
        dose_amt <- 0  # mg
        dose_frq <- 24  # hourly
    # GFR > 90
      } else if (subj_crcl > 90) {
        if (sample_conc < 10) {  # < 10 mg/L; 1 g 6-hrly
          dose_amt <- 1000  # mg
          dose_frq <- 6  # hourly
        } else if (sample_conc < 15) {  # 10 - 14.9 mg/L; 1.25 g 8-hrly
          dose_amt <- 1250  # mg
          dose_frq <- 8  # hourly
        } else if (sample_conc <= 20) {  # 15 - 20 mg/L; no change
          dose_amt <- dplyr::if_else(last_amt != 0, last_amt, 1000)  # mg
          dose_frq <- dplyr::if_else(last_amt != 0, last_frq, 12)  # hourly
          if (last_amt == 0) browser()
        } else if (sample_conc <= 25) {  # 21 - 25 mg/L; 1.25 g 12-hrly
          dose_amt <- 1250  # mg
          dose_frq <- 12  # hourly
        } else if (sample_conc <= 30) {  # 26 - 30 mg/L; 1 g 12-hrly
          dose_amt <- 1000  # mg
          dose_frq <- 12  # hourly
        } else stop("Critical Error")  # capture errors (sample_conc == NA)
    # GFR 60-90
      } else if (subj_crcl >= 60) {
        if (sample_conc < 10) {  # < 10 mg/L; 1.5 g 12-hrly
          dose_amt <- 1500  # mg
          dose_frq <- 12  # hourly
        } else if (sample_conc < 15) {  # 10 - 14.9 mg/L; 1.25 g 12-hrly
          dose_amt <- 1250  # mg
          dose_frq <- 12  # hourly
        } else if (sample_conc <= 20) {  # 15 - 20 mg/L; no change
          dose_amt <- dplyr::if_else(last_amt != 0, last_amt, 750)  # mg
          dose_frq <- dplyr::if_else(last_amt != 0, last_frq, 12)  # hourly
          if (last_amt == 0) browser()
        } else if (sample_conc <= 30) {  # 21 - 30 mg/L; 750 mg 12-hrly
          dose_amt <- 750  # mg
          dose_frq <- 12  # hourly
        } else stop("Critical Error")  # capture errors (sample_conc == NA)
    # GFR 40-59
      } else if (subj_crcl >= 40) {
        if (sample_conc < 15) {  # < 15 mg/L; 1 g 12-hrly
          dose_amt <- 1000  # mg
          dose_frq <- 12  # hourly
        } else if (sample_conc <= 20) {  # 15 - 20 mg/L; no change
          dose_amt <- dplyr::if_else(last_amt != 0, last_amt, 500)  # mg
          dose_frq <- dplyr::if_else(last_amt != 0, last_frq, 12)  # hourly
          if (last_amt == 0) browser()
        } else if (sample_conc <= 30) {  # 21 - 30 mg/L; 500 mg 12-hrly
          dose_amt <- 500  # mg
          dose_frq <- 12  # hourly
        } else stop("Critical Error")  # capture errors (sample_conc == NA)
    # GFR 20-39
      } else if (subj_crcl >= 20) {
        if (sample_conc < 15) {  # < 15 mg/L; 750 mg 12-hrly
          dose_amt <- 750  # mg
          dose_frq <- 12  # hourly
          warning("
            Subject with GFR 20-39 & trough < 15 mg/L; algorithm does not alter 
            magnitude of dose dependent on trough levels and GFR nor does it 
            consider changing to continuous infusion.
          ")
        } else if (sample_conc <= 20) {  # 15 - 20 mg/L; no change
          dose_amt <- dplyr::if_else(last_amt != 0, last_amt, 750)  # mg
          dose_frq <- dplyr::if_else(last_amt != 0, last_frq, 24)  # hourly
          if (last_amt == 0) browser()
        } else if (sample_conc <= 30) {  # 21 - 30 mg/L; 750 mg 24-hrly
          dose_amt <- 750  # mg
          dose_frq <- 24  # hourly
        } else stop("Critical Error")  # capture errors (sample_conc == NA)
      }  # end if else
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Set dose and dose frequency and simulate concentration that result from
    #   the dose adjustment.
      if (dose_frq == last_frq & dose_amt == last_amt) {
        mrginput_tb <- sahin_tb %>%
          dplyr::select(ID, time, evid, amt, cmt, rate, EPS1, EPS2)
      } else {
        dose_times <- seq(sample_time, final_time, by = dose_frq)
        adjusted_tb <- dplyr::filter(sahin_tb, time >= sample_time) %>%
          dplyr::mutate(
            evid = dplyr::if_else(time %in% dose_times, 1, 0),
            amt = dplyr::if_else(evid == 1, dose_amt, 0),
            rate = dplyr::if_else(amt != 0, -2, 0)
          )
        mrginput_tb <- dplyr::filter(sahin_tb, time < sample_time) %>%
          dplyr::bind_rows(adjusted_tb) %>%
          dplyr::select(ID, time, evid, amt, cmt, rate, EPS1, EPS2)
      }
    # Simulate new trough concentrations
      sahout_tb <- mod %>%
        mrgsolve::data_set(mrginput_tb) %>%
        mrgsolve::idata_set(pop_tb) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
        mrgsolve::mrgsim() %>%
        tibble::as_tibble() 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Determine when the next blood sample should be taken. If dose was held
    #   next blood sample should be in 24 hours. Otherwise, next sample should
    #   be in 72 hours.
      sample_next <- dplyr::if_else(dose_amt == 0, 24, 72)
      dose_num <- dose_num + sample_next/dose_frq
      if (dose_num >= dim(dplyr::filter(sahout_tb, evid != 0))[1]) break
      sahin_tb <- sahout_tb
    }  # end repeat
    sahout_tb
  }  # brackets closing "bayes_fn"
  
  tictoc::tic()
  output_tb <- input_tb %>%
    # dplyr::filter(ID %in% 1:100) %>%
  { tibble::add_column(., ID2 = .$ID) } %>%  # so that ID is carried inside of the nest structure
    dplyr::group_by(ID) %>% tidyr::nest() %>%  # create list column for ID data
    dplyr::mutate(bayes = purrr::map(data, SAHealth_fn))  # create new list column using SA_fn
  tictoc::toc() 
    
  readr::write_rds(output_tdmstep_df, path = "output/vanc_regimen_sagov.rds")
  
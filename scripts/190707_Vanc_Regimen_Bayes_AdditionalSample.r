# Vancomycin Dose-Adjustment Regimen - Bayesian Optimised Dosing Regimen
# -----------------------------------------------------------------------------
# Takes induction dataset and changes the dosage after the time when the first
#   blood sample would be taken. Uses empirical bayesian estimation to 
#   determine patient pharmacokinetic characteristics. This is then used to
#   determine the best dose for the patient.
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
    dose_num <- 2
    if (subj_crcl <= 39) {
      if (subj_crcl < 20) {
        stop("Function not designed for patients with CrCL < 20 mL/min.")
      }
      mnt_frq <- 24
    } else {
      mnt_frq <- 12
    }
    bayesin_tb <- induction_tb
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Begin loop for successive dose adjustment
    sample_all <- c(NULL)
    repeat {
    # Determine sample time and sample concentration based on the designated
    #   sample number assigned initially, or after dose adjustment.
      sample_time <- dplyr::filter(bayesin_tb, evid != 0) %>%
        dplyr::slice(dose_num) %>%
        dplyr::pull(time)
      sample_all <- c(sample_all, sample_time)
      # sample_conc <- dplyr::filter(bayesin_tb, time == sample_time) %>%
      #   dplyr::pull(DV)
    # Determine previous dose and frequency
      last_amt <- dplyr::filter(bayesin_tb, time == sample_time) %>%
        dplyr::pull(amt)
      last_frq <- dplyr::filter(bayesin_tb, evid != 0 & time >= sample_time) %>%
        dplyr::pull(time) %>% diff() %>% unique()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Empirical Bayesian Estimation
    # Loop until successful minimisation
      run_num <- 0
      sample_tb <- bayesin_tb %>%
        dplyr::filter(evid == 1) %>%
        dplyr::mutate(DV = dplyr::if_else(time %in% sample_all, DV, NA_real_))
      repeat {
      # Initial estimates for Bayes parameters
        ETABSV <- mrgsolve::omat(mod, make = TRUE) %>% diag()
        n_eta <- length(ETABSV)
        mod_etas <- 1:n_eta
        init_par <- exp(double(n_eta))
        if (run_num > 0) {
          init_par <- init_par*exp(runif(n_eta, min = -0.01, max = 0.01))
        }
      # Previous dependent variable values
        prev_DV <- dplyr::filter(bayesin_tb, time %in% sample_all) %>% 
          dplyr::pull(DV)
      # Define bayesian estimation function
        bayes_estimate_fn <- function(par) {
        # Describe parameters to be optimised within mrgsolve data set
          ETA <- log(par)
        # Define mrgsolve dataset
          estim_tb <- bayesin_tb %>%
            dplyr::select(ID, time, evid, amt, cmt, rate, WT, DIAL, CRCL) %>%
            dplyr::mutate_at(paste0("ETA", mod_etas), function(x) {
              eta <- as.numeric(substr(deparse(substitute(x)), 4, 4))
              which_eta <- which(mod_etas == eta)
              ETA[which_eta]
            }) %>%
          # Run data through mrgsolve, with idata using initial tumour size
            mrgsolve::data_set(x = mod, data = .) %>%
            mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
            mrgsolve::mrgsim() %>%
            tibble::as_tibble() %>% 
          # Ensure IPRED has finite values
            dplyr::mutate(IPRE = dplyr::if_else(
              !is.finite(IPRE) | IPRE < .Machine$double.eps,
              .Machine$double.eps,
              IPRE
            ))
        # Define yhat
          yhat <- dplyr::filter(estim_tb, time %in% sample_all) %>%
            dplyr::pull(IPRE)
        # Posterior log-likelihood
        # Error model: IPRE*(1 + ERR_PRO) + ERR_ADD
        # Can be simplified to: IPRE + W*ERR
        # Where W = sqrt(pow(IPRE*ERR_PRO, 2) + pow(ERR_ADD, 2))
          loglikpost_sd <- sqrt((yhat*mod$ERR_PRO)^2 + mod$ERR_ADD^2)
          loglikpost <- dnorm(prev_DV, mean = yhat, sd = loglikpost_sd, log = T)
        # Prior log-likelihood
          loglikprior <- dnorm(ETA, mean = 0, sd = sqrt(ETABSV[mod_etas]), log = T)
        # Return objective function value to be minimised
          return(-1*sum(loglikpost, loglikprior))
        }  # end bayes_estimate_fn
      # Run bayes_estimate_fn using optim()
        bayes_estimate <- try(optim(init_par, bayes_estimate_fn, method = "L-BFGS-B",
          lower = rep(0.001, times = n_eta), upper = rep(1000, times = n_eta),
          control = list(
            parscale = init_par, fnscale = bayes_estimate_fn(init_par),
            factr = 1e12, pgtol = 1e-8
          )
        ))
        if (class(bayes_estimate) == "try-error") browser()  # error catch
        minimised <- "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
        if (bayes_estimate$message == minimised) break
        run_num <- run_num + 1
      }  # end loop as successful minimisation has occurred
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculate concentrations according to new Bayes estimates.
      input_sim_tb <- bayesin_tb %>%
        dplyr::select(ID, time, evid, amt, cmt, rate, WT, DIAL, CRCL) %>%
        dplyr::mutate_at(paste0("ETA", mod_etas), function(x) {
          eta <- as.numeric(substr(deparse(substitute(x)), 4, 4))
          which_eta <- which(mod_etas == eta)
          log(bayes_estimate$par[which_eta])
        })
      output_sim_tb <- input_sim_tb %>%
        dplyr::filter(time <= sample_time) %>%
        mrgsolve::data_set(x = mod, data = .) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
        mrgsolve::mrgsim() %>%
        tibble::as_tibble()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Optimise dose for the individual using individual Bayes predicted 
    #   concentrations (and compartment amounts) at time of last sampling.
    #   Only optimise if the AUC0-24 for the last day is outside of the target
    #   range.
    # Determine AUC over the last 24 hours then define when next sample occurs
      last_bayes_auc <- output_sim_tb %>%
        dplyr::filter(time %in% c(sample_time - 24, sample_time)) %>%
        dplyr::pull(AUC)
      if (length(last_bayes_auc) != 1)  last_bayes_auc <- diff(last_bayes_auc)
      else last_bayes_auc <- last_bayes_auc*2
    # If optimising before maintenance dose
      if(!(exists("mnt_flag")) | last_bayes_auc >= 800) {
        sample_next <- 24
      } else {
        sample_next <- 72
      }
      next_time <- sample_time + sample_next
      if (next_time > final_time) next_time <- final_time
    # Define trough target and upper bound, begin if statement
      auc_lower <- 400
      auc_target <- 500
      auc_upper <- 550
      if (last_bayes_auc < auc_lower | last_bayes_auc >= auc_upper | last_amt == 0) {
      # Modify model code ready for simulation
        optim_mod <- mrgsolve::init(mod, list(
          CENT = dplyr::filter(output_sim_tb, time == sample_time) %>% 
            dplyr::pull(CENT),
          PERI = dplyr::filter(output_sim_tb, time == sample_time) %>% 
            dplyr::pull(PERI),
          AUC = dplyr::filter(output_sim_tb, time == sample_time) %>% 
            dplyr::pull(AUC)
        ))
      # Set initial dose and error estimates
        init_par <- c(1000, 0.01)
        if (dose_num == 2) {
          prev_frq <- mnt_frq
        } else {
          prev_frq <- dplyr::if_else(last_amt == 0, 12, last_frq)
        }
        repeat{
        # Subset input dataset so only future concentrations are predicted
          test_times <- seq(sample_time, final_time, by = prev_frq)
          input_optim_tb <- input_sim_tb %>% 
            dplyr::filter(time >= sample_time & time <= next_time) %>%
            dplyr::mutate(
              evid = dplyr::if_else(time %in% test_times, 1, 0),
              rate = dplyr::if_else(evid != 0, -2, 0)
            )
        # Find the doses that maximise the likelihood of trough concentrations
        #   being the target concentration
          optimise_dose_fn <- function(par) {
          # Add fitted parameters to the input data set, then simulate 
          #   concentration-time profiles with fitted doses
            output_optim_tb <- input_optim_tb %>%
              dplyr::mutate(amt = dplyr::if_else(
                evid == 1, par[1], amt
              )) %>%
              mrgsolve::data_set(x = optim_mod, data = .) %>%
              mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
              mrgsolve::mrgsim(start = sample_time, end = next_time) %>%
              tibble::as_tibble() %>%
              dplyr::mutate(IPRE = dplyr::if_else(
                !is.finite(IPRE) | IPRE < .Machine$double.eps,
                .Machine$double.eps,
                IPRE
              ))
          # Define yhat and the residual
            yhat <- output_optim_tb %>%
              dplyr::filter(time %in% c(next_time - 24, next_time)) %>%
              dplyr::pull(AUC) %>% diff()
            res <- dnorm(auc_target, yhat, yhat*par[2], log = T)
          # Objective function value to be minimised
            return(-1*sum(res))
          }
          optimise_dose <- try(optim(init_par, optimise_dose_fn, method = "L-BFGS-B",
            lower = c(0.0001, 0.0001),
            upper = c(5000, Inf),
            control = list(parscale = init_par, factr = 1e7)
          ))
          if (class(optimise_dose) == "try-error") browser()  # error catch
        # Administer the individual the optimised dose
          exact_amt <- optimise_dose$par[1]
          if (exact_amt >= 500 | exact_amt <= 2000) {
            dose_frq <- prev_frq
            if (exists("pls_stop_later")) browser()  # exception
            break
          }
          frq_list <- c(6, 8, 12, 24)
          frq_num <- which(frq_list == prev_frq) - ifelse(exact_amt < 500, -1, 1)
          prev_frq <- frq_list[frq_num]
          pls_stop_later <- 1  # test
          browser()  # exception
        }  # end repeat
        if (exact_amt < 1) exact_amt <- 0
        dose_amt <- ceiling(exact_amt/250)*250
      } else {
      # Give previous dose
        dose_amt <- last_amt
        dose_frq <- last_frq
      }  # end if at target auc
    # Check to see if dose change is required
      if (dose_frq == last_frq & dose_amt == last_amt) {
        input_final_tb <- bayesin_tb %>%
          dplyr::select(ID, time, evid, amt, cmt, rate, EPS1, EPS2)
      } else {
        dose_times <- seq(sample_time, final_time, by = dose_frq)
        adjusted_tb <- dplyr::filter(bayesin_tb, time >= sample_time) %>%
          dplyr::mutate(
            evid = dplyr::if_else(time %in% dose_times, 1, 0),
            amt = dplyr::if_else(evid == 1, dose_amt, 0),
            rate = dplyr::if_else(amt != 0, -2, 0)
          )
        input_final_tb <- dplyr::filter(bayesin_tb, time < sample_time) %>%
          dplyr::bind_rows(adjusted_tb) %>%
          dplyr::select(ID, time, evid, amt, cmt, rate, EPS1, EPS2)
      }
    # Simulate to represent time passing since last sample
      bayesout_tb <- mod %>%
        mrgsolve::data_set(data = input_final_tb) %>%
        mrgsolve::idata_set(data = pop_tb) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt) %>%
        mrgsolve::mrgsim() %>%
        tibble::as_tibble()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Determine when the next blood sample should be taken. If dose was held
    #   next blood sample should be in 24 hours. Otherwise, next sample should
    #   be in 72 hours.
    # This is done above as it impacts dosing optimisation
      # sample_next <- dplyr::if_else(auc24_val >= 800, 24, 72)
      if (dose_num == 2) {
        mnt_flag <- 1
        if (subj_crcl <= 39) {
          dose_num <- 3
        } else {
          dose_num <- 4
        }
      } else {
        dose_num <- dose_num + sample_next/dose_frq
      }
      if (dose_num >= dim(dplyr::filter(bayesout_tb, evid != 0))[1]) break
      bayesin_tb <- bayesout_tb
    }  # end repeat
    dplyr::select(bayesout_tb, -ID)
  }  # brackets closing "bayes_fn"
  
  tictoc::tic()
  output_tb <- input_tb %>%
    # dplyr::filter(ID %in% c(3001:3100)) %>%
  { tibble::add_column(., ID2 = .$ID) } %>%  # so that ID is carried inside of the nest structure
    dplyr::group_by(ID) %>% tidyr::nest() %>%  # create list column for ID data
    dplyr::mutate(bayes = purrr::map(data, bayes_fn))  # create new list column using bayes_fn
  tictoc::toc() 
    
  readr::write_rds(output_tb, path = "output/vanc_regimen_bayes_addsamp.rds")
  
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
  
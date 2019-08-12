# Vancomytin PopPK Model Simulation - Population
# ------------------------------------------------------------------------------
# Create population of representative patients for simulations
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define demographic data
  set.seed(123456)

# Reference population
# Weight - 79 kg (?) [33 - 255]
# Creatinine Clearance - 62 (?) [4 - 150]
# Haemodialysis - 18.5% on dialysis

# Simulated population
# Weight - 75 kg (fixed)
# Creatinine Clearance - 100, 75, 50, 25 ml/min (four fixed groups)
# Haemodialysis - 0% (all not on dialysis)
  cov_tb <- tibble::tibble(
    WT = rep(75, times = nid),
    CRCL = rep(c(100, 75, 50, 25), each = nid/4),
    DIAL = rep(0, times = nid)
  )
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Population parameter variability
# Extract omega block values from model and allocate individual distributions
  eta_tb <- mrgsolve::omat(mod, make = TRUE) %>%
    {MASS::mvrnorm(n = nid, mu = rep(0, dim(.)[1]), Sigma = .)} %>%
    tibble::as_tibble() %>% 
    dplyr::rename_all(function(x) paste0("ETA", readr::parse_number(x)))
  
# Create data frame of individuals with varying demographics and ETA values
  pop_tb <- dplyr::bind_cols(ID = ID, cov_tb, eta_tb)
  
# Vancomycin Population Pharmacokinetic Model - Goti et al.
# ------------------------------------------------------------------------------
# Model sourced from:
# Goti V, Chaturvedula A, Fossler MJ, Mok S, Jacob JT. Hospitalized Patients 
#   With and Without Hemodialysis Have Markedly Different Vancomycin 
#   Pharmacokinetics: A Population Pharmacokinetic Model-Based Analysis. Ther 
#   Drug Monit. 2018;40(2):212-221.
#
# Infusion duration as per "Clinical Practice Guideline for Dosing and 
#   Monitoring of Vancomycin in Adults" by SA Health of the Government of
#   South Australia.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load libraries
  # library(dplyr)
  # library(mrgsolve)

# Define model code
  code <- '
$INIT  // Initial Conditions for Compartments
  CENT = 0,  // Central compartment
  PERI = 0,  // Peripheral compartment
  AUC = 0,  // Area under the curve

$SET  // Set Differential Equation Solver Options			
  atol = 1e-8, rtol = 1e-8
  maxsteps = 100000

$PARAM  // Population parameters
  TVCL = 4.5,  // clearance (L/h)
  TVVC = 58.4,  // central volume of distribution (L)
  TVVP = 38.4,  // peripheral volume of distribution (L)
  TVQ = 6.5,  // inter-compartmental clearance (L/h)

  AVBWT = 70,  // average body weight (kg)
  AVCRCL = 120,  // average creatinine clearance (mL/min)

  // Covariate Effects
  CL_CRCL = 0.8,  // creatinine clearance on clearance
  CL_DIAL = 0.7,  // haemodialysis on clearance
  VC_DIAL = 0.5,  // haemodialysis on volume of distribution

  // Residual Variability
  ERR_PRO = 0.227,  // proportional error (fraction)
  ERR_ADD = 3.4,  // additive error (mg/L)

  // Default Covariate Values for Simulation
  WT = 79,  // body weight (kg)
  CRCL = 62,  // creatinine clearance (mL/min)
  DIAL = 0,  // haemodialysis status (0 - not on dialysis; 1 - on dialysis)

  // Default ETA Values for Simulation
  // Allocated in population so set to zero
  ETA1 = 0,  // PPVCL
  ETA2 = 0,  // PPVVC
  ETA3 = 0,  // PPVVP

  // Default EPS values for simulation
  // Allocated in population so set to zero
  EPS1 = 0,  // RUVPRO
  EPS2 = 0,  // RUVADD

$OMEGA  // Population parameter Variability
  name = "omega1"
  block = FALSE
  0.158404  // PPVCL
  0.665856  // PPVVC
  0.326041  // PPVVP

$SIGMA  // Residual Unexplained Variability	
  block = FALSE
  1  // RUVPRO
  1  // RUVADD

$PREAMBLE  // Set up C++ Environment to Capture Last Dose
  double last_dose = 0;

$MAIN  // Covariate Effects
  double COVCL = TVCL*pow(CRCL/120, CL_CRCL)*pow(CL_DIAL, DIAL);
  double COVVC = TVVC*pow(WT/70, 1)*pow(VC_DIAL, DIAL);

  // Individual Parameter Values
  double CL = COVCL*exp(ETA1);  // *exp(ETA(1))
  double VC = COVVC*exp(ETA2);  // *exp(ETA(2))
  double VP = TVVP*exp(ETA3);  // *exp(ETA(3))
  double Q = TVQ;

  // Infusion Duration (hours)
  if(EVID == 1) last_dose = self.amt;
  if(ceil(last_dose/1000) <= 1) {
    D_CENT = 60/60;
  } else if(ceil(last_dose/1500) <= 1) {
    D_CENT = 90/60;
  } else if(ceil(last_dose/2000) <= 1) {
    D_CENT = 120/60;
  } else {
    D_CENT = self.amt/1000;
  }

$ODE  // Differential Equations
  double C1 = CENT/VC;
  double C2 = PERI/VP;

  dxdt_CENT = C2*Q - C1*Q - C1*CL;
  dxdt_PERI = C1*Q - C2*Q;
  dxdt_AUC = C1;

$TABLE  // Determines Values and Includes in Output	
  double IPRE = C1;  // individual prediction
  double DV = IPRE*(1 + EPS1*ERR_PRO) + EPS2*ERR_ADD;  // observed concentration
  // double DV = IPRE*(1 + EPS(1)*ERR_PRO) + EPS(2)*ERR_ADD;  // debug

$CAPTURE 
  WT CRCL DIAL IPRE DV C2 CENT PERI AUC
  CL VC VP Q COVCL COVVC
  ETA1 ETA2 ETA3 EPS1 EPS2
'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile the model code
  mod <- mrgsolve::mcode("vancPK", code)
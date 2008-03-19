// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstring>
#include <cmath>
#include "iceModel.hh"


const PetscScalar DEFAULT_START_YEAR = 0;
const PetscScalar DEFAULT_RUN_YEARS = 1000.0;  // years

//used in iMutil.C
const PetscScalar DEFAULT_ADAPT_TIMESTEP_RATIO = 0.12;  // yes, I'm confident this is o.k.

//used in iMvelocity.C
const PetscScalar DEFAULT_MINH_SSA = 5.0;  // m; minimum thickness for SSA velocity computation
const PetscScalar DEFAULT_MIN_SHEET_TO_DRAGGING = 50.0;   // m/a; critical SIA speed for switch SIA --> SSA
const PetscScalar DEFAULT_MAX_SPEED_DRAGGING_TO_SHEET = 5.0;  // m/a; crit Mac speed for switch SSA --> SIA
const PetscScalar DEFAULT_MAX_SPEEDSIA_DRAGGING_TO_SHEET = 50.0;    // m/a; crit SIA speed for switch SSA --> SIA
const PetscScalar DEFAULT_MAXSLOPE_SSA = 1.0e-3; // no units/pure number; cap to avoid bad behavior
const PetscScalar DEFAULT_EPSILON_SSA = 1.0e15;  // kg m^-1 s^-1;  initial amount of (denominator) regularization in computation of effective viscosity
const PetscScalar DEFAULT_VERT_VEL_SSA = 0.0;  // temp evolution uses this value; incompressibility not satisfied
//const PetscScalar DEFAULT_BASAL_DRAG_COEFF_SSA = 2.0e9; // Pa s m^-1 Hulbe & MacAyeal (1999), p. 25,356
const PetscScalar DEFAULT_BASAL_DRAG_COEFF_SSA = 4.0e9; // seems to work better
const PetscScalar DEFAULT_TAUC = 1e4;  // 10^4 Pa = 0.1 bar
//used in iMvelocity.C and iMutil.C
const PetscScalar DEFAULT_MIN_TEMP_FOR_SLIDING = 273.0;  // note less than 
     // ice.meltingTemp; if above this value then decide to slide
const PetscScalar DEFAULT_INITIAL_AGE_YEARS = 1000.0;  // age to start age computation
const PetscScalar DEFAULT_OCEAN_HEAT_FLUX = 0.5;  // 0.5 W/m^2;
        // about 4 times more heating than peak of 
        // Shapiro&Ritzwoller geothermal fluxes (i.e. about 130 mW/m^2)
const PetscScalar DEFAULT_MAX_HMELT = 2.0;  // max of 2 m thick basal melt water layer

// see iMpdd.cc
const PetscScalar DEFAULT_PDD_STD_DEV = 5.0;  // K
const PetscScalar DEFAULT_PDD_FACTOR_SNOW = 0.003;  // (m ice-equivalent) day^-1 (deg C)^-1
const PetscScalar DEFAULT_PDD_FACTOR_ICE  = 0.008;  // (m ice-equivalent) day^-1 (deg C)^-1
const PetscScalar DEFAULT_PDD_REFREEZE_FRAC = 0.6;  // [pure fraction]
const PetscScalar DEFAULT_PDD_SUMMER_WARMING = 15.0;  //  K
     // re SUMMER_WARMING:  (30.38 - 0.006277 * 1000.0 - 0.3262 * 75.0)
     //                    - (49.13 - 0.007992 * 1000.0 -0.7576 * 75.0)
     //                   =  15.32   K
     // is result of EISMINT-GREENLAND formulas for h=1000.0 m and lat=75.0 deg N
const PetscScalar DEFAULT_PDD_SUMMER_PEAK_DAY = 243.0;  //  Julian day; August 1st

const PetscInt    DEFAULT_VERBOSITY_LEVEL = 2;

const PetscScalar DEFAULT_GRAIN_SIZE_INTERVAL_YEARS = 10.0;
const PetscScalar DEFAULT_MAX_TIME_STEP_YEARS = 60.0;  // years

const PetscScalar DEFAULT_ENHANCEMENT_FACTOR = 1.0;
const PetscTruth  DEFAULT_DO_MASS_CONSERVE = PETSC_TRUE;
const PetscTruth  DEFAULT_DO_TEMP = PETSC_TRUE;
const PetscTruth  DEFAULT_DO_TEMPSKIP = PETSC_FALSE; // off by default
const PetscInt    DEFAULT_TEMPSKIP_MAX = 10;
const PetscScalar DEFAULT_GLOBAL_MIN_ALLOWED_TEMP = 200.0;
const PetscInt    DEFAULT_MAX_LOW_TEMP_COUNT = 10;  // 

const PetscTruth  DEFAULT_INCLUDE_BMR_IN_CONTINUITY = PETSC_TRUE;
const PetscTruth  DEFAULT_REAL_AGE_FOR_GRAIN_SIZE = PETSC_FALSE;
const PetscTruth  DEFAULT_IS_DRY_SIMULATION = PETSC_FALSE;
const PetscTruth  DEFAULT_THERMAL_BEDROCK = PETSC_TRUE;
const PetscInt    DEFAULT_NOSPOKESLEVEL = 0;  // iterations of smoothing of Sigma
//const PetscScalar DEFAULT_MU_SLIDING = 3.17e-11;  // 100 m/a at 100kPa
const PetscScalar DEFAULT_MU_SLIDING = 0.0;

const PetscScalar DEFAULT_ISOTHERMAL_FLUX_N_EXPONENT = 3.0;

const PetscTruth  DEFAULT_DO_BED_DEF = PETSC_FALSE;
const PetscTruth  DEFAULT_DO_BED_ISO = PETSC_FALSE;
// model is so cheap you might as well update frequently:
const PetscScalar DEFAULT_BED_DEF_INTERVAL_YEARS = 10.0;  

const PetscTruth  DEFAULT_OCEAN_KILL = PETSC_FALSE;

const PetscTruth  DEFAULT_USE_SSA_VELOCITY = PETSC_FALSE;
const PetscTruth  DEFAULT_DO_PLASTIC_TILL = PETSC_FALSE;
const PetscTruth  DEFAULT_DO_SUPERPOSE = PETSC_FALSE;
const PetscTruth  DEFAULT_PURE_SUPERPOSE = PETSC_FALSE;
const PetscInt    DEFAULT_MAX_ITERATIONS_SSA = 300;
const PetscTruth  DEFAULT_USE_CONSTANT_NU_FOR_SSA = PETSC_FALSE;
const PetscTruth  DEFAULT_USE_CONSTANT_HARDNESS_FOR_SSA = PETSC_FALSE;
const PetscTruth  DEFAULT_COMPUTE_SURF_GRAD_INWARD_SSA = PETSC_FALSE;
// 30.0 * 1e6 * secpera = 9.45e14 is Ritz et al (2001) value of 30 MPa yr
const PetscScalar DEFAULT_CONSTANT_NU_FOR_SSA = 30.0 * 1.0e6 * secpera;
const PetscScalar DEFAULT_CONSTANT_HARDNESS_FOR_SSA = 1.9e8;  // Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
// for next constants, note (VELOCITY/LENGTH)^2  is very close to 10^-27; compare "\epsilon^2/L^2" which
// appears in formula (4.1) in C. Schoof 2006 "A variational approach to ice streams" J Fluid Mech 556 pp 227--251
const PetscScalar DEFAULT_PLASTIC_REGULARIZATION = 0.01 / secpera;
const PetscScalar DEFAULT_REGULARIZING_VELOCITY_SCHOOF = 1.0 / secpera;  // 1 m/a is small vel for stream/shelf
const PetscScalar DEFAULT_REGULARIZING_LENGTH_SCHOOF = 1000.0e3;         // 1000km is largish for dim of stream/shelf
const PetscScalar DEFAULT_SSA_RELATIVE_CONVERGENCE = 1.0e-4;

// pure number; pore water pressure is this fraction of overburden:
const PetscScalar DEFAULT_TILL_PW_FRACTION = 0.95;  
const PetscScalar DEFAULT_TILL_C_0 = 0.0;  // Pa; cohesion of till; 
            // note Schoof uses zero but Paterson pp 168--169 gives range 0 -- 40 kPa; but Paterson
            // notes that "... all the pairs c_0 and phi in the table would give a yield stress
            // for Ice Stream B that exceeds the basal shear stress there ..."
const PetscScalar DEFAULT_TILL_PHI = 30.0;  // pure number; tan(30^o) = ; till friction angle


//! Assigns default values to the many parameters and flags in IceModel.
/*!
The order of precedence for setting parameters in PISM is:
  - default values: Reasonable values to set up the model with are given in setDefaults()
    and in file pism/src/base/iMdefaults.  setDefaults() is called in the constructor for
    IceModel.  It would be reasonable to have setDefaults() read the defaults from a
    (default!) NetCDF file of a format so that others could be substituted.
  - derived class overrides:  The constructor of a derived class can choose its own 
    defaults for data members of IceModel (and its own data members).  These will override
    the above.
  - command line options:  The driver calls IceModel::setFromOptions() after the instance of
    IceModel (or a derived class) is constructed.  setFromOptions() is virtual but should 
    usually be called first if a derived class has a setFromOptions.

The input file (\c -if or \c -bif) will not contain (in Feb 2008 version of PISM) any values 
for the quantities which are set in setDefaults().  (There are parameters which can be set at
the command line or by the input file, like \c grid.Mx.  For \c -if the data file has the final
word but for -bif the command line options have the final word.)
 
The defaults should be reasonable values under all circumstances or they should indicate 
missing values in some manner.
 */
PetscErrorCode IceModel::setDefaults() {
  PetscErrorCode ierr;
  
  ierr = verbPrintf(3,grid.com, "setting IceModel defaults...\n"); CHKERRQ(ierr);

  initialized_p = PETSC_FALSE;

  // No X11 diagnostics by default, but allow them
  strcpy(diagnostic, "");
  strcpy(diagnosticBIG, "");
  showViewers = PETSC_TRUE;

  setVerbosityLevel(DEFAULT_VERBOSITY_LEVEL);
  ierr = setExecName("pism"); CHKERRQ(ierr);  // drivers typically override this

  enhancementFactor = DEFAULT_ENHANCEMENT_FACTOR;
  muSliding = DEFAULT_MU_SLIDING;
  thermalBedrock = DEFAULT_THERMAL_BEDROCK;
  doOceanKill = DEFAULT_OCEAN_KILL;
  
  ierr = grid.chooseEquallySpacedVertical(); CHKERRQ(ierr);
  
  computeSIAVelocities = PETSC_TRUE;
  transformForSurfaceGradient = PETSC_FALSE;

  useSSAVelocity = DEFAULT_USE_SSA_VELOCITY;
  ssaIntervalYears = -1.0;  // i.e. the default is to do an update every time step
  doPlasticTill = DEFAULT_DO_PLASTIC_TILL;
  doSuperpose = DEFAULT_DO_SUPERPOSE;
  pureSuperpose = DEFAULT_PURE_SUPERPOSE;
  ssaMaxIterations = DEFAULT_MAX_ITERATIONS_SSA;
  useConstantNuForSSA = DEFAULT_USE_CONSTANT_NU_FOR_SSA;
  useConstantHardnessForSSA = DEFAULT_USE_CONSTANT_HARDNESS_FOR_SSA;
  constantNuForSSA = DEFAULT_CONSTANT_NU_FOR_SSA;
  constantHardnessForSSA = DEFAULT_CONSTANT_HARDNESS_FOR_SSA;
  beta_default_drag_SSA = DEFAULT_BASAL_DRAG_COEFF_SSA;
  min_thickness_SSA = DEFAULT_MINH_SSA;
  regularizingVelocitySchoof = DEFAULT_REGULARIZING_VELOCITY_SCHOOF;
  regularizingLengthSchoof = DEFAULT_REGULARIZING_LENGTH_SCHOOF;
  ssaRelativeTolerance = DEFAULT_SSA_RELATIVE_CONVERGENCE;
  ssaEpsilon = DEFAULT_EPSILON_SSA;
  computeSurfGradInwardSSA = DEFAULT_COMPUTE_SURF_GRAD_INWARD_SSA;
  ssaSystemToASCIIMatlab = PETSC_FALSE;
  leaveNuAloneSSA = PETSC_FALSE;
  strcpy(ssaMatlabFilePrefix, "pism_SSA");

  plastic_till_pw_fraction = DEFAULT_TILL_PW_FRACTION;
  plastic_till_c_0 = DEFAULT_TILL_C_0;
  plastic_till_mu = tan((pi/180.0)*DEFAULT_TILL_PHI);
  plasticRegularization = DEFAULT_PLASTIC_REGULARIZATION;
  tauc_default_value = DEFAULT_TAUC;
  holdTillYieldStress = PETSC_FALSE;
  useConstantTillPhi = PETSC_TRUE;
  Hmelt_max = DEFAULT_MAX_HMELT;
  oceanHeatFlux = DEFAULT_OCEAN_HEAT_FLUX;

  pddFactorSnow = DEFAULT_PDD_FACTOR_SNOW;
  pddFactorIce = DEFAULT_PDD_FACTOR_ICE;
  pddRefreezeFrac = DEFAULT_PDD_REFREEZE_FRAC;
  pddSummerPeakDay = DEFAULT_PDD_SUMMER_PEAK_DAY;
  pddSummerWarming = DEFAULT_PDD_SUMMER_WARMING;
  pddStdDev = DEFAULT_PDD_STD_DEV;

  setMaxTimeStepYears(DEFAULT_MAX_TIME_STEP_YEARS);
  setAdaptTimeStepRatio(DEFAULT_ADAPT_TIMESTEP_RATIO);
  setAllGMaxVelocities(-1.0);

  run_year_default = DEFAULT_RUN_YEARS;
  setStartYear(DEFAULT_START_YEAR);
  setEndYear(DEFAULT_RUN_YEARS + DEFAULT_START_YEAR);
  yearsStartRunEndDetermined = PETSC_FALSE;
  initial_age_years_default = DEFAULT_INITIAL_AGE_YEARS;

  doMassConserve = DEFAULT_DO_MASS_CONSERVE;
  doTemp = DEFAULT_DO_TEMP;
  doTempSkip = DEFAULT_DO_TEMPSKIP;
  tempskipMax = DEFAULT_TEMPSKIP_MAX;
  reportHomolTemps = PETSC_TRUE;
  globalMinAllowedTemp = DEFAULT_GLOBAL_MIN_ALLOWED_TEMP;
  maxLowTempCount = DEFAULT_MAX_LOW_TEMP_COUNT;
  min_temperature_for_SIA_sliding = DEFAULT_MIN_TEMP_FOR_SLIDING;  
  includeBMRinContinuity = DEFAULT_INCLUDE_BMR_IN_CONTINUITY;
  isDrySimulation = DEFAULT_IS_DRY_SIMULATION;
  updateHmelt = PETSC_TRUE;

  flowLawUsesGrainSize = PETSC_FALSE;
  realAgeForGrainSize = DEFAULT_REAL_AGE_FOR_GRAIN_SIZE;

  doBedDef = DEFAULT_DO_BED_DEF;
  doBedIso = DEFAULT_DO_BED_ISO;
  bedDefIntervalYears = DEFAULT_BED_DEF_INTERVAL_YEARS;
  noSpokesLevel = DEFAULT_NOSPOKESLEVEL;
  doPDD = PETSC_FALSE;
  doPDDTrueRand = PETSC_FALSE;

  isothermalFlux_n_exponent = DEFAULT_ISOTHERMAL_FLUX_N_EXPONENT;
  
  // set default locations of soundings and slices
  id = (grid.Mx - 1)/2;
  jd = (grid.My - 1)/2;
  kd = 0;
  return 0;
}

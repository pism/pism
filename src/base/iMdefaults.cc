// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

// The order of precedence for setting parameters in PISM is:
//
//      Default values:
//              Reasonable values to set up the model with.
//              IceModel::setDefaults() below does this job and is virtual.  
//              Ideally all options will have defaults defined here.
//      Derived class overrides:
//              The derived class can override these options.
//      Command line options:
//              IceModel::setFromOptions().  This method is also virtual.  A derived
//              class may have its own setFromOptions(), at the beginning of which
//              IceModel::setFromOptions() should be called.
//      Data file (input) overrides:
//              The input (e.g. NetCDF) file may set quantities
//              The data always has highest precedence.  The user is
//              warned when their command line options are overridden.
//
// Thus the input data takes precedence over command line options which take 
// precedence over driver overrides which have precedence over the defaults below.
// The defaults should be reasonable values under all circumstances or they should
// indicate missing values in some manner.

const PetscScalar IceModel::DEFAULT_START_YEAR = 0;
const PetscScalar IceModel::DEFAULT_RUN_YEARS = 1000.0;  // years

//used in iMutil.C
const PetscScalar IceModel::DEFAULT_ADDED_TO_SLOPE_FOR_DIFF_IN_ADAPTIVE = 1.0e-4;
const PetscScalar IceModel::DEFAULT_ADDED_TO_GDMAX_ADAPT = 1.0e-2;
const PetscScalar IceModel::DEFAULT_ADAPT_TIMESTEP_RATIO = 0.12;  // yes, I'm confident this is o.k.

//used in iMIO.C
const PetscScalar IceModel::DEFAULT_h_VALUE_MISSING = 0.0;
const PetscScalar IceModel::DEFAULT_H_VALUE_MISSING = 0.0;
const PetscScalar IceModel::DEFAULT_BED_VALUE_MISSING = -5000.0;
const PetscScalar IceModel::DEFAULT_ACCUM_VALUE_MISSING = -0.5/ secpera;
const PetscScalar IceModel::DEFAULT_SURF_TEMP_VALUE_MISSING = 263.15;
const PetscScalar IceModel::DEFAULT_GEOTHERMAL_FLUX_VALUE_MISSING = 0.042; // J/m^2 s

//used in iMvelocity.C
const PetscScalar IceModel::DEFAULT_MINH_SSA = 5.0;  // m; minimum thickness for SSA velocity computation
const PetscScalar IceModel::DEFAULT_MIN_SHEET_TO_DRAGGING = 50.0;   // m/a; critical SIA speed for switch SIA --> SSA
const PetscScalar IceModel::DEFAULT_MAX_SPEED_DRAGGING_TO_SHEET = 5.0;  // m/a; crit Mac speed for switch SSA --> SIA
const PetscScalar IceModel::DEFAULT_MAX_SPEEDSIA_DRAGGING_TO_SHEET = 50.0;    // m/a; crit SIA speed for switch SSA --> SIA
const PetscScalar IceModel::DEFAULT_MAXSLOPE_SSA = 1.0e-3; // no units/pure number; cap to avoid bad behavior
const PetscScalar IceModel::DEFAULT_EPSILON_SSA = 1.0e15;  // kg m^-1 s^-1;  initial amount of (denominator) regularization in computation of effective viscosity
const PetscScalar IceModel::DEFAULT_EPSILON_MULTIPLIER_SSA = 4.0;  // no units/pure number; epsilon goes up by this ratio when
// previous value failed
const PetscScalar IceModel::DEFAULT_VERT_VEL_SSA = 0.0;  // temp evolution uses this value; incompressibility not satisfied
const PetscScalar IceModel::DEFAULT_MAX_VEL_FOR_CFL = 1000.0 / secpera;  // 10 km/a
//const PetscScalar IceModel::DEFAULT_BASAL_DRAG_COEFF_SSA = 2.0e9; // Pa s m^-1 Hulbe & MacAyeal (1999), p. 25,356
const PetscScalar IceModel::DEFAULT_BASAL_DRAG_COEFF_SSA = 4.0e9; // seems to work better
const PetscScalar IceModel::DEFAULT_TAUC = 1e4;  // 10^4 Pa = 0.1 bar
//used in iMvelocity.C and iMutil.C
const PetscScalar IceModel::DEFAULT_MIN_TEMP_FOR_SLIDING = 273.0;  // note less than 
     // ice.meltingTemp; if above this value then decide to slide
const PetscScalar IceModel::DEFAULT_INITIAL_AGE_YEARS = 1000.0;  // age to start age computation
const PetscScalar IceModel::DEFAULT_GRAIN_SIZE = 0.001;  // size of grains when assumed constant; for gk ice
const PetscScalar IceModel::DEFAULT_OCEAN_HEAT_FLUX = 0.5;  // 0.5 W/m^2;
        // about 4 times more heating than peak of 
        // Shapiro&Ritzwoller geothermal fluxes (i.e. about 130 mW/m^2)
const PetscScalar IceModel::DEFAULT_MAX_HMELT = 2.0;  // max of 2 m thick basal melt water layer

// see iMpdd.cc
const PetscScalar IceModel::DEFAULT_PDD_STD_DEV = 5.0;  // K
const PetscScalar IceModel::DEFAULT_PDD_FACTOR_SNOW = 0.003;  // (m ice-equivalent) day^-1 (deg C)^-1
const PetscScalar IceModel::DEFAULT_PDD_FACTOR_ICE  = 0.008;  // (m ice-equivalent) day^-1 (deg C)^-1
const PetscScalar IceModel::DEFAULT_PDD_REFREEZE_FRAC = 0.6;  // [pure fraction]
const PetscScalar IceModel::DEFAULT_PDD_SUMMER_WARMING = 15.0;  //  K
     // re SUMMER_WARMING:  (30.38 - 0.006277 * 1000.0 - 0.3262 * 75.0)
     //                    - (49.13 - 0.007992 * 1000.0 -0.7576 * 75.0)
     //                   =  15.32   K
     // is result of EISMINT-GREENLAND formulas for h=1000.0 m and lat=75.0 deg N
const PetscScalar IceModel::DEFAULT_PDD_SUMMER_PEAK_DAY = 243.0;  //  Julian day; August 1st

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
const PetscTruth  DEFAULT_DO_GRAIN_SIZE = PETSC_TRUE;
const PetscTruth  DEFAULT_REAL_AGE_FOR_GRAIN_SIZE = PETSC_FALSE;
const PetscTruth  DEFAULT_IS_DRY_SIMULATION = PETSC_FALSE;
const PetscTruth  DEFAULT_THERMAL_BEDROCK = PETSC_TRUE;
const PetscInt    DEFAULT_NOSPOKESLEVEL = 0;  // iterations of smoothing of Sigma
//const PetscScalar DEFAULT_MU_SLIDING = 3.17e-11;  // 100 m/a at 100kPa
const PetscScalar DEFAULT_MU_SLIDING = 0.0;

const PetscScalar DEFAULT_ISOTHERMAL_FLUX_N_EXPONENT = 3.0;

const PetscTruth  DEFAULT_DO_BED_DEF = PETSC_FALSE;
const PetscTruth  DEFAULT_DO_BED_ISO = PETSC_FALSE;
const PetscScalar DEFAULT_BED_DEF_INTERVAL_YEARS = 10.0;  // model is so cheap you might as well update frequently

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
PetscErrorCode IceModel::setDefaults() {

  //ierr = PetscPrintf(grid.com, "setting IceModel defaults...\n"); CHKERRQ(ierr);
  initialized_p = PETSC_FALSE;

  // No X11 diagnostics by default, but allow them
  strcpy(diagnostic, "");
  strcpy(diagnosticBIG, "");
  showViewers = PETSC_TRUE;

  setVerbosityLevel(DEFAULT_VERBOSITY_LEVEL);
  
  enhancementFactor = DEFAULT_ENHANCEMENT_FACTOR;
  muSliding = DEFAULT_MU_SLIDING;
  thermalBedrock = DEFAULT_THERMAL_BEDROCK;
  doOceanKill = DEFAULT_OCEAN_KILL;
  
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
  holdTillYieldStress = PETSC_FALSE;
  useConstantTillPhi = PETSC_TRUE;

  setMaxTimeStepYears(DEFAULT_MAX_TIME_STEP_YEARS);
  setAdaptTimeStepRatio(DEFAULT_ADAPT_TIMESTEP_RATIO);
  setAllGMaxVelocities(DEFAULT_MAX_VEL_FOR_CFL);

  setStartYear(DEFAULT_START_YEAR);
  setEndYear(DEFAULT_RUN_YEARS);
  yearsStartRunEndDetermined = PETSC_FALSE;

  doMassConserve = DEFAULT_DO_MASS_CONSERVE;
  doTemp = DEFAULT_DO_TEMP;
  doTempSkip = DEFAULT_DO_TEMPSKIP;
  tempskipMax = DEFAULT_TEMPSKIP_MAX;
  reportHomolTemps = PETSC_TRUE;
  globalMinAllowedTemp = DEFAULT_GLOBAL_MIN_ALLOWED_TEMP;
  maxLowTempCount = DEFAULT_MAX_LOW_TEMP_COUNT;
    
  includeBMRinContinuity = DEFAULT_INCLUDE_BMR_IN_CONTINUITY;
  doGrainSize = DEFAULT_DO_GRAIN_SIZE;
  realAgeForGrainSize = DEFAULT_REAL_AGE_FOR_GRAIN_SIZE;
  isDrySimulation = DEFAULT_IS_DRY_SIMULATION;
  updateHmelt = PETSC_TRUE;
  gsIntervalYears = DEFAULT_GRAIN_SIZE_INTERVAL_YEARS;
  doBedDef = DEFAULT_DO_BED_DEF;
  doBedIso = DEFAULT_DO_BED_ISO;
  bedDefIntervalYears = DEFAULT_BED_DEF_INTERVAL_YEARS;
  noSpokesLevel = DEFAULT_NOSPOKESLEVEL;
  doPDD = PETSC_FALSE;
  doPDDTrueRand = PETSC_FALSE;

  isothermalFlux_n_exponent = DEFAULT_ISOTHERMAL_FLUX_N_EXPONENT;
  
  // set default locations of soundings and slices
  id = (grid.p->Mx - 1)/2;
  jd = (grid.p->My - 1)/2;
  kd = 0;
  return 0;
}

// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstring>
#include "iceModel.hh"

// These are really default values.  They are set by IceModel::setDefaults()
// Ideally all options will have defaults defined here.

// The order of execution is:
//      Default values:
//              Reasonable values to set up the model with.
//              IceModel::setDefaults() does this job and is virtual
//      Driver overrides:
//              The driver program can override these options through
//              calls to IceModel::set[OPT]()
//      Command line options:
//              If the driver is willing to use command line options, it
//              should call IceModel::setFromOptions() _after_ any overrides
//              that it makes.  This method is also virtual.

// Thus command line options take precedence over driver overrides which
// have precedence over these defaults.  These should be reasonable value
// under all circumstances.

const PetscTruth  DEFAULT_BE_VERBOSE = PETSC_FALSE;
const PetscScalar DEFAULT_START_YEAR = 0;
const PetscScalar DEFAULT_RUN_YEARS = 1000.0;  // years
const PetscScalar DEFAULT_ENHANCEMENT_FACTOR = 1.0;
const PetscTruth  DEFAULT_DO_MASS_BAL = PETSC_TRUE;
const PetscTruth  DEFAULT_DO_TEMP = PETSC_TRUE;
const PetscTruth  DEFAULT_INCLUDE_BMR_IN_CONTINUITY = PETSC_TRUE;
const PetscTruth  DEFAULT_DO_GRAIN_SIZE = PETSC_TRUE;
const PetscTruth  DEFAULT_IS_DRY_SIMULATION = PETSC_FALSE;
const PetscScalar DEFAULT_GRAIN_SIZE_INTERVAL_YEARS = 60.0;
const PetscTruth  DEFAULT_DO_BED_DEF = PETSC_FALSE;
const PetscTruth  DEFAULT_DO_BED_ISO = PETSC_FALSE;
const PetscScalar DEFAULT_BED_DEF_INTERVAL_YEARS = 200.0;
const PetscScalar DEFAULT_MAX_TIME_STEP_YEARS = 60.0;  // years
const PetscInt    DEFAULT_TEMPSKIP = 1;  // every mass balance step is a temp step
const PetscInt    DEFAULT_NOSPOKESLEVEL = 0;  // iterations of smoothing of Sigma

const PetscScalar DEFAULT_ISOTHERMAL_FLUX_N_EXPONENT = 3;
const PetscScalar DEFAULT_ISOTHERMAL_FLUX_A_SOFTNESS = 1.0e-16 / secpera; // Pa^{-3} s^{-1}

const PetscTruth  DEFAULT_THERMAL_BEDROCK = PETSC_TRUE;
const PetscTruth  DEFAULT_OCEAN_KILL = PETSC_FALSE;
const PetscTruth  DEFAULT_USE_MACAYEAL_VELOCITY = PETSC_FALSE;
const PetscTruth  DEFAULT_USE_CONSTANT_NU_FOR_MACAYEAL = PETSC_FALSE;
const PetscScalar DEFAULT_CONSTANT_NU_FOR_MACAYEAL = 30.0 * 1.0e6 * secpera; // from Ritz et al (2001)
const PetscScalar DEFAULT_MACAYEAL_RELATIVE_CONVERGENCE = 1.0e-4;
const PetscScalar DEFAULT_MU_SLIDING = 3.17e-11;  // 100 m/a at 100kPa

PetscErrorCode IceModel::setDefaults() {
  PetscErrorCode ierr;
    
  //ierr = PetscPrintf(grid.com, "setting IceModel defaults...\n"); CHKERRQ(ierr);
  initialized_p = PETSC_FALSE;

  // No X11 diagnostics by default, but allow them
  strcpy(diagnostic, "");
  strcpy(diagnosticBIG, "");
  setShowViewers(PETSC_TRUE);

  setBeVerbose(DEFAULT_BE_VERBOSE);
  setEnhancementFactor(DEFAULT_ENHANCEMENT_FACTOR);
  setMuSliding(DEFAULT_MU_SLIDING);
  setThermalBedrock(DEFAULT_THERMAL_BEDROCK);
  setOceanKill(DEFAULT_OCEAN_KILL);
  setUseMacayealVelocity(DEFAULT_USE_MACAYEAL_VELOCITY);
  useConstantNuForMacAyeal = DEFAULT_USE_MACAYEAL_VELOCITY;
  constantNuForMacAyeal = DEFAULT_CONSTANT_NU_FOR_MACAYEAL;
  setMacayealRelativeTolerance(DEFAULT_MACAYEAL_RELATIVE_CONVERGENCE);
  setMacayealEpsilon(DEFAULT_EPSILON_MACAYEAL);
  
  setMaxTimeStepYears(DEFAULT_MAX_TIME_STEP_YEARS);
  setAdaptTimeStepRatio(DEFAULT_ADAPT_TIMESTEP_RATIO);
  setAllGMaxVelocities(DEFAULT_MAX_VEL_FOR_CFL);

  setStartYear(DEFAULT_START_YEAR);
  ierr = setRunYears(DEFAULT_RUN_YEARS); CHKERRQ(ierr);
  setDoMassBal(DEFAULT_DO_MASS_BAL);
  setDoTemp(DEFAULT_DO_TEMP);
  setIncludeBMRinContinuity(DEFAULT_INCLUDE_BMR_IN_CONTINUITY);
  setDoGrainSize(DEFAULT_DO_GRAIN_SIZE);
  setIsDrySimulation(DEFAULT_IS_DRY_SIMULATION);
  setGSIntervalYears(DEFAULT_GRAIN_SIZE_INTERVAL_YEARS);
  setDoBedDef(DEFAULT_DO_BED_DEF);
  setDoBedIso(DEFAULT_DO_BED_ISO);
  setBedDefIntervalYears(DEFAULT_BED_DEF_INTERVAL_YEARS);
  setTempskip(DEFAULT_TEMPSKIP);
  setAllowRegridding(PETSC_TRUE);
  setNoSpokes(DEFAULT_NOSPOKESLEVEL);
  
  setIsothermalFlux(PETSC_FALSE, DEFAULT_ISOTHERMAL_FLUX_N_EXPONENT,
                    DEFAULT_ISOTHERMAL_FLUX_A_SOFTNESS);
  return 0;
}

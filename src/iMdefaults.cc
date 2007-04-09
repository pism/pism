// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

const PetscInt    DEFAULT_VERBOSITY_LEVEL = 2;

const PetscScalar DEFAULT_START_YEAR = 0;
const PetscScalar DEFAULT_RUN_YEARS = 1000.0;  // years
const PetscScalar DEFAULT_GRAIN_SIZE_INTERVAL_YEARS = 60.0;
const PetscScalar DEFAULT_MAX_TIME_STEP_YEARS = 60.0;  // years

const PetscScalar DEFAULT_ENHANCEMENT_FACTOR = 1.0;
const PetscTruth  DEFAULT_DO_MASS_CONSERVE = PETSC_TRUE;
const PetscTruth  DEFAULT_DO_TEMP = PETSC_TRUE;
const PetscTruth  DEFAULT_DO_TEMPSKIP = PETSC_FALSE; // off by default
const PetscInt    DEFAULT_TEMPSKIP_MAX = 10;

const PetscTruth  DEFAULT_INCLUDE_BMR_IN_CONTINUITY = PETSC_TRUE;
const PetscTruth  DEFAULT_DO_GRAIN_SIZE = PETSC_TRUE;
const PetscTruth  DEFAULT_IS_DRY_SIMULATION = PETSC_FALSE;
const PetscTruth  DEFAULT_THERMAL_BEDROCK = PETSC_TRUE;
const PetscInt    DEFAULT_NOSPOKESLEVEL = 0;  // iterations of smoothing of Sigma
const PetscScalar DEFAULT_MU_SLIDING = 3.17e-11;  // 100 m/a at 100kPa

const PetscScalar DEFAULT_ISOTHERMAL_FLUX_N_EXPONENT = 3.0;
const PetscScalar DEFAULT_ISOTHERMAL_FLUX_A_SOFTNESS = 1.0e-16 / secpera; // Pa^{-3} s^{-1}

const PetscTruth  DEFAULT_DO_BED_DEF = PETSC_FALSE;
const PetscTruth  DEFAULT_DO_BED_ISO = PETSC_FALSE;
const PetscScalar DEFAULT_BED_DEF_INTERVAL_YEARS = 200.0;

const PetscTruth  DEFAULT_OCEAN_KILL = PETSC_FALSE;

const PetscTruth  DEFAULT_USE_MACAYEAL_VELOCITY = PETSC_FALSE;
const PetscTruth  DEFAULT_DO_SUPERPOSE = PETSC_FALSE;
const PetscInt    DEFAULT_MAX_ITERATIONS_MACAYEAL = 50;
const PetscTruth  DEFAULT_USE_CONSTANT_NU_FOR_MACAYEAL = PETSC_FALSE;
const PetscTruth  DEFAULT_USE_CONSTANT_HARDNESS_FOR_MACAYEAL = PETSC_FALSE;
const PetscTruth  DEFAULT_COMPUTE_SURF_GRAD_INWARD_MACAYEAL = PETSC_FALSE;
 // for next value, compare Ritz et al (2001)
const PetscScalar DEFAULT_CONSTANT_NU_FOR_MACAYEAL = 30.0 * 1.0e6 * secpera;
const PetscScalar DEFAULT_CONSTANT_HARDNESS_FOR_MACAYEAL = 1.9e8;  // Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
// for next constants, note (VELOCITY/LENGTH)^2  is very close to 10^-27; compare "\epsilon^2/L^2" which
// appears in formula (4.1) in C. Schoof 2006 "A variational approach to ice streams" J Fluid Mech 556 pp 227--251
const PetscScalar DEFAULT_REGULARIZING_VELOCITY_SCHOOF = 1.0 / secpera;  // 1 m/a is small vel for stream/shelf
const PetscScalar DEFAULT_REGULARIZING_LENGTH_SCHOOF = 1000.0e3;         // 1000km is largish for dim of stream/shelf
const PetscScalar DEFAULT_MACAYEAL_RELATIVE_CONVERGENCE = 1.0e-4;


PetscErrorCode IceModel::setDefaults() {
  PetscErrorCode ierr;

  //ierr = PetscPrintf(grid.com, "setting IceModel defaults...\n"); CHKERRQ(ierr);
  initialized_p = PETSC_FALSE;

  // No X11 diagnostics by default, but allow them
  strcpy(diagnostic, "");
  strcpy(diagnosticBIG, "");
  setShowViewers(PETSC_TRUE);

  setVerbosityLevel(DEFAULT_VERBOSITY_LEVEL);
  
  setEnhancementFactor(DEFAULT_ENHANCEMENT_FACTOR);
  setMuSliding(DEFAULT_MU_SLIDING);
  setThermalBedrock(DEFAULT_THERMAL_BEDROCK);
  setOceanKill(DEFAULT_OCEAN_KILL);
  
  setUseMacayealVelocity(DEFAULT_USE_MACAYEAL_VELOCITY);
  setDoSuperpose(DEFAULT_DO_SUPERPOSE);
  macayealMaxIterations = DEFAULT_MAX_ITERATIONS_MACAYEAL;
  useConstantNuForMacAyeal = DEFAULT_USE_CONSTANT_NU_FOR_MACAYEAL;
  useConstantHardnessForMacAyeal = DEFAULT_USE_CONSTANT_HARDNESS_FOR_MACAYEAL;
  constantNuForMacAyeal = DEFAULT_CONSTANT_NU_FOR_MACAYEAL;
  constantHardnessForMacAyeal = DEFAULT_CONSTANT_HARDNESS_FOR_MACAYEAL;
  setRegularizingVelocitySchoof(DEFAULT_REGULARIZING_VELOCITY_SCHOOF);
  setRegularizingLengthSchoof(DEFAULT_REGULARIZING_LENGTH_SCHOOF);
  setMacayealRelativeTolerance(DEFAULT_MACAYEAL_RELATIVE_CONVERGENCE);
  setMacayealEpsilon(DEFAULT_EPSILON_MACAYEAL);
  computeSurfGradInwardMacAyeal = DEFAULT_COMPUTE_SURF_GRAD_INWARD_MACAYEAL;

  setMaxTimeStepYears(DEFAULT_MAX_TIME_STEP_YEARS);
  setAdaptTimeStepRatio(DEFAULT_ADAPT_TIMESTEP_RATIO);
  setAllGMaxVelocities(DEFAULT_MAX_VEL_FOR_CFL);

  setStartYear(DEFAULT_START_YEAR);
  ierr = setRunYears(DEFAULT_RUN_YEARS); CHKERRQ(ierr);
  setDoMassConserve(DEFAULT_DO_MASS_CONSERVE);
  setDoTemp(DEFAULT_DO_TEMP);
  doTempSkip = DEFAULT_DO_TEMPSKIP;
  tempskipMax = DEFAULT_TEMPSKIP_MAX;
  
  setIncludeBMRinContinuity(DEFAULT_INCLUDE_BMR_IN_CONTINUITY);
  setDoGrainSize(DEFAULT_DO_GRAIN_SIZE);
  setIsDrySimulation(DEFAULT_IS_DRY_SIMULATION);
  setGSIntervalYears(DEFAULT_GRAIN_SIZE_INTERVAL_YEARS);
  setDoBedDef(DEFAULT_DO_BED_DEF);
  setDoBedIso(DEFAULT_DO_BED_ISO);
  setBedDefIntervalYears(DEFAULT_BED_DEF_INTERVAL_YEARS);
  setAllowRegridding(PETSC_TRUE);
  setNoSpokes(DEFAULT_NOSPOKESLEVEL);

  setIsothermalFlux(PETSC_FALSE, DEFAULT_ISOTHERMAL_FLUX_N_EXPONENT,
                    DEFAULT_ISOTHERMAL_FLUX_A_SOFTNESS);
  return 0;
}

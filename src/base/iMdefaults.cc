// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

The input file (\c -i or \c -boot_from) will not contain (in Feb 2008 version of PISM) any values 
for the quantities which are set in setDefaults().  (There are parameters which can be set at
the command line or by the input file, like \c grid.Mx.  For \c -i the data file has the final
word but for -boot_from the command line options have the final word.)
 
The defaults should be reasonable values under all circumstances or they should indicate 
missing values in some manner.
 */
PetscErrorCode IceModel::setDefaults() {
  PetscErrorCode ierr;
  
  ierr = verbPrintf(3,grid.com, "setting IceModel defaults...\n"); CHKERRQ(ierr);

  char alt_config[PETSC_MAX_PATH_LEN];
  PetscTruth use_alt_config;
  ierr = PetscOptionsGetString(PETSC_NULL, "-config", alt_config, PETSC_MAX_PATH_LEN, &use_alt_config);
  if (use_alt_config) {
    ierr = config.read(alt_config); CHKERRQ(ierr);
  } else {
    ierr = config.read(PISM_DEFAULT_CONFIG_FILE); CHKERRQ(ierr);
  }
  config.print();

  // No X11 diagnostics by default, but allow them
  strcpy(diagnostic, "");
  strcpy(diagnosticBIG, "");
  showViewers = PETSC_TRUE;

  ierr = setExecName("pism"); CHKERRQ(ierr);  // drivers typically override this

  enhancementFactor = config.get("enhancement_factor");
  muSliding = config.get("mu_sliding");
  thermalBedrock = config.get_flag("thermal_bedrock");
  doOceanKill = config.get_flag("ocean_kill");
  floatingIceKilled = config.get_flag("floating_ice_killed");

  grid.vertical_spacing = EQUAL;
  
  computeSIAVelocities = PETSC_TRUE;
  transformForSurfaceGradient = PETSC_FALSE;

  useSSAVelocity           = config.get_flag("use_ssa_velocity");
  doPlasticTill            = config.get_flag("do_plastic_till");
  doPseudoPlasticTill      = config.get_flag("do_pseudo_plastic_till");
  doSuperpose              = config.get_flag("do_superpose");
  ssaMaxIterations         = config.get("max_iterations_ssa");
  useConstantNuHForSSA     = config.get_flag("use_constant_nuh_for_ssa");
  ssaRelativeTolerance     = config.get("ssa_relative_convergence");
  ssaEpsilon               = config.get("epsilon_ssa");
  computeSurfGradInwardSSA = config.get_flag("compute_surf_grad_inward_ssa");
  ssaSystemToASCIIMatlab   = PETSC_FALSE;
  leaveNuHAloneSSA         = false;

  strcpy(ssaMatlabFilePrefix, "pism_SSA");

  plastic_till_pw_fraction  = config.get("till_pw_fraction");
  plastic_till_c_0          = config.get("till_c_0");
  plastic_till_mu           = tan((pi/180.0)*config.get("till_phi"));
  plasticRegularization     = config.get("plastic_regularization") / secpera;
  tauc_default_value        = config.get("tauc");
  pseudo_plastic_q          = config.get("pseudo_plastic_q");
  pseudo_plastic_uthreshold = config.get("pseudo_plastic_uthreshold") / secpera;
  holdTillYieldStress = PETSC_FALSE;
  useConstantTillPhi = PETSC_FALSE;
  
  shelvesDragToo = PETSC_FALSE;
  betaShelvesDragToo = config.get("beta_shelves_drag_too");
  
  Hmelt_max = config.get("max_hmelt");

  setMaxTimeStepYears(config.get("maximum_time_step_years"));
  setAdaptTimeStepRatio(config.get("adaptive_timestepping_ratio"));
  setAllGMaxVelocities(-1.0);

  run_year_default = config.get("run_length_years");
  setStartYear(config.get("start_year"));
  setEndYear(run_year_default + config.get("start_year"));
  yearsStartRunEndDetermined = PETSC_FALSE;
  initial_age_years_default       = config.get("initial_age_of_ice_years");

  doMassConserve                  = config.get_flag("do_mass_conserve");
  doTemp                          = config.get_flag("do_temp");
  doSkip                          = config.get_flag("do_skip");
  skipMax                         = config.get("skip_max");
  reportHomolTemps = PETSC_TRUE;
  globalMinAllowedTemp            = config.get("global_min_allowed_temp");
  maxLowTempCount                 = config.get("max_low_temp_count");
  min_temperature_for_SIA_sliding = config.get("minimum_temperature_for_sliding");
  includeBMRinContinuity          = config.get_flag("include_bmr_in_continuity");

  isDrySimulation = config.get_flag("is_dry_simulation");
  
  updateHmelt = PETSC_TRUE;

  realAgeForGrainSize = PETSC_FALSE;
  constantGrainSize   = config.get("constant_grain_size");

  doBedDef            = config.get_flag("do_bed_deformation");
  doBedIso            = config.get_flag("do_bed_iso");
  bedDefIntervalYears = config.get("bed_def_interval_years");
  noSpokesLevel       = config.get("no_spokes_level");

  // set default locations of soundings and slices
  id = (grid.Mx - 1)/2;
  jd = (grid.My - 1)/2;
  kd = 0;

  // default polar stereographic projection settings: South Pole
  psParams.svlfp = 0.0;
  psParams.lopo = -90.0;
  psParams.sp = -71.0;

  return 0;
}

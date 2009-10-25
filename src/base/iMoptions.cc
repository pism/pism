// Copyright (C) 2004--2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

//! Read runtime (command line) options and alter the corresponding parameters or flags as appropriate.
/*!

A critical principle of this procedure is that it will not alter IceModel parameters and flags
\e unless the user sets an option to do so.  Thus this base class setFromOptions() can be
called by a derived class after the derived class has set its own defaults for base class
parameters and flags.

Also, this procedure should not allocate memory or create new objects using the
new operator.

In fact this procedure only reads the majority of the options.  Some are read in 
initFromOptions(), writeFiles(), and setStartRunEndYearsFromOptions(), among other places.

Note there are no options to directly set \c dx, \c dy, \c dz, \c Lbz, and \c year as the user 
should not directly set these grid parameters.  There are, however, options for directly 
setting \c Mx, \c My, \c Mz, \c Mbz and also \c Lx, \c Ly, \c Lz.

Note that additional options are read by PISMClimateCoupler instances, including -pdd... 
and -d?forcing options.
 */
PetscErrorCode  IceModel::setFromOptions() {
  PetscErrorCode ierr;

  PetscTruth flag;

  PetscTruth  my_useConstantNuH, 
              myssaSystemToASCIIMatlab,
              myholdTillYieldStress, realageSet,
              etaSet, doShelvesDragToo;
  PetscScalar my_nuH;

  ierr = verbPrintf(3, grid.com,
		    "Processing physics-related command-line options...\n"); CHKERRQ(ierr);

  /* 
  note on pass-by-reference for options:
     For the last argument "flag" to PetscOptionsXXXX(....,&flag), the flag always indicates
     whether the option has been set.  Therefore "flag" is altered by this function call.
     For other arguments "value" to PetscOptionsXXXX(....,&value,&flag), the value of "value"
     is only set if the user specified the option.  Therefore "flag" should always be given a
     local PetscTruth variable if we want to preserve previously set IceModel flags.  By 
     contrast, for various parameters "value" we can use the IceModel parameter itself
     without fear of overwriting defaults unless, of course, the user wants them overwritten.
     It is also o.k. to have a local variable for "value", and then proceed to set the IceModel
     member accordingly.
  */

  ierr = config.scalar_from_option("adapt_ratio",
				   "adaptive_timestepping_ratio"); CHKERRQ(ierr);

  ierr = check_option("-bed_def_iso", flag); CHKERRQ(ierr);
  if (flag) {
    config.set_flag("do_bed_deformation", true);
    config.set_flag("do_bed_iso", true);
  }

  bool bed_def_iso = flag;

  ierr = check_option("-bed_def_lc", flag); CHKERRQ(ierr);  
  if (flag) {
    config.set_flag("do_bed_deformation", true);
    config.set_flag("do_bed_iso", false);
  }

  if (bed_def_iso && flag)  {
    ierr = verbPrintf(1,grid.com,
       "WARNING: both options -bed_def_iso and -bed_def_lc set; using Lingle & Clark model\n");
       CHKERRQ(ierr);
  }

  ierr = config.flag_from_option("bmr_in_cont", "include_bmr_in_continuity"); CHKERRQ(ierr);

// "-cbar_to_till" read in invertVelocitiesFromNetCDF() in iMinverse.cc

// "-chebZ" read in IceGrid::set_grid_from_options()

  ierr = PetscOptionsGetReal(PETSC_NULL, "-constant_nuH", &my_nuH, &my_useConstantNuH); CHKERRQ(ierr);
  // user gives nu*H in MPa yr m (e.g. Ritz et al 2001 value is 30.0 * 1.0)
  if (my_useConstantNuH == PETSC_TRUE) {
    setConstantNuHForSSA(my_nuH  * 1.0e6 * secpera); // convert to Pa s m
  }

// "-csurf_to_till" read in invertVelocitiesFromNetCDF() in iMinverse.cc

  /* FIXME: this is definitely needed, but needs work, too
  // regular size viewers
  ierr = PetscOptionsGetString(PETSC_NULL, "-d", diagnostic, PETSC_MAX_PATH_LEN, PETSC_NULL); 
            CHKERRQ(ierr);
  if (showViewers == PETSC_FALSE) {
    ierr = verbPrintf(1,grid.com,
       "WARNING: viewers requested with -d but showViewers is false, so none shown\n");
       CHKERRQ(ierr);
    strcpy(diagnostic, "\0");
  }
  
  // big viewers (which have higher priority than regular viewers)
  ierr = PetscOptionsGetString(PETSC_NULL, "-dbig", diagnosticBIG, PETSC_MAX_PATH_LEN, PETSC_NULL); 
            CHKERRQ(ierr);
  if (showViewers == PETSC_FALSE)  {
    ierr = verbPrintf(1,grid.com,
       "WARNING: viewers requested with -dbig but showViewers is false, so none shown\n");
       CHKERRQ(ierr);
    strcpy(diagnosticBIG, "\0");
  }
  */

  ierr = config.scalar_from_option("e", "enhancement_factor"); CHKERRQ(ierr);

  ierr = check_option("-eta", etaSet); CHKERRQ(ierr);
  if (etaSet == PETSC_TRUE)  transformForSurfaceGradient = PETSC_TRUE;

  ierr = config.flag_from_option("f3d", "force_full_diagnostics"); CHKERRQ(ierr);

  // whether or not to kill ice (zero thickness) if it is (or becomes) floating
  ierr = config.flag_from_option("float_kill", "floating_ice_killed"); CHKERRQ(ierr);

  // note "-gk" is used for specifying Goldsby-Kohlstedt ice
  //   this form allows a constant value of grain size to be input in mm
  ierr = config.scalar_from_option("gk", "constant_grain_size"); CHKERRQ(ierr);

  ierr = check_option("-gk", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = iceFactory.setType(ICE_HYBRID);CHKERRQ(ierr);
  }

  // note "-gk_age" is also used for specifying Goldsby-Kohlstedt ice;
  ierr = check_option("-gk_age", realAgeForGrainSize); CHKERRQ(ierr);
  if (realAgeForGrainSize) {
    ierr = iceFactory.setType(ICE_HYBRID);CHKERRQ(ierr);
  }

  ierr = check_option("-hold_tauc", myholdTillYieldStress); CHKERRQ(ierr);
  if (myholdTillYieldStress == PETSC_TRUE)    holdTillYieldStress = PETSC_TRUE;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-id", &id, PETSC_NULL); CHKERRQ(ierr);

// note "-i" is in use for input file name; see initFromOptions() in iMutil.cc

  ierr = PetscOptionsGetInt(PETSC_NULL, "-jd", &jd, PETSC_NULL); CHKERRQ(ierr);

  ierr = config.scalar_from_option("low_temp", "global_min_allowed_temp"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("max_dt",        "maximum_time_step_years"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("mu_sliding",    "mu_sliding");              CHKERRQ(ierr);
  ierr = config.scalar_from_option("max_low_temps", "max_low_temp_count");      CHKERRQ(ierr);

  ierr = config.flag_from_option("mass", "do_mass_conserve"); CHKERRQ(ierr);
  ierr = config.flag_from_option("temp", "do_temp"); CHKERRQ(ierr);

// note "-o" is in use for output file name; see iMIO.cc

  // whether or not to kill ice at locations where mask=FLOATING_OCEAN0;
  //   also determines if mask=FLOATING_OCEAN0 or mask=FLOATING
  //   at bootstrapping (-boot_from), if original condition was ice-free ocean
  ierr = config.flag_from_option("ocean_kill", "ocean_kill"); CHKERRQ(ierr);

  // use a plastic basal till mechanical model
  ierr = config.flag_from_option("plastic", "do_plastic_till"); CHKERRQ(ierr);

  // plastic_till_c_0 is a parameter in the computation of the till yield stress tau_c
  // from the thickness of the basal melt water; see updateYieldStressFromHmelt()
  // Note: option is given in kPa.
  ierr = config.scalar_from_option("plastic_c0", "till_c_0");      CHKERRQ(ierr);

  // till_pw_fraction is a parameter in the computation of the till yield stress tau_c
  // from the thickness of the basal melt water; see updateYieldStressFromHmelt()
  // option a pure number (a fraction); no conversion
  ierr = config.scalar_from_option("plastic_pwfrac", "till_pw_fraction"); CHKERRQ(ierr);

  // controls regularization of plastic basal sliding law
  ierr = config.scalar_from_option("plastic_reg", "plastic_regularization"); CHKERRQ(ierr);

  // "friction angle" in degrees
  ierr = config.scalar_from_option("plastic_phi", "default_till_phi"); CHKERRQ(ierr);

  // use pseudo plastic instead of pure plastic; see iMbasal.cc
  ierr = config.flag_from_option("pseudo_plastic", "do_pseudo_plastic_till"); CHKERRQ(ierr);

  // power in denominator on pseudo_plastic_uthreshold; typical is q=0.25; q=0 is pure plastic
  ierr = config.scalar_from_option("pseudo_plastic_q", "pseudo_plastic_q"); CHKERRQ(ierr);
  ierr = check_option("-pseudo_plastic_q", flag);  CHKERRQ(ierr);
  if (flag)
    config.set_flag("do_pseudo_plastic_till", true);

  // threshold; at this velocity tau_c is basal shear stress
  ierr = config.scalar_from_option("pseudo_plastic_uthreshold", "pseudo_plastic_uthreshold"); CHKERRQ(ierr);
  ierr = check_option("-pseudo_plastic_uthreshold", flag);  CHKERRQ(ierr);
  if (flag)
    config.set_flag("do_pseudo_plastic_till", true);

  // see updateGrainSizeNow(); option to choose modeled age vtau instead of pseudo age in
  // computing grainsize through Vostok core correlation
  ierr = check_option("-real_age_grainsize", realageSet); CHKERRQ(ierr);
  //if (realageSet == PETSC_TRUE)   realAgeForGrainSize = PETSC_TRUE;
  if (realageSet == PETSC_TRUE) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: -real_age_grainsize (sets realAgeForGrainSize) not implemented\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

// note "-regrid_from" is in use for regrid file name; see iMregrid.cc

// note "-regrid_vars" is in use for regrid variable names; see iMregrid.cc

  // see assembleSSAMatrix(); used in MISMIP
  ierr = check_option("-shelves_drag_too", doShelvesDragToo); CHKERRQ(ierr);
  if (doShelvesDragToo == PETSC_TRUE)   shelvesDragToo = PETSC_TRUE;
  
  ierr = config.flag_from_option("ssa", "use_ssa_velocity"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("ssa_eps",  "epsilon_ssa"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_maxi", "max_iterations_ssa"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_rtol", "ssa_relative_convergence"); CHKERRQ(ierr);
  
  // option to save linear system in Matlab-readable ASCII format at end of each
  // numerical solution of SSA equations; can be given with or without filename prefix
  // (i.e. "-ssa_matlab " or "-ssa_matlab foo" are both legal; in former case get 
  // "pism_SSA_[year].m" if "pism_SSA" is default prefix, and in latter case get "foo_[year].m")
  ierr = check_option("-ssa_matlab", myssaSystemToASCIIMatlab); CHKERRQ(ierr);
  if (myssaSystemToASCIIMatlab == PETSC_TRUE)   ssaSystemToASCIIMatlab = PETSC_TRUE;
  if (ssaSystemToASCIIMatlab == PETSC_TRUE) {  // now get the prefix if it was given by user
    char tempPrefix[PETSC_MAX_PATH_LEN];
    ierr = PetscOptionsGetString(PETSC_NULL, "-ssa_matlab", tempPrefix, 
             PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
    if (strlen(tempPrefix) > 0) {
      strcpy(ssaMatlabFilePrefix, tempPrefix);
    } // otherwise keep default prefix, whatever it was
  }

// -ssaBC used in IceROSSModel
  
  // apply "glaciological superposition to low order", i.e. add SIA results to those of 
  // SSA equations where DRAGGING; this version is  U = f(|v|) u + v   where u is SIA and v is SSA
  ierr = config.flag_from_option("super", "do_superpose"); CHKERRQ(ierr);
  
  /* This allows more than one mass continuity step per temperature/age and SSA
     computation */
  ierr = config.scalar_from_option("skip", "skip_max"); CHKERRQ(ierr);
  ierr = config.flag_from_option("skip",   "do_skip");  CHKERRQ(ierr);

  // Process -y, -ys, -ye. We are reading these options here because couplers
  // might need to know what year it is.
  ierr = set_time_from_options();

  return 0;
}

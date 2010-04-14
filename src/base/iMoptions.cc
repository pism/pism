// Copyright (C) 2004--2010 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <sstream>
#include <set>
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

Note that additional options are read by PISM{Atmosphere|Surface|Ocean}Model
instances, including -pdd... and -d?forcing options.
 */
PetscErrorCode  IceModel::setFromOptions() {
  PetscErrorCode ierr;

  bool flag;

  bool  my_useConstantNuH, 
    myssaSystemToASCIIMatlab,
    myholdTillYieldStress, realageSet,
    doShelvesDragToo;
  PetscReal my_nuH = 0;

  ierr = verbPrintf(3, grid.com,
		    "Processing physics-related command-line options...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Options overriding config flags and parameters", ""); CHKERRQ(ierr);

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

  ierr = config.flag_from_option("age", "do_age"); CHKERRQ(ierr);

  // FIXME: -bed_def options should be handled by the bed deformation module.
  // I'm leaving it as it is so far; will fix once that code is re-factored (CK).

  ierr = check_old_option_and_stop(grid.com, "-bed_def_iso", "-bed_def"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-bed_def_lc", "-bed_def"); CHKERRQ(ierr);

  // see getBasalWaterPressure()
  ierr = config.flag_from_option("bmr_enhance", "bmr_enhance_basal_water_pressure");
     CHKERRQ(ierr);
  // in units m a-1 :
  ierr = config.scalar_from_option("bmr_enhance_scale", "bmr_enhance_scale"); CHKERRQ(ierr);

  ierr = config.flag_from_option("bmr_in_cont", "include_bmr_in_continuity"); CHKERRQ(ierr);

// "-cbar_to_till" read in invertVelocitiesFromNetCDF() in iMinverse.cc

  // if set, use old IceModel::temperatureStep(), and set enthalpy as though
  //   ice is cold
  ierr = config.flag_from_option("cold", "do_cold_ice_methods"); CHKERRQ(ierr);

  ierr = PISMOptionsReal("-constant_nuH",
			 "Sets a constant value for the product of viscosity and thickness used in the SSA velocity computation",
			 my_nuH, my_useConstantNuH); CHKERRQ(ierr);
  // user gives nu*H in MPa yr m (e.g. Ritz et al 2001 value is 30.0 * 1.0)
  if (my_useConstantNuH) {
    setConstantNuHForSSA(my_nuH  * 1.0e6 * secpera); // convert to Pa s m
  }

// "-csurf_to_till" read in invertVelocitiesFromNetCDF() in iMinverse.cc

  ierr = config.scalar_from_option("e", "enhancement_factor"); CHKERRQ(ierr);

  ierr = config.flag_from_option("eta", "use_eta_transformation"); CHKERRQ(ierr);

  ierr = config.flag_from_option("f3d", "force_full_diagnostics"); CHKERRQ(ierr);

  // whether or not to kill ice (zero thickness) if it is (or becomes) floating
  ierr = config.flag_from_option("float_kill", "floating_ice_killed"); CHKERRQ(ierr);

  // note "-gk" is used for specifying Goldsby-Kohlstedt ice
  //   this form allows a constant value of grain size to be input in mm
  ierr = config.scalar_from_option("gk", "constant_grain_size"); CHKERRQ(ierr);

  ierr = PISMOptionsIsSet("-gk", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = iceFactory.setType(ICE_HYBRID);CHKERRQ(ierr);
  }

  // note "-gk_age" is also used for specifying Goldsby-Kohlstedt ice;
  ierr = PISMOptionsIsSet("-gk_age", flag); CHKERRQ(ierr);
  if (flag) {
    realAgeForGrainSize = PETSC_TRUE;
    ierr = iceFactory.setType(ICE_HYBRID);CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-hold_tauc", myholdTillYieldStress); CHKERRQ(ierr);
  if (myholdTillYieldStress == PETSC_TRUE)    holdTillYieldStress = PETSC_TRUE;

  ierr = PISMOptionsInt("-id", "Specifies the sounding row", id, flag); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-jd", "Specifies the sounding column", jd, flag); CHKERRQ(ierr);

  ierr = config.scalar_from_option("low_temp", "global_min_allowed_temp"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("max_low_temps", "max_low_temp_count"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("max_dt",        "maximum_time_step_years"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("mu_sliding",    "mu_sliding");              CHKERRQ(ierr);

  ierr = config.flag_from_option("mass", "do_mass_conserve"); CHKERRQ(ierr);

  ierr = config.flag_from_option("temp", "do_temp"); CHKERRQ(ierr);

// note "-o" is in use for output file name; see iMIO.cc

  // whether or not to kill ice at locations where mask=FLOATING_OCEAN0;
  //   also determines if mask=FLOATING_OCEAN0 or mask=FLOATING
  //   at bootstrapping (-boot_from), if original condition was ice-free ocean
  ierr = config.flag_from_option("ocean_kill", "ocean_kill"); CHKERRQ(ierr);


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

  ierr = PISMOptionsIsSet("-pseudo_plastic_q", flag);  CHKERRQ(ierr);
  if (flag)
    config.set_flag("do_pseudo_plastic_till", true);

  // threshold; at this velocity tau_c is basal shear stress
  ierr = config.scalar_from_option("pseudo_plastic_uthreshold", "pseudo_plastic_uthreshold"); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-pseudo_plastic_uthreshold", flag);  CHKERRQ(ierr);
  if (flag)
    config.set_flag("do_pseudo_plastic_till", true);

  // see updateGrainSizeNow(); option to choose modeled age vtau instead of pseudo age in
  // computing grainsize through Vostok core correlation
  ierr = PISMOptionsIsSet("-real_age_grainsize", realageSet); CHKERRQ(ierr);
  //if (realageSet == PETSC_TRUE)   realAgeForGrainSize = PETSC_TRUE;
  if (realageSet == PETSC_TRUE) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: -real_age_grainsize (sets realAgeForGrainSize) not implemented\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  // check -ssa_floating_only
  ierr = PISMOptionsIsSet("-ssa_floating_only", flag);  CHKERRQ(ierr);
  if (flag) {
    config.set_flag("use_ssa_velocity", true);
    config.set_flag("use_ssa_when_grounded", false);
  }

  // check -ssa_sliding
  ierr = PISMOptionsIsSet("-ssa_sliding", flag);  CHKERRQ(ierr);
  if (flag) {
    config.set_flag("use_ssa_velocity", true);
    config.set_flag("use_ssa_when_grounded", true);
  }

  // check if -super is set and warn if it has no effect:
  ierr = PISMOptionsIsSet("-super", flag); CHKERRQ(ierr);
  if (flag && (config.get_flag("use_ssa_when_grounded") == false)) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: option -super has no effect "
		      "if use_ssa_when_grounded is not set.\n"); CHKERRQ(ierr);
  }
  ierr = PISMOptionsIsSet("-no_super", flag); CHKERRQ(ierr);
  if (flag && (config.get_flag("use_ssa_when_grounded") == false)) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: option -no_super has no effect "
		      "if use_ssa_when_grounded is not set.\n"); CHKERRQ(ierr);
  }
  ierr = config.flag_from_option("super", "do_superpose"); CHKERRQ(ierr);

  ierr = check_old_option_and_stop(grid.com, "-ssa",
				   "-ssa_sliding' or '-ssa_floating_only"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-plastic", "-ssa_sliding"); CHKERRQ(ierr);

  // see assembleSSAMatrix(); used in MISMIP
  ierr = PISMOptionsIsSet("-shelves_drag_too", doShelvesDragToo); CHKERRQ(ierr);
  if (doShelvesDragToo == PETSC_TRUE)   shelvesDragToo = PETSC_TRUE;

  ierr = config.scalar_from_option("ssa_eps",  "epsilon_ssa"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_maxi", "max_iterations_ssa"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_rtol", "ssa_relative_convergence"); CHKERRQ(ierr);
  
  // option to save linear system in Matlab-readable ASCII format at end of each
  // numerical solution of SSA equations; can be given with or without filename prefix
  // (i.e. "-ssa_matlab " or "-ssa_matlab foo" are both legal; in former case get 
  // "pism_SSA_[year].m" if "pism_SSA" is default prefix, and in latter case get "foo_[year].m")

  string tempPrefix;
  ierr = PISMOptionsString("-ssa_matlab", "Save linear system in Matlab-readable ASCII format",
			   tempPrefix, myssaSystemToASCIIMatlab); CHKERRQ(ierr);
  if (myssaSystemToASCIIMatlab == PETSC_TRUE)   ssaSystemToASCIIMatlab = PETSC_TRUE;
  if (ssaSystemToASCIIMatlab == PETSC_TRUE) {  // now get the prefix if it was given by user
    if (tempPrefix.size() > 0) {
      strcpy(ssaMatlabFilePrefix, tempPrefix.c_str());
    } // otherwise keep default prefix, whatever it was
  }

// -ssaBC used in IceROSSModel
  
  
  /* This allows more than one mass continuity step per temperature/age and SSA
     computation */
  ierr = config.scalar_from_option("skip", "skip_max"); CHKERRQ(ierr);
  ierr = config.flag_from_option("skip",   "do_skip");  CHKERRQ(ierr);

  ierr = config.scalar_from_option("summary_volarea_scale_factor_log10",
                                   "summary_volarea_scale_factor_log10"); CHKERRQ(ierr);

  // if set, makes the thickness affect the pore_pressure; near margin there
  //   is a reduction in basal water pressure, a conceptual drainage mechanism
  ierr = config.flag_from_option("thk_eff", "thk_eff_basal_water_pressure"); CHKERRQ(ierr);
  // next two in  m  :
  ierr = config.scalar_from_option("thk_eff_H_high","thk_eff_H_high");  CHKERRQ(ierr);
  ierr = config.scalar_from_option("thk_eff_H_low","thk_eff_H_low");  CHKERRQ(ierr);
  // pure number :
  ierr = config.scalar_from_option("thk_eff_reduced","thk_eff_reduced");  CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  return 0;
}

//! Assembles a list of variables corresponding to an output file size.
PetscErrorCode IceModel::set_output_size(string option,
					 string description,
					 string default_value,
					 set<string> &result) {
  PetscErrorCode ierr;
  set<string> choices;
  string keyword;
  bool flag;

  choices.insert("small");
  choices.insert("medium");
  choices.insert("big");

  ierr = PISMOptionsList(grid.com, option,
			 description, choices,
			 default_value, keyword, flag); CHKERRQ(ierr);

  result.clear();

  // Add all the model-state variables:
  set<IceModelVec*> vars = variables.get_variables();
  set<IceModelVec*>::iterator i = vars.begin();
  while (i != vars.end()) {

    string intent = (*i)->string_attr("pism_intent");
    if ( (intent == "model_state") || (intent == "mapping") ||
	 (intent == "climate_steady") )
      result.insert((*i)->string_attr("name"));

    ++i;
  }

  // add {u,v}bar_ssa if SSA is "on":
  if (config.get_flag("use_ssa_velocity")) {
    result.insert("uvbar_ssa");	// will write both ubar_ssa and vbar_ssa
  }

  if (config.get_flag("do_age"))
    result.insert("age");

  if (config.get_flag("force_full_diagnostics"))
    keyword = "big";

  if (keyword == "small") {
    // only model-state variables are saved; we're done
    return 0;

  } else if (keyword == "medium") {
    // add all the variables listed in the config file ("medium" size):
    string tmp = config.get_string("output_medium");
    istringstream list(tmp);
  
    // split the list; note that this also removes any duplicate entries
    while (getline(list, tmp, ' ')) {
      if (!tmp.empty())		// this ignores multiple spaces separating variable names
	result.insert(tmp);
    }
    return 0;

  } else if (keyword == "big") {
    // add all the variables listed in the config file ("big" size):
    string tmp = config.get_string("output_big");
    istringstream list(tmp);
  
    // split the list; note that this also removes any duplicate entries
    while (getline(list, tmp, ' ')) {
      if (!tmp.empty())		// this ignores multiple spaces separating variable names
	result.insert(tmp);
    }
    return 0;

  } else {
    SETERRQ(1, "can't happen");
  }

  return 0;
}

// Copyright (C) 2004--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

//! \file iMoptions.cc Reading runtime options and setting configuration parameters.


//! Read runtime (command line) options and alter the corresponding parameters or flags as appropriate.
/*!
A critical principle of this procedure is that it will not alter IceModel parameters and flags
\e unless the user sets an option to do so.  This base class setFromOptions() can be
called by an IceModel-derived class after the it has set its own defaults.

In fact this procedure only reads the majority of the options.  Some are read in 
initFromOptions(), writeFiles(), and setStartRunEndYearsFromOptions(), among other places.

There are no options to directly set \c dx, \c dy, \c dz, \c Lbz, and \c year as the user 
should not directly set these grid parameters.  There are, however, options for directly 
setting \c Mx, \c My, \c Mz, \c Mbz and also \c Lx, \c Ly, \c Lz.

Note that additional options are read by PISM{Atmosphere|Surface|Ocean}Model
instances, including -pdd... and -d?forcing options.
 */
PetscErrorCode  IceModel::setFromOptions() {
  PetscErrorCode ierr;

  bool flag, myssaSystemToASCIIMatlab, realageSet;

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

  ierr = config.scalar_from_option("bed_smoother_range", "bed_smoother_range"); CHKERRQ(ierr);

  ierr = config.flag_from_option("blatter", "do_blatter"); CHKERRQ(ierr);

  // see getBasalWaterPressure()
  ierr = config.flag_from_option("bmr_enhance", "bmr_enhance_basal_water_pressure");
     CHKERRQ(ierr);
  // in units m a-1 :
  ierr = config.scalar_from_option("bmr_enhance_scale", "bmr_enhance_scale"); CHKERRQ(ierr);

  ierr = config.flag_from_option("bmr_in_cont", "include_bmr_in_continuity"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("calving_at_thickness", "calving_at_thickness"); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-calving_at_thickness", flag);  CHKERRQ(ierr);
  if (flag)  config.set_flag("do_thickness_calving", true);

  ierr = config.flag_from_option("cfbc", "calving_front_stress_boundary_condition"); CHKERRQ(ierr);

  // if set, use old IceModel::temperatureStep(), and set enthalpy as though
  //   ice is cold
  ierr = config.flag_from_option("cold", "do_cold_ice_methods"); CHKERRQ(ierr);

  ierr = config.flag_from_option("count_steps", "count_time_steps"); CHKERRQ(ierr);

  ierr = config.flag_from_option("diffuse_bwat", "do_diffuse_bwat"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("e", "enhancement_factor"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("e_ssa", "ssa_enhancement_factor"); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-e_ssa", flag);  CHKERRQ(ierr);
  if (flag)  config.set_flag("do_ssa_enhancement", true);

  ierr = config.scalar_from_option("eigen_calving", "eigen_calving"); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-eigen_calving", flag);  CHKERRQ(ierr);
  if (flag)  config.set_flag("do_eigen_calving", true);

  ierr = config.flag_from_option("dirichlet_bc", "dirichlet_bc"); CHKERRQ(ierr);


  bool gradient_set;
  string keyword;
  set<string> choices;
  choices.insert("eta");
  choices.insert("haseloff");
  choices.insert("mahaffy");
  ierr = PISMOptionsList(grid.com, "-gradient", "Surface gradient computation method.",
			 choices, "haseloff", keyword, gradient_set); CHKERRQ(ierr);
  if (gradient_set)  config.set_string("surface_gradient_method", keyword);

  // related old options
  ierr = check_old_option_and_stop(grid.com, "-eta", "-gradient"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-no_eta", "-gradient"); CHKERRQ(ierr);

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
    config.set_flag("compute_grain_size_using_age", true);
    ierr = iceFactory.setType(ICE_HYBRID);CHKERRQ(ierr);
  }

  ierr = PISMOptionsInt("-id", "Specifies the sounding row", id, flag); CHKERRQ(ierr);

  bool initfromT, initfromTandOm;
  ierr = PISMOptionsIsSet("-init_from_temp", initfromT); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-init_from_temp_and_liqfrac", initfromTandOm); CHKERRQ(ierr);

  ierr = config.string_from_option("institution", "institution"); CHKERRQ(ierr);

  ierr = PISMOptionsInt("-jd", "Specifies the sounding column", jd, flag); CHKERRQ(ierr);

  ierr = config.flag_from_option("kill_icebergs", "kill_icebergs"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("low_temp", "global_min_allowed_temp"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("max_low_temps", "max_low_temp_count"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("max_dt",        "maximum_time_step_years"); CHKERRQ(ierr);

  if (config.get("maximum_time_step_years") <= 0) {
    PetscPrintf(grid.com, "PISM ERROR: maximum_time_step_years has to be greater than 0.\n");
    PISMEnd();
  }

  ierr = config.scalar_from_option("mu_sliding",    "mu_sliding");              CHKERRQ(ierr);

  ierr = config.flag_from_option("mass", "do_mass_conserve"); CHKERRQ(ierr);

  // implements an option e.g. described in \ref Greve that is the
  // enhancement factor is coupled to the age of the ice with
  // e = 1 (A < 11'000 years), e = 3 otherwise
  ierr = PISMOptionsIsSet("-e_age_coupling", flag); CHKERRQ(ierr);
  if (flag) {
    config.set_flag("do_age", true);
    config.set_flag("do_e_age_coupling", true);
    ierr = verbPrintf(2, grid.com,
		      "  setting age-dependent enhancement factor: "
                      "e=1 if A<11'000 years, e=3 otherwise\n"); CHKERRQ(ierr);

  } else {
    config.set_flag("do_e_age_coupling", false);
  }

// note "-o" is in use for output file name; see iMIO.cc

  // whether or not to kill ice at locations where mask=FLOATING_OCEAN0;
  //   also determines if mask=FLOATING_OCEAN0 or mask=FLOATING
  //   at bootstrapping (-boot_file), if original condition was ice-free ocean
  ierr = config.flag_from_option("ocean_kill", "ocean_kill"); CHKERRQ(ierr);

  ierr = config.flag_from_option("part_grid", "part_grid"); CHKERRQ(ierr);

  ierr = config.flag_from_option("part_redist", "part_redist"); CHKERRQ(ierr);

  // option "-pik" turns on a suite of PISMPIK effects (but not -eigen_calving)
  ierr = PISMOptionsIsSet("-pik", "enable suite of PISM-PIK mechanisms", flag); CHKERRQ(ierr);
  if (flag) {
    config.set_flag("calving_front_stress_boundary_condition", true);
    config.set_flag("part_grid", true);
    config.set_flag("part_redist", true);
    config.set_flag("kill_icebergs", true);
  }

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
  if (flag)  config.set_flag("do_pseudo_plastic_till", true);

  // threshold; at this velocity tau_c is basal shear stress
  ierr = config.scalar_from_option("pseudo_plastic_uthreshold", "pseudo_plastic_uthreshold"); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-pseudo_plastic_uthreshold", flag);  CHKERRQ(ierr);
  if (flag)  config.set_flag("do_pseudo_plastic_till", true);

  // see updateGrainSizeNow(); option to choose modeled age vtau instead of pseudo age in
  // computing grainsize through Vostok core correlation
  ierr = PISMOptionsIsSet("-real_age_grainsize", realageSet); CHKERRQ(ierr);
  //if (realageSet == PETSC_TRUE)   realAgeForGrainSize = PETSC_TRUE;
  if (realageSet == PETSC_TRUE) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: -real_age_grainsize (sets realAgeForGrainSize) not implemented\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = config.flag_from_option("sia", "do_sia"); CHKERRQ(ierr);
  
  ierr = config.scalar_from_option("sliding_scale_brutal","sliding_scale_brutal"); CHKERRQ(ierr);  
  ierr = PISMOptionsIsSet("-sliding_scale_brutal", flag);  CHKERRQ(ierr);
  if (flag)  config.set_flag("scalebrutalSet", true);
  
  // check -ssa_floating_only
  ierr = PISMOptionsIsSet("-ssa_floating_only", flag);  CHKERRQ(ierr);
  if (flag) {
    config.set_flag("use_ssa_velocity", true);
    config.set_flag("use_ssa_when_grounded", false);
  }

  // Decide on the algorithm for solving the SSA
  set<string> ssa_choices;
  ssa_choices.insert("fem");
  ssa_choices.insert("fd");
  ssa_choices.insert("fd_pik");
  string ssa_method;
  bool ssa_method_set;
  ierr = PISMOptionsList(grid.com, "-ssa_method", "Algorithm for computing the SSA solution",
                         ssa_choices, ssa_method, ssa_method, ssa_method_set); CHKERRQ(ierr);
  if (ssa_method_set) {
    config.set_string("ssa_method",ssa_method);
  }

  // check -ssa_sliding
  ierr = PISMOptionsIsSet("-ssa_sliding", flag);  CHKERRQ(ierr);
  if (flag) {
    config.set_flag("use_ssa_velocity", true);
    config.set_flag("use_ssa_when_grounded", true);
  }

  ierr = check_old_option_and_stop(grid.com, "-ssa",
				   "-ssa_sliding' or '-ssa_floating_only"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-plastic", "-ssa_sliding"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("ssa_eps",  "epsilon_ssafd"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_maxi", "max_iterations_ssafd"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_rtol", "ssafd_relative_convergence"); CHKERRQ(ierr);
  
  // option to save linear system in Matlab-readable ASCII format at end of each
  // numerical solution of SSA equations; can be given with or without filename prefix
  // (i.e. "-ssa_matlab " or "-ssa_matlab foo" are both legal; in former case get 
  // "pism_SSA_[year].m" if "pism_SSA" is default prefix, and in latter case get "foo_[year].m")
  string tempPrefix;
  ierr = PISMOptionsString("-ssa_matlab", "Save linear system in Matlab-readable ASCII format",
			   tempPrefix, myssaSystemToASCIIMatlab); CHKERRQ(ierr);
  if (myssaSystemToASCIIMatlab == PETSC_TRUE)
    config.set_flag("write_ssa_system_to_matlab", true);
  
  /* This allows more than one mass continuity step per temperature/age and SSA
     computation */
  ierr = config.scalar_from_option("skip", "skip_max"); CHKERRQ(ierr);
  ierr = config.flag_from_option("skip",   "do_skip");  CHKERRQ(ierr);

  ierr = config.scalar_from_option("summary_volarea_scale_factor_log10",
                                   "summary_volarea_scale_factor_log10"); CHKERRQ(ierr);

  ierr = config.flag_from_option("energy", "do_energy"); CHKERRQ(ierr);

  // if set, makes the thickness affect the pore_pressure; near margin there
  //   is a reduction in basal water pressure, a conceptual drainage mechanism
  ierr = config.flag_from_option("thk_eff", "thk_eff_basal_water_pressure"); CHKERRQ(ierr);
  // next two in  m  :
  ierr = config.scalar_from_option("thk_eff_H_high","thk_eff_H_high");  CHKERRQ(ierr);
  ierr = config.scalar_from_option("thk_eff_H_low","thk_eff_H_low");  CHKERRQ(ierr);
  // pure number :
  ierr = config.scalar_from_option("thk_eff_reduced","thk_eff_reduced");  CHKERRQ(ierr);
  
  ierr = config.string_from_option("title", "run_title"); CHKERRQ(ierr);

  ierr = config.flag_from_option("vpik", "verbose_pik_messages");  CHKERRQ(ierr);
  if (getVerbosityLevel() > 2)  config.set_flag("verbose_pik_messages", true);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  global_attributes.set_string("title", config.get_string("run_title"));
  global_attributes.set_string("institution", config.get_string("institution"));
  global_attributes.set_string("command", pism_args_string());

  // warn about some option combinations
  
  if (!config.get_flag("do_mass_conserve") && config.get_flag("do_skip")) {
    ierr = verbPrintf(2, grid.com,
      "PISM WARNING: Both -skip and -no_mass are set.\n"
      "              -skip only makes sense in runs updating ice geometry.\n"); CHKERRQ(ierr);
  }

  if (config.get_flag("do_thickness_calving") && !config.get_flag("part_grid")) {
    ierr = verbPrintf(2, grid.com,
      "PISM WARNING: Calving at certain terminal ice thickness (-calving_at_thickness)\n"
      "              without application of partially filled grid cell scheme (-part_grid)\n"
      "              may lead to (incorrect) non-moving ice shelf front.\n"); CHKERRQ(ierr);
  }

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
  set<string> vars = variables.keys();

  set<string>::iterator i = vars.begin();
  while (i != vars.end()) {
    IceModelVec *var = variables.get(*i);

    string intent = var->string_attr("pism_intent");
    if ((intent == "model_state") || (intent == "mapping") || (intent == "climate_steady")) {
      result.insert(*i);
    }
    i++;
  }

  if (config.get_flag("do_age"))
    result.insert("age");

  if (config.get_flag("ocean_kill"))
    result.insert("ocean_kill_mask");

  if (config.get_flag("force_full_diagnostics"))
    keyword = "big";

  if (beddef != NULL)
    beddef->add_vars_to_output(keyword, result);

  if (btu != NULL)
    btu->add_vars_to_output(keyword, result);

  if (basal_yield_stress != NULL)
    basal_yield_stress->add_vars_to_output(keyword, result);

  // Ask the stress balance module to add more variables:
  if (stress_balance != NULL)
    stress_balance->add_vars_to_output(keyword, result);

  // Ask ocean and surface models to add more variables to the list:
  if (ocean != NULL)
    ocean->add_vars_to_output(keyword, result);

  if (surface != NULL)
    surface->add_vars_to_output(keyword, result);

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

    if (!config.get_flag("do_age"))
      result.erase("age");

    return 0;
  } else {
    SETERRQ(1, "can't happen");
  }

  return 0;
}


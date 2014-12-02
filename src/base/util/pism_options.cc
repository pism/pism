// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include <sstream>
#include <algorithm>
#include <string.h>
#include <assert.h>

#include "pism_options.hh"
#include "NCVariable.hh"
#include "PISMConfig.hh"

#include "error_handling.hh"

namespace pism {

//! \brief Stop if -version is set.
PetscErrorCode stop_on_version_option() {
  PetscErrorCode ierr;

  bool vSet = false;
  ierr = OptionsIsSet("-version", vSet); CHKERRQ(ierr);
  if (vSet == false) {
    return 0;
  }

  // FIXME!
  // PISMEndQuiet();
  return 0;
}

//! Determine verbosity level from user options.
/*!
\verbatim
   level  option        meaning
   -----  ------        -------
   0      -verbose 0    never print to std out AT ALL!
   1      -verbose 1    less verbose than default: thresh must be 1 to print
   2      -verbose 2    DEFAULT
   3      -verbose 3    somewhat verbose
          -verbose      same as "-verbose 3"
   4      -verbose 4    fairly verbose
   5      -verbose 5    very verbose: print everything
\endverbatim
See verbPrintf().
 */
PetscErrorCode verbosityLevelFromOptions() {
  PetscErrorCode ierr;
  PetscInt       myLevel;
  PetscBool     verbose, levelSet;

  ierr = setVerbosityLevel(2);
  ierr = PetscOptionsGetInt(NULL, "-verbose", &myLevel, &levelSet);
  PISM_PETSC_CHK(ierr, "PetscOptionsGetInt");
  if (levelSet == PETSC_TRUE) {
    ierr = setVerbosityLevel(myLevel);
  } else {
    ierr = PetscOptionsHasName(NULL, "-verbose", &verbose);
    PISM_PETSC_CHK(ierr, "PetscOptionsHasName");
    if (verbose == PETSC_TRUE) {
      ierr = setVerbosityLevel(3);
    }
  }
  return 0;
}

//! Print a warning telling the user that an option was ignored.
PetscErrorCode ignore_option(MPI_Comm com, std::string name) {
  PetscErrorCode ierr;
  PetscBool option_is_set;

  char tmp[1]; // dummy string
  ierr = PetscOptionsGetString(NULL, name.c_str(), tmp, 1, &option_is_set);
  PISM_PETSC_CHK(ierr, "PetscOptionsGetString");

  if (option_is_set) {
    verbPrintf(1, com, "PISM WARNING: ignoring command-line option '%s'.\n",
               name.c_str());
  }

  return 0;
}

//! Stop if an option `old_name` is set, printing a message that `new_name` should be used instead.
PetscErrorCode check_old_option_and_stop(std::string old_name, std::string new_name) {
  PetscErrorCode ierr;
  PetscBool option_is_set;

  char tmp[1]; // dummy string
  ierr = PetscOptionsGetString(NULL, old_name.c_str(), tmp, 1, &option_is_set);
  PISM_PETSC_CHK(ierr, "PetscOptionsGetString");

  if (option_is_set) {
    throw RuntimeError::formatted("command-line option '%s' is deprecated. Please use '%s' instead.",
                                  old_name.c_str(), new_name.c_str());
  }

  return 0;
}

//!Stop if an option `name` is set.
PetscErrorCode stop_if_set(std::string name) {
  PetscErrorCode ierr;
  PetscBool option_is_set;

  char tmp[1]; // dummy string
  ierr = PetscOptionsGetString(NULL, name.c_str(), tmp, 1, &option_is_set);
  PISM_PETSC_CHK(ierr, "PetscOptionsGetString");

  if (option_is_set) {
    throw RuntimeError::formatted("command-line option '%s' is not allowed.",
                                  name.c_str());
  }

  return 0;
}

//! \brief Print a usage message.
PetscErrorCode just_show_usage(MPI_Comm com, std::string execname, std::string usage) {
  verbPrintf(1,com,
             "%s is a PISM (http://www.pism-docs.org) executable.\nOptions cheat-sheet:\n",
             execname.c_str());
  verbPrintf(1,com,usage.c_str());
  verbPrintf(1,com,
             "Parallel run using N processes (typical case):  mpiexec -n N %s ...\n"
             "For more help with PISM:\n"
             "  1. download PDF User's Manual:\n"
             "       http://www.pism-docs.org/wiki/lib/exe/fetch.php?media=manual.pdf\n"
             "  2. read browser for technical details:\n"
             "       http://www.pism-docs.org/doxy/html/index.html\n"
             "  3. view issues/bugs at source host: https://github.com/pism/pism/issues\n"
             "  4. do '%s -help | grep foo' to see PISM and PETSc options with 'foo'.\n"
             "  5. email for help:  help@pism-docs.org\n",
             execname.c_str(), execname.c_str());
  return 0;
}


//! @brief Show provided usage message and quit. (Consider using
//! show_usage_check_req_opts() in preference to this one.)
PetscErrorCode show_usage_and_quit(MPI_Comm com, std::string execname, std::string usage) {

  stop_on_version_option();

  just_show_usage(com, execname, usage);

  // FIXME!
  // PISMEndQuiet();
  return 0;
}


//! @brief In a single call a driver program can provide a usage string to
//! the user and check if required options are given, and if not, end.
void show_usage_check_req_opts(MPI_Comm com, std::string execname,
                               std::vector<std::string> required_options,
                               std::string usage) {
  stop_on_version_option();

  bool usageSet = false;
  OptionsIsSet("-usage", usageSet);
  if (usageSet == true) {
    show_usage_and_quit(com, execname, usage);
  }

  // go through list of required options, and if not given, fail
  bool req_absent = false;
  for (size_t ii=0; ii < required_options.size(); ii++) {
    bool set = false;
    OptionsIsSet(required_options[ii], set);
    if (set == false) {
      req_absent = true;
      verbPrintf(1,com,
                 "PISM ERROR: option %s required\n",required_options[ii].c_str());
    }
  }
  if (req_absent == true) {
    verbPrintf(1,com,"\n");
    show_usage_and_quit(com, execname, usage);
  }

  // show usage message with -help, but don't fail
  bool helpSet = false;
  OptionsIsSet("-help", helpSet);
  if (helpSet == true) {
    just_show_usage(com, execname, usage);
  }
}


/*
   note on pass-by-reference for options: For the last argument "flag" to
   PetscOptionsXXXX(....,&flag), the flag always indicates whether the option
   has been set. Therefore "flag" is altered by this function call. For other
   arguments "value" to PetscOptionsXXXX(....,&value,&flag), the value of
   "value" is only set if the user specified the option. Therefore "flag"
   should always be given a local PetscBool variable if we want to preserve
   previously set IceModel flags. By contrast, for various parameters "value"
   we can use the IceModel parameter itself without fear of overwriting
   defaults unless, of course, the user wants them overwritten. It is also o.k.
   to have a local variable for "value", and then proceed to set the IceModel
   member accordingly.
*/

//! PISM wrapper replacing PetscOptionsEList.
/*
  Ignores everything after the first colon, i.e. if "-foo bar" is allowed, then
  "-foo bar:baz" is also allowed.

  This is to make it possible to pass a parameter to a module selected using a
  command-line option without adding one mode option.
 */
PetscErrorCode OptionsList(std::string opt, std::string description,
                           std::set<std::string> choices, std::string default_value,
                           std::string &result, bool &flag) {
  PetscErrorCode ierr;
  char tmp[TEMPORARY_STRING_LENGTH];
  std::string list, descr;
  PetscBool opt_set = PETSC_FALSE;

  if (choices.empty()) {
    throw RuntimeError::formatted("OptionsList: empty choices argument");
  }

  std::set<std::string>::iterator j = choices.begin();
  list = "[" + *j++;
  while (j != choices.end()) {
    list += ", " + (*j++);
  }
  list += "]";

  descr = description + " Choose one of " + list;

  ierr = PetscOptionsString(opt.c_str(), descr.c_str(), "", default_value.c_str(),
                            tmp, TEMPORARY_STRING_LENGTH, &opt_set);
  PISM_PETSC_CHK(ierr, "PetscOptionsString");

  // return the default value if the option was not set
  if (!opt_set) {
    flag = false;
    result = default_value;
    return 0;
  }

  std::string keyword = tmp;
  // find ":" and discard everything that goes after
  size_t n = keyword.find(":");
  if (n != std::string::npos) {
    keyword.resize(n);
  }

  // return the choice if it is valid and stop if it is not
  if (choices.find(keyword) != choices.end()) {
    flag = true;
    result = tmp;
  } else {
    throw RuntimeError::formatted("invalid %s argument: '%s'. Please choose one of %s.\n",
                                  opt.c_str(), tmp, list.c_str());
  }

  return 0;
}

//! \brief Process a command-line option taking a string as an argument.
PetscErrorCode OptionsString(std::string option, std::string text,
                             std::string &result, bool &is_set, bool allow_empty_arg) {
  PetscErrorCode ierr;
  char tmp[TEMPORARY_STRING_LENGTH];
  PetscBool flag;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "",
                            result.c_str(), tmp,
                            TEMPORARY_STRING_LENGTH, &flag);
  PISM_PETSC_CHK(ierr, "PetscOptionsString");

  is_set = (flag == PETSC_TRUE);

  if (is_set) {
    if (strlen(tmp) == 0) {
      if (allow_empty_arg) {
        result.clear();
      } else {
        throw RuntimeError::formatted("command line option '%s' requires an argument.\n",
                                      option.c_str());
      }
    } else {
      result = tmp;
    }
  }

  return 0;
}

//! PISM wrapper replacing PetscOptionsStringArray.
PetscErrorCode OptionsStringArray(std::string opt, std::string text, std::string default_value,
				      std::vector<std::string>& result, bool &flag) {
  PetscErrorCode ierr;
  char tmp[TEMPORARY_STRING_LENGTH];
  PetscBool opt_set = PETSC_FALSE;

  ierr = PetscOptionsString(opt.c_str(), text.c_str(), "", default_value.c_str(),
                            tmp, TEMPORARY_STRING_LENGTH, &opt_set);
  PISM_PETSC_CHK(ierr, "PetscOptionsString");

  result.clear();

  std::string word;
  if (opt_set) {
    std::istringstream arg(tmp);
    while (getline(arg, word, ',')) {
      result.push_back(word);
    }

    if (result.empty()) {
      throw RuntimeError::formatted("command line option '%s' requires an argument.\n",
                                    opt.c_str());
    }

    flag = true;
  } else {                      // parse the default list given
    std::istringstream arg(default_value);
    while (getline(arg, word, ',')) {
      result.push_back(word);
    }

    flag = false;
  }

  return 0;
}

//! Process a command-line option and return a set of strings.
PetscErrorCode OptionsStringSet(std::string opt, std::string text, std::string default_value,
				    std::set<std::string>& result, bool &flag) {
  std::vector<std::string> tmp;

  OptionsStringArray(opt, text, default_value, tmp, flag);

  result.clear();
  std::vector<std::string>::iterator j = tmp.begin();
  while(j != tmp.end()) {
    result.insert(*j);
    ++j;
  }

  return 0;
}

//! \brief Process a command-line option taking an integer as an argument.
PetscErrorCode OptionsInt(std::string option, std::string text,
                          int &result, bool &is_set) {
  PetscErrorCode ierr;
  char str[TEMPORARY_STRING_LENGTH];
  PetscBool flag;
  char *endptr;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "", "none", str,
                            TEMPORARY_STRING_LENGTH, &flag);
  PISM_PETSC_CHK(ierr, "PetscOptionsString");

  is_set = (flag == PETSC_TRUE);

  if (is_set == false) {
    return 0;
  }

  if (strlen(str) == 0) {
    throw RuntimeError::formatted("command line option '%s' requires an argument.",
                                  option.c_str());
  }

  result = (int) strtol(str, &endptr, 10);
  if (*endptr != '\0') {
    throw RuntimeError::formatted("Can't parse '%s %s': (%s is not an integer).\n",
                                  option.c_str(), str, str);
  }

  return 0;
}

//! \brief Process a command-line option taking a real number as an argument.
PetscErrorCode OptionsReal(std::string option, std::string text,
			       double &result, bool &is_set) {
  PetscErrorCode ierr;
  char str[TEMPORARY_STRING_LENGTH];
  PetscBool flag;
  char *endptr;

  memset(str, 0, TEMPORARY_STRING_LENGTH);
  snprintf(str, TEMPORARY_STRING_LENGTH, "%f", result);

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "", str, str,
                            TEMPORARY_STRING_LENGTH, &flag);
  PISM_PETSC_CHK(ierr, "PetscOptionsString");

  is_set = (flag == PETSC_TRUE);

  if (is_set == false) {
    return 0;
  }

  if (strlen(str) == 0) {
    throw RuntimeError::formatted("command line option '%s' requires an argument.",
                                  option.c_str());
  }

  result = strtod(str, &endptr);
  if (*endptr != '\0') {
    throw RuntimeError::formatted("Can't parse '%s %s': (%s is not a number).",
                                  option.c_str(), str, str);
  }

  return 0;
}
//! \brief Process a command-line option taking a comma-separated list of reals
//! as an argument.
PetscErrorCode OptionsRealArray(std::string option, std::string text,
				    std::vector<double> &result, bool &is_set) {
  PetscErrorCode ierr;
  char str[TEMPORARY_STRING_LENGTH];
  PetscBool flag;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "",
                            "none", str,
                            TEMPORARY_STRING_LENGTH, &flag);
  PISM_PETSC_CHK(ierr, "PetscOptionsString");

  is_set = (flag == PETSC_TRUE);

  if (is_set) {
    std::istringstream arg(str);
    std::string tmp;

    result.clear();
    while(getline(arg, tmp, ',')) {
      double d;
      char *endptr;

      d = strtod(tmp.c_str(), &endptr);
      if (*endptr != '\0') {
        throw RuntimeError::formatted("Can't parse %s (%s is not a number).",
                                      tmp.c_str(), tmp.c_str());
      }
      else {
        result.push_back(d);
      }
    }
  }

  return 0;
}

//! \brief Process a command-line option taking a comma-separated list of
//! integers as an argument.
PetscErrorCode OptionsIntArray(std::string option, std::string text,
				   std::vector<int> &result, bool &is_set) {
  std::vector<double> tmp;

  OptionsRealArray(option, text, tmp, is_set);

  result.clear();
  for (unsigned int j = 0; j < tmp.size(); ++j) {
    result.push_back(static_cast<int>(tmp[j]));
  }

  return 0;
}

//! Checks if an option is present in the PETSc option database.
/*!

  This is (essentially) a reimplementation of PetscOptionsHasName, except that
  this *always* sets `flag` to PETSC_TRUE if an option is present.

  PetscOptionsHasName, on the other hand, sets `flag` to PETSC_FALSE if an
  option was set as "-foo FALSE", "-foo NO" or "-foo 0". Note that if one uses
  "-foo 0.0", PetscOptionsHasName will set `flag` to PETSC_TRUE.

  This unpredictability is bad. We want a function that does not depend on the
  argument given with an option.
 */
PetscErrorCode OptionsIsSet(std::string option, bool &result) {
  PetscErrorCode ierr;
  char tmp[1];
  PetscBool flag;

  ierr = PetscOptionsGetString(NULL, option.c_str(), tmp, 1, &flag);
  PISM_PETSC_CHK(ierr, "PetscOptionsGetString");

  result = (flag == PETSC_TRUE);

  return 0;
}

//! A version of OptionsIsSet that prints a -help message.
PetscErrorCode OptionsIsSet(std::string option, std::string text,
				bool &result) {
  PetscErrorCode ierr;
  char tmp[1];
  PetscBool flag;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "",
                            "", tmp, 1, &flag);
  PISM_PETSC_CHK(ierr, "PetscOptionsString");

  result = (flag == PETSC_TRUE);

  return 0;
}


/**
 * Checks if a command-line option `option` is set and is given an argument.
 *
 * @param option name of the option
 * @param result true if `option` is set and has an argument, `false` otherwise
 *
 * @return 0 on success
 */
PetscErrorCode OptionsHasArgument(std::string option, bool &result) {
  std::string tmp;

  OptionsString(option, "", tmp, result, true);

  if (result == false || tmp.empty() == true) {
    result = false;
  }

  return 0;
}

//! Initializes the config parameter and flag database.
/*!
  Processes -config and -config_override command line options.
 */
PetscErrorCode init_config(MPI_Comm com,
			   Config &config, Config &overrides,
			   bool process_options) {
  PetscErrorCode ierr;

  std::string alt_config = PISM_DefaultConfigFile,
    override_config;
  bool use_alt_config, use_override_config;

  ierr = PetscOptionsBegin(com, "", "PISM config file options", "");
  PISM_PETSC_CHK(ierr, "PetscOptionsBegin");
  {
    OptionsString("-config", "Specifies the name of an alternative config file",
                  alt_config, use_alt_config);
    OptionsString("-config_override", "Specifies a config override file name",
                  override_config, use_override_config);
  }
  ierr = PetscOptionsEnd();
  PISM_PETSC_CHK(ierr, "PetscOptionsEnd");

  config.read(alt_config);

  if (use_override_config) {
    overrides.read(override_config);
    config.import_from(overrides);
    verbPrintf(2, com, "CONFIG OVERRIDES read from file '%s'.\n",
               override_config.c_str());
  }

  if (process_options) {
    set_config_from_options(config);
  }

  config.print_to_stdout();

  return 0;
}

PetscErrorCode set_config_from_options(Config &config) {
  bool flag;

  config.keyword_from_option("periodicity", "grid_periodicity", "none,x,y,xy");

  // Energy modeling
  config.flag_from_option("varc", "use_linear_in_temperature_heat_capacity");
  config.flag_from_option("vark",
                          "use_temperature_dependent_thermal_conductivity");

  config.flag_from_option("bmr_in_cont", "include_bmr_in_continuity");

  {
    bool energy_set;
    std::string energy;
    std::set<std::string> choices;
    choices.insert("none");
    choices.insert("cold");
    choices.insert("enthalpy");

    OptionsList("-energy",
                "choose the energy model (one of 'none', 'cold', 'enthalpy')",
                choices, "enthalpy", energy, energy_set);

    if (energy_set == true) {
      if (energy == "none") {
        config.set_flag_from_option("do_energy", false);
        // Allow selecting cold ice flow laws in isothermal mode. 
        config.set_flag_from_option("do_cold_ice_methods", true);
      } else if (energy == "cold") {
        config.set_flag_from_option("do_energy", true);
        config.set_flag_from_option("do_cold_ice_methods", true);
      } else if (energy == "enthalpy") {
        config.set_flag_from_option("do_energy", true);
        config.set_flag_from_option("do_cold_ice_methods", false);
      } else {
        // can't happen (OptionsList validates its input)
        assert(false);
      }
    }
  }

  // at bootstrapping, choose whether the method uses smb as upper boundary for
  // vertical velocity
  config.keyword_from_option("boot_temperature_heuristic",
                             "bootstrapping_temperature_heuristic", "smb,quartic_guess");

  config.scalar_from_option("low_temp", "global_min_allowed_temp");
  config.scalar_from_option("max_low_temps", "max_low_temp_count");

  // Sub-models
  config.flag_from_option("age", "do_age");
  config.flag_from_option("mass", "do_mass_conserve");

  // hydrology
  config.keyword_from_option("hydrology", "hydrology_model",
                             "null,routing,distributed");
  config.flag_from_option("hydrology_use_const_bmelt",
                          "hydrology_use_const_bmelt");
  config.scalar_from_option("hydrology_const_bmelt",
                            "hydrology_const_bmelt");
  config.scalar_from_option("hydrology_tillwat_max",
                            "hydrology_tillwat_max");
  config.scalar_from_option("hydrology_tillwat_decay_rate",
                            "hydrology_tillwat_decay_rate");
  config.scalar_from_option("hydrology_hydraulic_conductivity",
                            "hydrology_hydraulic_conductivity");
  config.scalar_from_option("hydrology_thickness_power_in_flux",
                            "hydrology_thickness_power_in_flux");
  config.scalar_from_option("hydrology_gradient_power_in_flux",
                            "hydrology_gradient_power_in_flux");
  // additional to RoutingHydrology, these apply to DistributedHydrology:
  config.scalar_from_option("hydrology_roughness_scale",
                            "hydrology_roughness_scale");
  config.scalar_from_option("hydrology_cavitation_opening_coefficient",
                            "hydrology_cavitation_opening_coefficient");
  config.scalar_from_option("hydrology_creep_closure_coefficient",
                            "hydrology_creep_closure_coefficient");
  config.scalar_from_option("hydrology_regularizing_porosity",
                            "hydrology_regularizing_porosity");

  // Time-stepping
  config.keyword_from_option("calendar", "calendar",
                             "standard,gregorian,proleptic_gregorian,noleap,365_day,360_day,julian,none");

  config.scalar_from_option("adapt_ratio",
                            "adaptive_timestepping_ratio");

  config.scalar_from_option("timestep_hit_multiples",
                            "timestep_hit_multiples");

  config.flag_from_option("count_steps", "count_time_steps");
  config.scalar_from_option("max_dt", "maximum_time_step_years");


  // SIA-related
  config.scalar_from_option("bed_smoother_range", "bed_smoother_range");

  config.keyword_from_option("gradient", "surface_gradient_method",
                             "eta,haseloff,mahaffy");

  // rheology-related
  config.scalar_from_option("sia_n", "sia_Glen_exponent");
  config.scalar_from_option("ssa_n", "ssa_Glen_exponent");

  config.scalar_from_option("sia_e", "sia_enhancement_factor");
  config.scalar_from_option("ssa_e", "ssa_enhancement_factor");

  config.flag_from_option("e_age_coupling", "e_age_coupling");

  // This parameter is used by the Goldsby-Kohlstedt flow law.
  config.scalar_from_option("ice_grain_size", "ice_grain_size");

  config.flag_from_option("grain_size_age_coupling",
                          "compute_grain_size_using_age");

  // SSA
  // Decide on the algorithm for solving the SSA
  config.keyword_from_option("ssa_method", "ssa_method", "fd,fem");

  config.scalar_from_option("ssa_eps",  "epsilon_ssa");
  config.scalar_from_option("ssa_maxi", "max_iterations_ssafd");
  config.scalar_from_option("ssa_rtol", "ssafd_relative_convergence");

  config.scalar_from_option("ssafd_nuH_iter_failure_underrelaxation", "ssafd_nuH_iter_failure_underrelaxation");

  config.flag_from_option("ssa_dirichlet_bc", "ssa_dirichlet_bc");
  config.flag_from_option("cfbc", "calving_front_stress_boundary_condition");

  // Basal sliding fiddles
  config.flag_from_option("brutal_sliding", "brutal_sliding");
  config.scalar_from_option("brutal_sliding_scale","brutal_sliding_scale");

  config.scalar_from_option("sliding_scale_factor_reduces_tauc",
                            "sliding_scale_factor_reduces_tauc");

  // SSA Inversion

  config.keyword_from_option("inv_method","inv_ssa_method",
                             "sd,nlcg,ign,tikhonov_lmvm,tikhonov_cg,tikhonov_blmvm,tikhonov_lcl,tikhonov_gn");

  config.keyword_from_option("inv_design_param",
                             "inv_design_param","ident,trunc,square,exp");

  config.scalar_from_option("inv_target_misfit","inv_target_misfit");

  config.scalar_from_option("tikhonov_penalty","tikhonov_penalty_weight");
  config.scalar_from_option("tikhonov_atol","tikhonov_atol");
  config.scalar_from_option("tikhonov_rtol","tikhonov_rtol");
  config.scalar_from_option("tikhonov_ptol","tikhonov_ptol");

  config.keyword_from_option("inv_state_func",
                             "inv_state_func",
                             "meansquare,log_ratio,log_relative");
  config.keyword_from_option("inv_design_func",
                             "inv_design_func","sobolevH1,tv");

  config.scalar_from_option("inv_design_cL2","inv_design_cL2");
  config.scalar_from_option("inv_design_cH1","inv_design_cH1");
  config.scalar_from_option("inv_ssa_tv_exponent","inv_ssa_tv_exponent");
  config.scalar_from_option("inv_log_ratio_scale","inv_log_ratio_scale");

  // Basal strength
  config.scalar_from_option("till_cohesion", "till_cohesion");
  config.scalar_from_option("till_reference_void_ratio",
                            "till_reference_void_ratio");
  config.scalar_from_option("till_compressibility_coefficient",
                            "till_compressibility_coefficient");
  config.scalar_from_option("till_effective_fraction_overburden",
                            "till_effective_fraction_overburden");
  config.scalar_from_option("till_log_factor_transportable_water",
                            "till_log_factor_transportable_water");

  bool topg_to_phi_set;
  std::vector<double> inarray(4);
  // read the comma-separated list of four values
  OptionsRealArray("-topg_to_phi", "phi_min, phi_max, topg_min, topg_max",
                   inarray, topg_to_phi_set);
  if (topg_to_phi_set == true) {
    if (inarray.size() != 4) {
      throw RuntimeError::formatted("option -topg_to_phi requires a comma-separated list with 4 numbers; got %d",
                                    inarray.size());
    }
    config.set_flag("till_use_topg_to_phi", true);
    config.set_double("till_topg_to_phi_phi_min", inarray[0]);
    config.set_double("till_topg_to_phi_phi_max", inarray[1]);
    config.set_double("till_topg_to_phi_topg_min",inarray[2]);
    config.set_double("till_topg_to_phi_topg_max",inarray[3]);
  }

  config.flag_from_option("tauc_slippery_grounding_lines",
                          "tauc_slippery_grounding_lines");
  config.flag_from_option("tauc_add_transportable_water",
                          "tauc_add_transportable_water");

  config.keyword_from_option("yield_stress", "yield_stress_model",
                             "constant,mohr_coulomb");

  // all basal strength models use this in ice-free areas
  config.scalar_from_option("high_tauc", "high_tauc");

  // controls regularization of plastic basal sliding law
  config.scalar_from_option("plastic_reg", "plastic_regularization");

  // "friction angle" in degrees. We allow -plastic_phi without an
  // argument: MohrCoulombYieldStress interprets that as "set
  // constant till friction angle using the default read from a config
  // file or an override file".
  bool plastic_phi_set = false;
  OptionsHasArgument("-plastic_phi", plastic_phi_set);
  if (plastic_phi_set == true) {
    config.scalar_from_option("plastic_phi", "default_till_phi");
  }

  // use pseudo plastic instead of pure plastic; see iMbasal.cc
  config.flag_from_option("pseudo_plastic", "do_pseudo_plastic_till");

  // power in denominator on pseudo_plastic_uthreshold; typical is q=0.25; q=0 is pure plastic
  config.scalar_from_option("pseudo_plastic_q", "pseudo_plastic_q");

  // threshold; at this velocity tau_c is basal shear stress
  config.scalar_from_option("pseudo_plastic_uthreshold",
                            "pseudo_plastic_uthreshold");

  config.flag_from_option("subgl", "sub_groundingline");

  // Ice shelves
  config.flag_from_option("part_grid", "part_grid");

  config.flag_from_option("part_grid_reduce_frontal_thickness",
                          "part_grid_reduce_frontal_thickness");

  config.flag_from_option("part_redist", "part_redist");

  config.scalar_from_option("nu_bedrock", "nu_bedrock");
  OptionsIsSet("-nu_bedrock", flag);
  if (flag) {
    config.set_flag_from_option("nu_bedrock_set", true);
  }

  // fracture density
  config.flag_from_option("fractures", "do_fracture_density");
  config.flag_from_option("write_fd_fields", "write_fd_fields");
  config.scalar_from_option("fracture_softening",
                            "fracture_density_softening_lower_limit");

  // Calving
  config.string_from_option("calving", "calving_methods");

  config.scalar_from_option("thickness_calving_threshold", "thickness_calving_threshold");

  // evaluates the adaptive timestep based on a CFL criterion with respect to the eigenCalving rate
  config.flag_from_option("cfl_eigen_calving", "cfl_eigen_calving");
  config.scalar_from_option("eigen_calving_K", "eigen_calving_K");

  config.flag_from_option("kill_icebergs", "kill_icebergs");

  // Output
  config.keyword_from_option("o_order", "output_variable_order",
                             "xyz,yxz,zyx");

  config.keyword_from_option("o_format", "output_format",
                             "netcdf3,quilt,netcdf4_parallel,pnetcdf,hdf5");

  config.scalar_from_option("summary_vol_scale_factor_log10",
                            "summary_vol_scale_factor_log10");
  config.scalar_from_option("summary_area_scale_factor_log10",
                            "summary_area_scale_factor_log10");

  // Metadata
  config.string_from_option("title", "run_title");
  config.string_from_option("institution", "institution");

  // Skipping
  config.flag_from_option("skip", "do_skip");
  config.scalar_from_option("skip_max", "skip_max");

  // Shortcuts

  // option "-pik" turns on a suite of PISMPIK effects (but NOT a calving choice,
  // and in particular NOT  "-calving eigen_calving")
  OptionsIsSet("-pik", "enable suite of PISM-PIK mechanisms", flag);
  if (flag) {
    config.set_flag_from_option("calving_front_stress_boundary_condition", true);
    config.set_flag_from_option("part_grid", true);
    config.set_flag_from_option("part_redist", true);
    config.set_flag_from_option("kill_icebergs", true);
    config.set_flag_from_option("sub_groundingline", true);
  }

  if (config.get_string("calving_methods").find("eigen_calving") != std::string::npos) {
    config.set_flag_from_option("part_grid", true);
    // eigen-calving requires a wider stencil:
    config.set_double("grid_max_stencil_width", 3);
  }

  // all calving mechanisms require iceberg removal
  if (config.get_string("calving_methods").empty() == false) {
    config.set_flag_from_option("kill_icebergs", true);
  }

  // kill_icebergs requires part_grid
  if (config.get_flag("kill_icebergs")) {
    config.set_flag_from_option("part_grid", true);
  }

  config.keyword_from_option("stress_balance", "stress_balance_model",
                             "none,prescribed_sliding,sia,ssa,prescribed_sliding+sia,ssa+sia");

  bool test_climate_models = false;
  OptionsIsSet("-test_climate_models", "Disable ice dynamics to test climate models",
               test_climate_models);
  if (test_climate_models) {
    config.set_string_from_option("stress_balance_model", "none");
    config.set_flag_from_option("do_energy", false);
    config.set_flag_from_option("do_age", false);
    // let the user decide if they want to use "-no_mass" or not
  }

  config.keyword_from_option("bed_def",
                             "bed_deformation_model", "none,iso,lc");
  config.flag_from_option("bed_def_lc_elastic_model", "bed_def_lc_elastic_model");

  config.flag_from_option("dry", "is_dry_simulation");

  config.flag_from_option("clip_shelf_base_salinity",
                          "ocean_three_equation_model_clip_salinity");

  config.scalar_from_option("meltfactor_pik", "ocean_pik_melt_factor");

  // old options
  check_old_option_and_stop("-sliding_scale_brutal",
                            "-brutal_sliding' and '-brutal_sliding_scale");
  check_old_option_and_stop("-ssa_sliding", "-stress_balance ...");
  check_old_option_and_stop("-ssa_floating_only", "-stress_balance ...");
  check_old_option_and_stop("-sia", "-stress_balance ...");
  check_old_option_and_stop("-no_sia", "-stress_balance ...");
  check_old_option_and_stop("-hold_tauc", "-yield_stress constant");
  check_old_option_and_stop("-ocean_kill", "-calving ocean_kill -ocean_kill_file foo.nc");
  check_old_option_and_stop("-eigen_calving", "-calving eigen_calving -eigen_calving_K XXX");
  check_old_option_and_stop("-calving_at_thickness",
                            "-calving thickness_calving -thickness_calving_threshold XXX");
  check_old_option_and_stop("-float_kill", "-calving float_kill");
  check_old_option_and_stop("-no_energy", "-energy none");
  check_old_option_and_stop("-cold", "-energy cold");

  return 0;
}

} // end of namespace pism

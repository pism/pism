// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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
  int       myLevel;
  PetscBool     verbose, levelSet;

  ierr = setVerbosityLevel(2);
  ierr = PetscOptionsGetInt(PETSC_NULL, "-verbose", &myLevel, &levelSet); CHKERRQ(ierr);
  if (levelSet == PETSC_TRUE) {
    ierr = setVerbosityLevel(myLevel);
  } else {
    ierr = PetscOptionsHasName(PETSC_NULL, "-verbose", &verbose); CHKERRQ(ierr);
    if (verbose == PETSC_TRUE)   ierr = setVerbosityLevel(3);
  }
  return 0;
}

//! Print a warning telling the user that an option was ignored.
PetscErrorCode ignore_option(MPI_Comm com, std::string name) {
  PetscErrorCode ierr;
  PetscBool option_is_set;

  char tmp[1]; // dummy string
  ierr = PetscOptionsGetString(PETSC_NULL, name.c_str(), tmp, 1, &option_is_set); CHKERRQ(ierr);

  if (option_is_set) {
    ierr = verbPrintf(1, com, "PISM WARNING: ignoring command-line option '%s'.\n",
                      name.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

//! Stop if an option `old_name` is set, printing a message that `new_name` should be used instead.
PetscErrorCode check_old_option_and_stop(MPI_Comm com, std::string old_name, std::string new_name) {
  PetscErrorCode ierr;
  PetscBool option_is_set;

  char tmp[1]; // dummy string
  ierr = PetscOptionsGetString(PETSC_NULL, old_name.c_str(), tmp, 1, &option_is_set); CHKERRQ(ierr);

  if (option_is_set) {
    ierr = PetscPrintf(com,
                       "PISM ERROR: command-line option '%s' is deprecated. Please use '%s' instead.\n",
                       old_name.c_str(), new_name.c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  return 0;
}

//!Stop if an option `name` is set.
PetscErrorCode stop_if_set(MPI_Comm com, std::string name) {
  PetscErrorCode ierr;
  PetscBool option_is_set;

  char tmp[1]; // dummy string
  ierr = PetscOptionsGetString(PETSC_NULL, name.c_str(), tmp, 1, &option_is_set); CHKERRQ(ierr);

  if (option_is_set) {
    ierr = PetscPrintf(com, "PISM ERROR: command-line option '%s' is not allowed.\n",
                       name.c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  return 0;
}

//! \brief Stop if -version is set.
PetscErrorCode stop_on_version_option() {
  PetscErrorCode ierr;

  bool vSet = false;
  ierr = PISMOptionsIsSet("-version", vSet); CHKERRQ(ierr);
  if (vSet == false)
    return 0;

  PISMEndQuiet();
  return 0;
}

//! \brief Print a usage message.
PetscErrorCode just_show_usage(MPI_Comm com, std::string execname, std::string usage) {
  PetscErrorCode ierr;
  ierr = verbPrintf(1,com,
      "%s is a PISM (http://www.pism-docs.org) executable.\nOptions cheat-sheet:\n",
      execname.c_str()); CHKERRQ(ierr);
  ierr = verbPrintf(1,com,usage.c_str()); CHKERRQ(ierr);
  ierr = verbPrintf(1,com,
      "Parallel run using N processes (typical case):  mpiexec -n N %s ...\n"
      "For more help with PISM:\n"
      "  1. download PDF User's Manual:\n"
      "       http://www.pism-docs.org/wiki/lib/exe/fetch.php?media=manual.pdf\n"
      "  2. read browser for technical details:\n"
      "       http://www.pism-docs.org/doxy/html/index.html\n"
      "  3. view issues/bugs at source host: https://github.com/pism/pism/issues\n"
      "  4. do '%s -help | grep foo' to see PISM and PETSc options with 'foo'.\n"
      "  5. email for help:  help@pism-docs.org\n",
      execname.c_str(), execname.c_str());  CHKERRQ(ierr);
  return 0;
}


//! @brief Show provided usage message and quit. (Consider using
//! show_usage_check_req_opts() in preference to this one.)
PetscErrorCode show_usage_and_quit(MPI_Comm com, std::string execname, std::string usage) {
  PetscErrorCode ierr;

  ierr = stop_on_version_option(); CHKERRQ(ierr);

  ierr = just_show_usage(com, execname, usage); CHKERRQ(ierr);

  PISMEndQuiet();
  return 0;
}


//! @brief In a single call a driver program can provide a usage string to
//! the user and check if required options are given, and if not, end.
PetscErrorCode show_usage_check_req_opts(MPI_Comm com, std::string execname,
                                         std::vector<std::string> required_options,
                                         std::string usage) {
  PetscErrorCode ierr;

  ierr = stop_on_version_option(); CHKERRQ(ierr);

  bool usageSet = false;
  ierr = PISMOptionsIsSet("-usage", usageSet); CHKERRQ(ierr);
  if (usageSet == true) {
    ierr = show_usage_and_quit(com, execname, usage); CHKERRQ(ierr);
  }

  // go through list of required options, and if not given, fail
  bool req_absent = false;
  for (size_t ii=0; ii < required_options.size(); ii++) {
    bool set = false;
    ierr = PISMOptionsIsSet(required_options[ii], set); CHKERRQ(ierr);
    if (set == PETSC_FALSE) {
      req_absent = true;
      ierr = verbPrintf(1,com,
        "PISM ERROR: option %s required\n",required_options[ii].c_str());
        CHKERRQ(ierr);
    }
  }
  if (req_absent == PETSC_TRUE) {
    ierr = verbPrintf(1,com,"\n"); CHKERRQ(ierr);
    ierr = show_usage_and_quit(com, execname, usage); CHKERRQ(ierr);
  }

  // show usage message with -help, but don't fail
  bool helpSet = false;
  ierr = PISMOptionsIsSet("-help", helpSet); CHKERRQ(ierr);
  if (helpSet == true) {
    ierr = just_show_usage(com, execname, usage); CHKERRQ(ierr);
  }

  return 0;
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
PetscErrorCode PISMOptionsList(MPI_Comm com, std::string opt, std::string description,
                               std::set<std::string> choices, std::string default_value,
                               std::string &result, bool &flag) {
  PetscErrorCode ierr;
  char tmp[TEMPORARY_STRING_LENGTH];
  std::string list, descr;
  PetscBool opt_set = PETSC_FALSE;

  if (choices.empty()) {
    SETERRQ(com, 1, "PISMOptionsList: empty choices argument");
  }

  std::set<std::string>::iterator j = choices.begin();
  list = "[" + *j++;
  while (j != choices.end()) {
    list += ", " + (*j++);
  }
  list += "]";

  descr = description + " Choose one of " + list;

  ierr = PetscOptionsString(opt.c_str(), descr.c_str(), "", default_value.c_str(),
                            tmp, TEMPORARY_STRING_LENGTH, &opt_set); CHKERRQ(ierr);

  // return the default value if the option was not set
  if (!opt_set) {
    flag = false;
    result = default_value;
    return 0;
  }

  std::string keyword = tmp;
  // find ":" and discard everything that goes after
  size_t n = keyword.find(":");
  if (n != std::string::npos)
    keyword.resize(n);

  // return the choice if it is valid and stop if it is not
  if (choices.find(keyword) != choices.end()) {
    flag = true;
    result = tmp;
  } else {
    ierr = PetscPrintf(com, "ERROR: invalid %s argument: \"%s\". Please choose one of %s.\n",
                       opt.c_str(), tmp, list.c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  return 0;
}

//! \brief Process a command-line option taking a string as an argument.
PetscErrorCode PISMOptionsString(std::string option, std::string text,
                                 std::string &result, bool &is_set, bool allow_empty_arg) {
  PetscErrorCode ierr;
  char tmp[TEMPORARY_STRING_LENGTH];
  PetscBool flag;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "",
                            result.c_str(), tmp,
                            TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);

  is_set = (flag == PETSC_TRUE);

  if (is_set) {
    if (strlen(tmp) == 0) {
      if (allow_empty_arg)
        result.clear();
      else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                           "ERROR: command line option '%s' requires an argument.\n",
                           option.c_str()); CHKERRQ(ierr);
        PISMEnd();
      }
    } else
      result = tmp;
  }

  return 0;
}

//! PISM wrapper replacing PetscOptionsStringArray.
PetscErrorCode PISMOptionsStringArray(std::string opt, std::string text, std::string default_value,
                                      std::vector<std::string>& result, bool &flag) {
  PetscErrorCode ierr;
  char tmp[TEMPORARY_STRING_LENGTH];
  PetscBool opt_set = PETSC_FALSE;

  ierr = PetscOptionsString(opt.c_str(), text.c_str(), "", default_value.c_str(),
                            tmp, TEMPORARY_STRING_LENGTH, &opt_set); CHKERRQ(ierr);

  result.clear();

  std::string word;
  if (opt_set) {
    std::istringstream arg(tmp);
    while (getline(arg, word, ','))
      result.push_back(word);

    if (result.empty()) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "ERROR: command line option '%s' requires an argument.\n",
                         opt.c_str()); CHKERRQ(ierr);
      PISMEnd();
    }

    flag = true;
  } else {                      // parse the default list given
    std::istringstream arg(default_value);
    while (getline(arg, word, ','))
      result.push_back(word);

    flag = false;
  }

  return 0;
}

//! Process a command-line option and return a set of strings.
PetscErrorCode PISMOptionsStringSet(std::string opt, std::string text, std::string default_value,
                                    std::set<std::string>& result, bool &flag) {
  std::vector<std::string> tmp;
  PetscErrorCode ierr;

  ierr = PISMOptionsStringArray(opt, text, default_value, tmp, flag); CHKERRQ(ierr);

  result.clear();
  std::vector<std::string>::iterator j = tmp.begin();
  while(j != tmp.end()) {
    result.insert(*j);
    ++j;
  }

  return 0;
}

//! \brief Process a command-line option taking an integer as an argument.
PetscErrorCode PISMOptionsInt(std::string option, std::string text,
                              int &result, bool &is_set) {
  PetscErrorCode ierr;
  char str[TEMPORARY_STRING_LENGTH];
  PetscBool flag;
  char *endptr;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "", "none", str,
                            TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);

  is_set = (flag == PETSC_TRUE);

  if (is_set == false)
    return 0;

  if (strlen(str) == 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "ERROR: command line option '%s' requires an argument.\n",
                       option.c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  result = (int) strtol(str, &endptr, 10);
  if (*endptr != '\0') {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "PISM ERROR: Can't parse \"%s %s\": (%s is not an integer).\n",
                       option.c_str(), str, str); CHKERRQ(ierr);
    PISMEnd();
  }

  return 0;
}

//! \brief Process a command-line option taking a real number as an argument.
PetscErrorCode PISMOptionsReal(std::string option, std::string text,
                               double &result, bool &is_set) {
  PetscErrorCode ierr;
  char str[TEMPORARY_STRING_LENGTH];
  PetscBool flag;
  char *endptr;

  memset(str, 0, TEMPORARY_STRING_LENGTH);
  snprintf(str, TEMPORARY_STRING_LENGTH, "%f", result);

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "", str, str,
                            TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);

  is_set = (flag == PETSC_TRUE);

  if (is_set == false)
    return 0;

  if (strlen(str) == 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "ERROR: command line option '%s' requires an argument.\n",
                       option.c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  result = strtod(str, &endptr);
  if (*endptr != '\0') {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "PISM ERROR: Can't parse \"%s %s\": (%s is not a number).\n",
                       option.c_str(), str, str); CHKERRQ(ierr);
    PISMEnd();
  }

  return 0;
}
//! \brief Process a command-line option taking a comma-separated list of reals
//! as an argument.
PetscErrorCode PISMOptionsRealArray(std::string option, std::string text,
                                    std::vector<double> &result, bool &is_set) {
  PetscErrorCode ierr;
  char str[TEMPORARY_STRING_LENGTH];
  PetscBool flag;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "",
                            "none", str,
                            TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);

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
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                           "PISM ERROR: Can't parse %s (%s is not a number).\n",
                           tmp.c_str(), tmp.c_str()); CHKERRQ(ierr);
        PISMEnd();
      }
      else
        result.push_back(d);
    }
  }

  return 0;
}

//! \brief Process a command-line option taking a comma-separated list of
//! integers as an argument.
PetscErrorCode PISMOptionsIntArray(std::string option, std::string text,
                                    std::vector<int> &result, bool &is_set) {
  PetscErrorCode ierr;
  std::vector<double> tmp;

  ierr = PISMOptionsRealArray(option, text, tmp, is_set); CHKERRQ(ierr);

  result.clear();
  for (unsigned int j = 0; j < tmp.size(); ++j)
    result.push_back(static_cast<int>(tmp[j]));

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
PetscErrorCode PISMOptionsIsSet(std::string option, bool &result) {
  PetscErrorCode ierr;
  char tmp[1];
  PetscBool flag;

  ierr = PetscOptionsGetString(PETSC_NULL, option.c_str(), tmp, 1, &flag); CHKERRQ(ierr);

  result = (flag == PETSC_TRUE);

  return 0;
}

//! A version of PISMOptionsIsSet that prints a -help message.
PetscErrorCode PISMOptionsIsSet(std::string option, std::string text,
                                bool &result) {
  PetscErrorCode ierr;
  char tmp[1];
  PetscBool flag;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "",
                            "", tmp, 1, &flag); CHKERRQ(ierr);

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
PetscErrorCode PISMOptionsHasArgument(std::string option, bool &result) {
  PetscErrorCode ierr;
  std::string tmp;

  ierr = PISMOptionsString(option, "", tmp, result, true); CHKERRQ(ierr);

  if (result == false || tmp.empty() == true)
    result = false;

  return 0;
}

//! Initializes the config parameter and flag database.
/*!
  Processes -config and -config_override command line options.
 */
PetscErrorCode init_config(MPI_Comm com,
                           PISMConfig &config, PISMConfig &overrides,
                           bool process_options) {
  PetscErrorCode ierr;

  std::string alt_config = PISM_DefaultConfigFile,
    override_config;
  bool use_alt_config, use_override_config;

  ierr = PetscOptionsBegin(com, "", "PISM config file options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-config", "Specifies the name of an alternative config file",
                             alt_config, use_alt_config); CHKERRQ(ierr);
    ierr = PISMOptionsString("-config_override", "Specifies a config override file name",
                             override_config, use_override_config); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = config.read(alt_config); CHKERRQ(ierr);

  if (use_override_config) {
    ierr = overrides.read(override_config); CHKERRQ(ierr);
    config.import_from(overrides);
    ierr = verbPrintf(2, com, "CONFIG OVERRIDES read from file '%s'.\n",
                      override_config.c_str()); CHKERRQ(ierr);
  }

  if (process_options) {
    ierr =  set_config_from_options(com, config); CHKERRQ(ierr);
  }

  config.print_to_stdout();

  return 0;
}

PetscErrorCode set_config_from_options(MPI_Comm com, PISMConfig &config) {
  PetscErrorCode ierr;
  bool flag;

  // Energy modeling
  ierr = config.flag_from_option("varc", "use_linear_in_temperature_heat_capacity");  CHKERRQ(ierr);
  ierr = config.flag_from_option("vark",
                                 "use_temperature_dependent_thermal_conductivity");  CHKERRQ(ierr);

  ierr = config.flag_from_option("bmr_in_cont", "include_bmr_in_continuity"); CHKERRQ(ierr);

  {
    bool energy_set;
    std::string energy;
    std::set<std::string> choices;
    choices.insert("none");
    choices.insert("cold");
    choices.insert("enthalpy");

    ierr = PISMOptionsList(com, "-energy",
                           "choose the energy model (one of 'none', 'cold', 'enthalpy')",
                           choices, "enthalpy", energy, energy_set); CHKERRQ(ierr);

    if (energy_set == true) {
      if (energy == "none") {
        ierr = config.set_flag_from_option("do_energy", false); CHKERRQ(ierr);
        // Allow selecting cold ice flow laws in isothermal mode. 
        ierr = config.set_flag_from_option("do_cold_ice_methods", true); CHKERRQ(ierr);
      } else if (energy == "cold") {
        ierr = config.set_flag_from_option("do_energy", true); CHKERRQ(ierr);
        ierr = config.set_flag_from_option("do_cold_ice_methods", true); CHKERRQ(ierr);
      } else if (energy == "enthalpy") {
        ierr = config.set_flag_from_option("do_energy", true); CHKERRQ(ierr);
        ierr = config.set_flag_from_option("do_cold_ice_methods", false); CHKERRQ(ierr);
      } else {
        // can't happen (PISMOptionsList validates its input)
        assert(false);
      }
    }
  }

  // at bootstrapping, choose whether the method uses smb as upper boundary for
  // vertical velocity
  ierr = config.keyword_from_option("boot_temperature_heuristic",
                                    "bootstrapping_temperature_heuristic", "smb,quartic_guess"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("low_temp", "global_min_allowed_temp"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("max_low_temps", "max_low_temp_count"); CHKERRQ(ierr);

  // Sub-models
  ierr = config.flag_from_option("age", "do_age"); CHKERRQ(ierr);
  ierr = config.flag_from_option("mass", "do_mass_conserve"); CHKERRQ(ierr);

  // hydrology
  ierr = config.keyword_from_option("hydrology", "hydrology_model",
                                    "null,routing,distributed"); CHKERRQ(ierr);
  ierr = config.flag_from_option("hydrology_use_const_bmelt",
                                 "hydrology_use_const_bmelt"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("hydrology_const_bmelt",
                                   "hydrology_const_bmelt"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("hydrology_tillwat_max",
                                   "hydrology_tillwat_max"); CHKERRQ(ierr);
  // this only applies to PISMNullTransportHydrology:
  ierr = config.scalar_from_option("hydrology_tillwat_decay_rate_null",
                                   "hydrology_tillwat_decay_rate_null"); CHKERRQ(ierr);
  // these only apply to PISMRoutingHydrology:
  ierr = config.scalar_from_option("hydrology_tillwat_rate",
                                   "hydrology_tillwat_rate"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("hydrology_tillwat_transfer_proportion",
                                   "hydrology_tillwat_transfer_proportion"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("hydrology_hydraulic_conductivity",
                                   "hydrology_hydraulic_conductivity"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("hydrology_thickness_power_in_flux",
                                   "hydrology_thickness_power_in_flux"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("hydrology_gradient_power_in_flux",
                                   "hydrology_gradient_power_in_flux"); CHKERRQ(ierr);
  // additional to PISMRoutingHydrology, these apply to PISMDistributedHydrology:
  ierr = config.scalar_from_option("hydrology_roughness_scale",
                                   "hydrology_roughness_scale"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("hydrology_cavitation_opening_coefficient",
                                   "hydrology_cavitation_opening_coefficient"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("hydrology_creep_closure_coefficient",
                                   "hydrology_creep_closure_coefficient"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("hydrology_regularizing_porosity",
                                   "hydrology_regularizing_porosity"); CHKERRQ(ierr);

  // Time-stepping
  ierr = config.keyword_from_option("calendar", "calendar",
                                    "standard,gregorian,proleptic_gregorian,noleap,365_day,360_day,julian,none"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("adapt_ratio",
                                   "adaptive_timestepping_ratio"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("timestep_hit_multiples",
                                   "timestep_hit_multiples"); CHKERRQ(ierr);

  ierr = config.flag_from_option("count_steps", "count_time_steps"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("max_dt", "maximum_time_step_years"); CHKERRQ(ierr);


  // SIA
  ierr = config.scalar_from_option("bed_smoother_range", "bed_smoother_range"); CHKERRQ(ierr);

  ierr = config.keyword_from_option("gradient", "surface_gradient_method",
                                    "eta,haseloff,mahaffy"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("sia_e", "sia_enhancement_factor"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_e", "ssa_enhancement_factor"); CHKERRQ(ierr);

  ierr = config.flag_from_option("e_age_coupling", "e_age_coupling"); CHKERRQ(ierr);

  // This parameter is used by the Goldsby-Kohlstedt flow law.
  ierr = config.scalar_from_option("ice_grain_size", "ice_grain_size"); CHKERRQ(ierr);

  ierr = config.flag_from_option("grain_size_age_coupling",
                                 "compute_grain_size_using_age"); CHKERRQ(ierr);

  // SSA
  // Decide on the algorithm for solving the SSA
  ierr = config.keyword_from_option("ssa_method", "ssa_method", "fd,fem"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("ssa_eps",  "epsilon_ssa"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_maxi", "max_iterations_ssafd"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_rtol", "ssafd_relative_convergence"); CHKERRQ(ierr);

  ierr = config.flag_from_option("ssa_dirichlet_bc", "ssa_dirichlet_bc"); CHKERRQ(ierr);
  ierr = config.flag_from_option("cfbc", "calving_front_stress_boundary_condition"); CHKERRQ(ierr);

  // Basal sliding fiddles
  ierr = config.flag_from_option("brutal_sliding", "brutal_sliding"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("brutal_sliding_scale","brutal_sliding_scale"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("sliding_scale_factor_reduces_tauc",
                                   "sliding_scale_factor_reduces_tauc"); CHKERRQ(ierr);

  // SSA Inversion

  ierr = config.keyword_from_option("inv_method","inv_ssa_method",
                                    "sd,nlcg,ign,tikhonov_lmvm,tikhonov_cg,tikhonov_blmvm,tikhonov_lcl,tikhonov_gn");
  CHKERRQ(ierr);

  ierr = config.keyword_from_option("inv_design_param",
                                    "inv_design_param","ident,trunc,square,exp"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("inv_target_misfit","inv_target_misfit"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("tikhonov_penalty","tikhonov_penalty_weight"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("tikhonov_atol","tikhonov_atol"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("tikhonov_rtol","tikhonov_rtol"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("tikhonov_ptol","tikhonov_ptol"); CHKERRQ(ierr);

  ierr = config.keyword_from_option("inv_state_func",
                                    "inv_state_func",
                                    "meansquare,log_ratio,log_relative"); CHKERRQ(ierr);
  ierr = config.keyword_from_option("inv_design_func",
                                    "inv_design_func","sobolevH1,tv"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("inv_design_cL2","inv_design_cL2"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("inv_design_cH1","inv_design_cH1"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("inv_ssa_tv_exponent","inv_ssa_tv_exponent"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("inv_log_ratio_scale","inv_log_ratio_scale"); CHKERRQ(ierr);

  // Basal strength
  ierr = config.scalar_from_option("till_c_0", "till_c_0");      CHKERRQ(ierr);
  ierr = config.scalar_from_option("till_reference_void_ratio",
                                   "till_reference_void_ratio");      CHKERRQ(ierr);
  ierr = config.scalar_from_option("till_compressibility_coefficient",
                                   "till_compressibility_coefficient");      CHKERRQ(ierr);
  ierr = config.scalar_from_option("till_effective_fraction_overburden",
                                   "till_effective_fraction_overburden");      CHKERRQ(ierr);
  ierr = config.flag_from_option("tauc_slippery_grounding_lines",
                                 "tauc_slippery_grounding_lines"); CHKERRQ(ierr);

  ierr = config.keyword_from_option("yield_stress", "yield_stress_model",
                                    "constant,mohr_coulomb"); CHKERRQ(ierr);

  // all basal strength models use this in ice-free areas
  ierr = config.scalar_from_option("high_tauc", "high_tauc");      CHKERRQ(ierr);

  // controls regularization of plastic basal sliding law
  ierr = config.scalar_from_option("plastic_reg", "plastic_regularization"); CHKERRQ(ierr);

  // "friction angle" in degrees. We allow -plastic_phi without an
  // argument: PISMMohrCoulombYieldStress interprets that as "set
  // constant till friction angle using the default read from a config
  // file or an override file".
  bool plastic_phi_set = false;
  ierr = PISMOptionsHasArgument("-plastic_phi", plastic_phi_set); CHKERRQ(ierr);
  if (plastic_phi_set == true) {
    ierr = config.scalar_from_option("plastic_phi", "default_till_phi"); CHKERRQ(ierr);
  }

  // use pseudo plastic instead of pure plastic; see iMbasal.cc
  ierr = config.flag_from_option("pseudo_plastic", "do_pseudo_plastic_till"); CHKERRQ(ierr);

  // power in denominator on pseudo_plastic_uthreshold; typical is q=0.25; q=0 is pure plastic
  ierr = config.scalar_from_option("pseudo_plastic_q", "pseudo_plastic_q"); CHKERRQ(ierr);

  // threshold; at this velocity tau_c is basal shear stress
  ierr = config.scalar_from_option("pseudo_plastic_uthreshold",
                                   "pseudo_plastic_uthreshold"); CHKERRQ(ierr);

  ierr = config.flag_from_option("subgl", "sub_groundingline"); CHKERRQ(ierr);
  ierr = config.flag_from_option("subgl_basal_melt", "sub_groundingline_basal_melt"); CHKERRQ(ierr);

  // sub_groundingline is implied by sub_groundingline_basal_melt
  if (config.get_flag("sub_groundingline_basal_melt")) {
    ierr = config.set_flag_from_option("sub_groundingline", true); CHKERRQ(ierr);
  }

  // Ice shelves

  ierr = config.flag_from_option("part_grid", "part_grid"); CHKERRQ(ierr);

  ierr = config.flag_from_option("part_grid_reduce_frontal_thickness",
                                 "part_grid_reduce_frontal_thickness"); CHKERRQ(ierr);

  ierr = config.flag_from_option("part_redist", "part_redist"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("nu_bedrock", "nu_bedrock"); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-nu_bedrock", flag);  CHKERRQ(ierr);
  if (flag) {
    config.set_flag_from_option("nu_bedrock_set", true);
  }

  // fracture density
  ierr = config.flag_from_option("fractures", "do_fracture_density"); CHKERRQ(ierr);
  ierr = config.flag_from_option("write_fd_fields", "write_fd_fields"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("fracture_softening",
                                   "fracture_density_softening_lower_limit"); CHKERRQ(ierr);


  // Calving

  ierr = config.string_from_option("calving", "calving_methods"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("thickness_calving_threshold", "thickness_calving_threshold"); CHKERRQ(ierr);

  // evaluates the adaptive timestep based on a CFL criterion with respect to the eigenCalving rate
  ierr = config.flag_from_option("cfl_eigen_calving", "cfl_eigen_calving"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("eigen_calving_K", "eigen_calving_K"); CHKERRQ(ierr);

  ierr = config.flag_from_option("kill_icebergs", "kill_icebergs"); CHKERRQ(ierr);

  // Output
  ierr = config.keyword_from_option("o_order", "output_variable_order",
                                    "xyz,yxz,zyx"); CHKERRQ(ierr);

  ierr = config.keyword_from_option("o_format", "output_format",
                                    "netcdf3,quilt,netcdf4_parallel,pnetcdf"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("summary_vol_scale_factor_log10",
                                   "summary_vol_scale_factor_log10"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("summary_area_scale_factor_log10",
                                   "summary_area_scale_factor_log10"); CHKERRQ(ierr);

  // Metadata
  ierr = config.string_from_option("title", "run_title"); CHKERRQ(ierr);
  ierr = config.string_from_option("institution", "institution"); CHKERRQ(ierr);

  // Skipping
  ierr = config.flag_from_option("skip", "do_skip"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("skip_max", "skip_max"); CHKERRQ(ierr);

  // Shortcuts

  // option "-pik" turns on a suite of PISMPIK effects (but not "-calving eigen_calving")
  ierr = PISMOptionsIsSet("-pik", "enable suite of PISM-PIK mechanisms", flag); CHKERRQ(ierr);
  if (flag) {
    config.set_flag_from_option("calving_front_stress_boundary_condition", true);
    config.set_flag_from_option("part_grid", true);
    config.set_flag_from_option("part_redist", true);
    config.set_flag_from_option("kill_icebergs", true);
  }

  if (config.get_string("calving_methods").find("eigen_calving") != std::string::npos) {
    config.set_flag_from_option("part_grid", true);
  }

  // all calving mechanisms require iceberg removal
  if (config.get_string("calving_methods").empty() == false) {
    config.set_flag_from_option("kill_icebergs", true);
  }

  // kill_icebergs requires part_grid
  if (config.get_flag("kill_icebergs")) {
    config.set_flag_from_option("part_grid", true);
  }

  ierr = config.keyword_from_option("stress_balance", "stress_balance_model",
                                    "none,prescribed_sliding,sia,ssa,prescribed_sliding+sia,ssa+sia"); CHKERRQ(ierr);

  bool test_climate_models = false;
  ierr = PISMOptionsIsSet("-test_climate_models", "Disable ice dynamics to test climate models",
                          test_climate_models); CHKERRQ(ierr);
  if (test_climate_models) {
    config.set_string_from_option("stress_balance_model", "none");
    config.set_flag_from_option("do_energy", false);
    config.set_flag_from_option("do_age", false);
    // let the user decide if they want to use "-no_mass" or not
  }

  ierr = config.flag_from_option("bed_def_lc_elastic_model", "bed_def_lc_elastic_model"); CHKERRQ(ierr);

  ierr = config.flag_from_option("dry", "is_dry_simulation"); CHKERRQ(ierr);

  ierr = config.flag_from_option("clip_shelf_base_salinity",
                                 "ocean_three_equation_model_clip_salinity"); CHKERRQ(ierr);

  // old options
  ierr = check_old_option_and_stop(com, "-sliding_scale_brutal",
                                   "-brutal_sliding' and '-brutal_sliding_scale"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-ssa_sliding", "-stress_balance ..."); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-ssa_floating_only", "-stress_balance ..."); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-sia", "-stress_balance ..."); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-no_sia", "-stress_balance ..."); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-hold_tauc", "-yield_stress constant"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-ocean_kill", "-calving ocean_kill -ocean_kill_file foo.nc"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-eigen_calving", "-calving eigen_calving -eigen_calving_K XXX"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-calving_at_thickness",
                                   "-calving thickness_calving -thickness_calving_threshold XXX"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-float_kill", "-calving float_kill"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-no_energy", "-energy none"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(com, "-cold", "-energy cold"); CHKERRQ(ierr);

  return 0;
}

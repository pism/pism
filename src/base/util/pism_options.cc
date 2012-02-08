// Copyright (C) 2011, 2012 PISM Authors
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

#include <sstream>
#include <algorithm>

#include "pism_options.hh"
#include "NCVariable.hh"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

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
PetscErrorCode ignore_option(MPI_Comm com, const char name[]) {
  PetscErrorCode ierr;
  PetscBool option_is_set;

  char tmp[1]; // dummy string
  ierr = PetscOptionsGetString(PETSC_NULL, name, tmp, 1, &option_is_set); CHKERRQ(ierr);

  if (option_is_set) {
    ierr = verbPrintf(1, com, "PISM WARNING: ignoring command-line option '%s'.\n",
		      name); CHKERRQ(ierr);
  }

  return 0;
}

//! Stop if an option \c old_name is set, printing a message that \c new_name should be used instead.
PetscErrorCode check_old_option_and_stop(MPI_Comm com, const char old_name[], const char new_name[]) {
  PetscErrorCode ierr;
  PetscBool option_is_set;

  char tmp[1]; // dummy string
  ierr = PetscOptionsGetString(PETSC_NULL, old_name, tmp, 1, &option_is_set); CHKERRQ(ierr);

  if (option_is_set) {
    ierr = PetscPrintf(com, "PISM ERROR: command-line option '%s' is deprecated. Please use '%s' instead.\n",
		       old_name, new_name); CHKERRQ(ierr);
    PISMEnd();
  }

  return 0;
}

//!Stop if an option \c name is set.
PetscErrorCode stop_if_set(MPI_Comm com, const char name[]) {
  PetscErrorCode ierr;
  PetscBool option_is_set;

  char tmp[1]; // dummy string
  ierr = PetscOptionsGetString(PETSC_NULL, name, tmp, 1, &option_is_set); CHKERRQ(ierr);

  if (option_is_set) {
    ierr = PetscPrintf(com, "PISM ERROR: command-line option '%s' is not allowed.\n",
		       name); CHKERRQ(ierr);
    PISMEnd();
  }

  return 0;
}


//! Parses a MATLAB-style range (a:delta:b).
PetscErrorCode parse_range(MPI_Comm com, string str, double *a, double *delta, double *b, string &keyword) {
  PetscErrorCode ierr;
  istringstream arg(str);
  vector<string> numbers;
  double doubles[3];

  // Split the string:
  string tmp;
  while (getline(arg, tmp, ':'))
    numbers.push_back(tmp);

  // Check if we have 3 numbers:
  if (numbers.size() != 3) {
      ierr = PetscPrintf(com,
			 "PISM ERROR: A range has to consist of exactly three numbers, separated by colons.\n");
      CHKERRQ(ierr);
      return 1;
  }

  keyword = "simple";

  // Convert each number from a string to double:
  for (int j = 0; j < 3; ++j) {
    double d;
    char *endptr;

    // take care of daily, monthly and yearly keywords:
    if (j == 1 && (numbers[j] == "daily" ||
                   numbers[j] == "monthly" ||
                   numbers[j] == "yearly")) {
      keyword = numbers[j];
      doubles[j] = 0;
      continue;
    }

    d = strtod(numbers[j].c_str(), &endptr);
    if (*endptr != '\0') {
      ierr = PetscPrintf(com, "PISM ERROR: Can't parse %s (%s is not a number).\n",
			 str.c_str(), numbers[j].c_str()); CHKERRQ(ierr);
      return 1;
    }
    else
      doubles[j] = d;
  }

  if (a) *a = doubles[0];
  if (delta) *delta = doubles[1];
  if (b) *b = doubles[2];

  return 0;
}

vector<double> compute_times(MPI_Comm com, const NCConfigVariable &config,
                             int a, int b, string keyword) {
  utUnit unit;
  string unit_str = "seconds since " + config.get_string("reference_date");
  vector<double> result;
  double a_offset, b_offset;

  // scan the units:
  int err = utScan(unit_str.c_str(), &unit);
  if (err != 0) {
    PetscPrintf(com, "PISM ERROR: invalid units specification: %s\n",
                unit_str.c_str());
    PISMEnd();
  }

  // get the 'year' out of the reference date:
  int reference_year, tmp1;
  float tmp2;
  utCalendar(0, &unit, &reference_year,
             &tmp1, &tmp1, &tmp1, &tmp1, &tmp2);

  // compute the number of seconds-since-the-reference date 'a' corresponds to:
  utInvCalendar(reference_year + a, // year
                1, 1,               // month, day
                0, 0, 0,            // hour, minute, second
                &unit,
                &a_offset);

  // compute the number of seconds-since-the-reference date 'b' corresponds to:
  utInvCalendar(reference_year + b, // year
                1, 1,               // month, day
                0, 0, 0,            // hour, minute, second
                &unit,
                &b_offset);

  if (keyword == "daily") {

    double t = a_offset, delta = 60*60*24; // seconds per day
    int year;

    do {
      result.push_back(t);
      t += delta;
      utCalendar(t, &unit, &year, &tmp1, &tmp1, &tmp1, &tmp1, &tmp2);
    } while (year <= reference_year + b);

    // add the last record:
    result.push_back(t);

  } else if (keyword == "monthly") {

    double t;
    int y, m;
    for (y = a; y <= b; y++) {
      for (m = 1; m <= 12; m++) {
        utInvCalendar(reference_year + y,   // year
                      m, 1,                 // month, day
                      0, 0, 0,              // hour, minute, second
                      &unit,
                      &t);
        result.push_back(t);
      }
    }

    // add the last record:
    utInvCalendar(reference_year + b + 1,   // year
                  1, 1,                     // month, day
                  0, 0, 0,                  // hour, minute, second
                  &unit,
                  &t);
    result.push_back(t);

  } else if (keyword == "yearly") {

    double t;
    for (int y = a; y <= b+1; y++) {    // note the "b + 1"
      utInvCalendar(reference_year + y,   // year
                    1, 1,                 // month, day
                    0, 0, 0,              // hour, minute, second
                    &unit,
                    &t);
      result.push_back(t);
    }

  } else {
    PetscPrintf(com,
                "PISM ERROR: unknown time-step keyword: %s\n"
                "            (only 'daily', 'monthly' and 'yearly' are implemented).\n",
                keyword.c_str());
    PISMEnd();
  }
  
  return result;
}

//! Parses a time specification.
/*!
  If it is a MATLAB-style range, then calls parse_range and computes all the points.

  If it is a comma-separated list, converts to double (with error-checking).
 */
PetscErrorCode parse_times(MPI_Comm com, const NCConfigVariable &config, string str, vector<double> &result) {
  PetscErrorCode ierr;
  int N;

  if (str.find(':') != string::npos) { // it's a range specification
    
    double a, delta, b;
    string keyword;
    ierr = parse_range(com, str, &a, &delta, &b, keyword);
    if (ierr != 0) return 1;

    if (a >= b) {
      ierr = PetscPrintf(com, "PISM ERROR: a >= b in the range specification %s.\n",
			 str.c_str()); CHKERRQ(ierr);
      return 1;
    }

    if (keyword != "simple") {

      result = compute_times(com, config, (int)a, (int)b, keyword);

    } else {
      if (delta <= 0) {
        ierr = PetscPrintf(com, "PISM ERROR: delta <= 0 in the range specification %s.\n",
                           str.c_str()); CHKERRQ(ierr);
        return 1;
      }

      N = (int)floor((b - a)/delta) + 1; // number of points in the range
      result.resize(N);

      for (int j = 0; j < N; ++j)
        result[j] = convert(a + delta*j, "years", "seconds");
    }

  } else {			// it's a list of times
    istringstream arg(str);
    string tmp;

    result.clear();
    while(getline(arg, tmp, ',')) {
      double d;
      char *endptr;

      d = strtod(tmp.c_str(), &endptr);
      if (*endptr != '\0') {
	ierr = PetscPrintf(com, "PISM ERROR: Can't parse %s (%s is not a number).\n",
			   str.c_str(), tmp.c_str()); CHKERRQ(ierr);
	return 1;
      }
      else
	result.push_back(convert(d, "years", "seconds"));
    }
    sort(result.begin(), result.end());
  }

  return 0;
}

//! \brief Stop if at least one of -version and -pismversion is set.
PetscErrorCode stop_on_version_option() {
  PetscErrorCode ierr;

  bool vSet = false, pvSet = false;
  ierr = PISMOptionsIsSet("-version", vSet); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-pismversion", pvSet); CHKERRQ(ierr);
  if ((vSet == false) && (pvSet == false))
    return 0;

  PISMEndQuiet();
  return 0;
}

//! \brief Print a usage message.
PetscErrorCode just_show_usage(
    MPI_Comm com, const char execname[], const char usage[]) {
  PetscErrorCode ierr;
  ierr = verbPrintf(1,com,
      "%s is a PISM (http://www.pism-docs.org) executable.\nOptions cheat-sheet:\n",
      execname); CHKERRQ(ierr);
  ierr = verbPrintf(1,com,usage); CHKERRQ(ierr);
  ierr = verbPrintf(1,com,
      "Parallel run using N processes (typical case):  mpiexec -n N %s ...\n"
      "For more help with PISM:\n"
      "  1. download PDF User's Manual:\n"
      "       http://www.pism-docs.org/wiki/lib/exe/fetch.php?media=manual.pdf\n"
      "  2. read browser for technical details:\n"
      "       http://www.pism-docs.org/doxy/html/index.html\n"
      "  3. search bugs and tasks source host: https://gna.org/projects/pism\n"
      "  4. do '%s -help | grep foo' to see PISM and PETSc options with 'foo'.\n"
      "  5. email:  help@pism-docs.org\n", 
      execname,execname);  CHKERRQ(ierr);
  return 0;
}


//! Show provided usage message and quit.  (Consider using show_usage_check_req_opts() in preference to this one.)
PetscErrorCode show_usage_and_quit(
    MPI_Comm com, const char execname[], const char usage[]) {
  PetscErrorCode ierr;

  ierr = stop_on_version_option(); CHKERRQ(ierr);

  ierr = just_show_usage(com, execname, usage); CHKERRQ(ierr);

  PISMEndQuiet();
  return 0;
}


//! In a single call a driver program can provide a usage string to the user and check if required options are given, and if not, end.
PetscErrorCode show_usage_check_req_opts(
    MPI_Comm com, const char execname[], vector<string> required_options,
    const char usage[]) {
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
PetscErrorCode PISMOptionsList(MPI_Comm com, string opt, string description, set<string> choices,
			       string default_value, string &result, bool &flag) {
  PetscErrorCode ierr;
  char tmp[TEMPORARY_STRING_LENGTH];
  string list, descr;
  PetscBool opt_set = PETSC_FALSE;

  if (choices.empty()) {
    SETERRQ(com, 1, "PISMOptionsList: empty choices argument");
  }

  set<string>::iterator j = choices.begin();
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

  // return the choice if it is valid and stop if it is not
  if (choices.find(tmp) != choices.end()) {
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
PetscErrorCode PISMOptionsString(string option, string text,
				 string &result, bool &is_set, bool allow_empty_arg) {
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
PetscErrorCode PISMOptionsStringArray(string opt, string text, string default_value,
				      vector<string>& result, bool &flag) {
  PetscErrorCode ierr;
  char tmp[TEMPORARY_STRING_LENGTH];
  PetscBool opt_set = PETSC_FALSE;

  ierr = PetscOptionsString(opt.c_str(), text.c_str(), "", default_value.c_str(),
			    tmp, TEMPORARY_STRING_LENGTH, &opt_set); CHKERRQ(ierr);

  result.clear();

  string word;
  if (opt_set) {
    istringstream arg(tmp);
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
    istringstream arg(default_value);
    while (getline(arg, word, ','))
      result.push_back(word);

    flag = false;
  }

  return 0;
}

//! \brief Process a command-line option taking an integer as an argument.
PetscErrorCode PISMOptionsInt(string option, string text,
			      PetscInt &result, bool &is_set) {
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
                       "PISM ERROR: Can't parse \"%s %s\": (%s is not a number).\n",
                       option.c_str(), str, str); CHKERRQ(ierr);
    PISMEnd();
  }

  return 0;
}

//! \brief Process a command-line option taking a real number as an argument.
PetscErrorCode PISMOptionsReal(string option, string text,
			       PetscReal &result, bool &is_set) {
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
PetscErrorCode PISMOptionsRealArray(string option, string text,
				    vector<PetscReal> &result, bool &is_set) {
  PetscErrorCode ierr;
  char str[TEMPORARY_STRING_LENGTH];
  PetscBool flag;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "",
			    "none", str,
			    TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);

  is_set = (flag == PETSC_TRUE);

  if (is_set) {
    istringstream arg(str);
    string tmp;

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
PetscErrorCode PISMOptionsIntArray(string option, string text,
				    vector<PetscInt> &result, bool &is_set) {
  PetscErrorCode ierr;
  vector<PetscReal> tmp;

  ierr = PISMOptionsRealArray(option, text, tmp, is_set); CHKERRQ(ierr);

  result.clear();
  for (unsigned int j = 0; j < tmp.size(); ++j)
    result.push_back(static_cast<PetscInt>(tmp[j]));

  return 0;
}

//! Checks if an option is present in the PETSc option database.
/*!

  This is (essentially) a reimplementation of PetscOptionsHasName, except that
  this *always* sets \c flag to PETSC_TRUE if an option is present.

  PetscOptionsHasName, on the other hand, sets \c flag to PETSC_FALSE if an
  option was set as "-foo FALSE", "-foo NO" or "-foo 0". Note that if one uses
  "-foo 0.0", PetscOptionsHasName will set \c flag to PETSC_TRUE.

  This unpredictability is bad. We want a function that does not depend on the
  argument given with an option.
 */
PetscErrorCode PISMOptionsIsSet(string option, bool &result) {
  PetscErrorCode ierr;
  char tmp[1];
  PetscBool flag;

  ierr = PetscOptionsGetString(PETSC_NULL, option.c_str(), tmp, 1, &flag); CHKERRQ(ierr);

  result = (flag == PETSC_TRUE);

  return 0;
}

//! A version of PISMOptionsIsSet that prints a -help message.
PetscErrorCode PISMOptionsIsSet(string option, string text,
				bool &result) {
  PetscErrorCode ierr;
  char tmp[1];
  PetscBool flag;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "",
			    "", tmp, 1, &flag); CHKERRQ(ierr);

  result = (flag == PETSC_TRUE);

  return 0;
}

//! Initializes the config parameter and flag database.
/*!
  Processes -config and -config_override command line options.
 */
PetscErrorCode init_config(MPI_Comm com, PetscMPIInt rank,
			   NCConfigVariable &config, NCConfigVariable &overrides,
                           bool process_options) {
  PetscErrorCode ierr;

  config.init("pism_config", com, rank);
  overrides.init("pism_overrides", com, rank);

  string alt_config = PISM_DefaultConfigFile,
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

  ierr = config.read(alt_config.c_str()); CHKERRQ(ierr);

  if (use_override_config) {
    ierr = overrides.read(override_config.c_str()); CHKERRQ(ierr);
    config.import_from(overrides);
    ierr = verbPrintf(2, com, "CONFIG OVERRIDES read from file '%s'.\n",
		      override_config.c_str()); CHKERRQ(ierr);
  }

  if (process_options) {
    ierr =  set_config_from_options(com, config); CHKERRQ(ierr);
  }

  config.print();

  return 0;
}

PetscErrorCode set_config_from_options(MPI_Comm /*com*/, NCConfigVariable &config) {
  PetscErrorCode ierr;
  bool flag;

  // Energy modeling
  ierr = config.flag_from_option("varc", "use_linear_in_temperature_heat_capacity");  CHKERRQ(ierr);
  ierr = config.flag_from_option("vark", "use_temperature_dependent_thermal_conductivity");  CHKERRQ(ierr);

  // see getBasalWaterPressure()
  ierr = config.flag_from_option("bmr_enhance",
                                 "bmr_enhance_basal_water_pressure"); CHKERRQ(ierr);
  // in units m a-1 :
  ierr = config.scalar_from_option("bmr_enhance_scale", "bmr_enhance_scale"); CHKERRQ(ierr);

  ierr = config.flag_from_option("bmr_in_cont", "include_bmr_in_continuity"); CHKERRQ(ierr);

  // if set, use old IceModel::temperatureStep(), and set enthalpy as though
  //   ice is cold
  ierr = config.flag_from_option("cold", "do_cold_ice_methods"); CHKERRQ(ierr);

  ierr = config.flag_from_option("diffuse_bwat", "do_diffuse_bwat"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("low_temp", "global_min_allowed_temp"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("max_low_temps", "max_low_temp_count"); CHKERRQ(ierr);

  // Sub-models
  ierr = config.flag_from_option("blatter", "do_blatter"); CHKERRQ(ierr);
  ierr = config.flag_from_option("age", "do_age"); CHKERRQ(ierr);
  ierr = config.flag_from_option("mass", "do_mass_conserve"); CHKERRQ(ierr);
  ierr = config.flag_from_option("energy", "do_energy"); CHKERRQ(ierr);
  ierr = config.flag_from_option("sia", "do_sia"); CHKERRQ(ierr);

  // Time-stepping
  ierr = config.scalar_from_option("adapt_ratio",
				   "adaptive_timestepping_ratio"); CHKERRQ(ierr);

  ierr = config.flag_from_option("count_steps", "count_time_steps"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("max_dt", "maximum_time_step_years"); CHKERRQ(ierr);

	// evaluates the adaptive timestep based on a CFL criterion with respect to the eigenCalving rate
  ierr = config.flag_from_option("cfl_eigencalving", "cfl_eigencalving"); CHKERRQ(ierr);


  // SIA
  ierr = config.scalar_from_option("bed_smoother_range", "bed_smoother_range"); CHKERRQ(ierr);

  ierr = config.keyword_from_option("gradient", "surface_gradient_method",
                                    "eta,haseloff,mahaffy"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("e", "enhancement_factor"); CHKERRQ(ierr);

  ierr = config.flag_from_option("e_age_coupling", "do_e_age_coupling"); CHKERRQ(ierr);

  // note "-gk" is used for specifying Goldsby-Kohlstedt ice
  //   this form allows a constant value of grain size to be input in mm
  ierr = config.scalar_from_option("gk", "constant_grain_size"); CHKERRQ(ierr);

  // SSA
  // Decide on the algorithm for solving the SSA
  ierr = config.keyword_from_option("ssa_method", "ssa_method", "fd,fem"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("ssa_eps",  "epsilon_ssafd"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_maxi", "max_iterations_ssafd"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_rtol", "ssafd_relative_convergence"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("e_ssa", "ssa_enhancement_factor"); CHKERRQ(ierr);

  ierr = config.flag_from_option("ssa_dirichlet_bc", "ssa_dirichlet_bc"); CHKERRQ(ierr);
  ierr = config.flag_from_option("cfbc", "calving_front_stress_boundary_condition"); CHKERRQ(ierr);
  ierr = config.flag_from_option("brutal_sliding", "scalebrutalSet"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("brutal_sliding_scale","sliding_scale_brutal"); CHKERRQ(ierr); 
 

  // Basal strength

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

  // threshold; at this velocity tau_c is basal shear stress
  ierr = config.scalar_from_option("pseudo_plastic_uthreshold", "pseudo_plastic_uthreshold"); CHKERRQ(ierr);

  // If set, makes the thickness affect the pore_pressure; near margin there is
  // a reduction in basal water pressure, a conceptual drainage mechanism
  ierr = config.flag_from_option("thk_eff", "thk_eff_basal_water_pressure"); CHKERRQ(ierr);
  // next two in  m  :
  ierr = config.scalar_from_option("thk_eff_H_high","thk_eff_H_high");  CHKERRQ(ierr);
  ierr = config.scalar_from_option("thk_eff_H_low","thk_eff_H_low");  CHKERRQ(ierr);
  // pure number :
  ierr = config.scalar_from_option("thk_eff_reduced","thk_eff_reduced");  CHKERRQ(ierr);

  // Ice shelves

  ierr = config.flag_from_option("part_grid", "part_grid"); CHKERRQ(ierr);

  ierr = config.flag_from_option("part_redist", "part_redist"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("nuBedrock", "nuBedrock"); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-nuBedrock", flag);  CHKERRQ(ierr);
  if (flag)  config.set_flag("nuBedrockSet", true);


  // Calving

  // whether or not to kill ice at locations that were ice-free at
  // bootstrapping
  ierr = config.flag_from_option("ocean_kill", "ocean_kill"); CHKERRQ(ierr);

  // whether or not to kill ice (zero thickness) if it is (or becomes) floating
  ierr = config.flag_from_option("float_kill", "floating_ice_killed"); CHKERRQ(ierr);

  ierr = config.flag_from_option("thickness_calving", "do_thickness_calving"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("calving_at_thickness", "calving_at_thickness"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("eigen_calving_K", "eigen_calving_K"); CHKERRQ(ierr);
  ierr = config.flag_from_option("eigen_calving", "do_eigen_calving"); CHKERRQ(ierr);

  ierr = config.flag_from_option("kill_icebergs", "kill_icebergs"); CHKERRQ(ierr);

  // Output
  ierr = config.flag_from_option("acab_cumulative", "compute_cumulative_acab"); CHKERRQ(ierr);
  ierr = config.flag_from_option("f3d", "force_full_diagnostics"); CHKERRQ(ierr);

  ierr = config.keyword_from_option("o_order", "output_variable_order",
                                    "xyz,yxz,zyx"); CHKERRQ(ierr);

  ierr = config.keyword_from_option("o_format", "output_format",
                                    "netcdf3,netcdf4_parallel"); CHKERRQ(ierr);

  ierr = config.scalar_from_option("summary_volarea_scale_factor_log10",
                                   "summary_volarea_scale_factor_log10"); CHKERRQ(ierr);

  ierr = config.flag_from_option("vpik", "verbose_pik_messages");  CHKERRQ(ierr);

  // Metadata
  ierr = config.string_from_option("title", "run_title"); CHKERRQ(ierr);
  ierr = config.string_from_option("institution", "institution"); CHKERRQ(ierr);

  // Skipping
  ierr = config.flag_from_option("skip", "do_skip"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("skip_max", "skip_max"); CHKERRQ(ierr);

  // Shortcuts

  if (getVerbosityLevel() > 2)  config.set_flag("verbose_pik_messages", true);

  // option "-pik" turns on a suite of PISMPIK effects (but not -eigen_calving)
  ierr = PISMOptionsIsSet("-pik", "enable suite of PISM-PIK mechanisms", flag); CHKERRQ(ierr);
  if (flag) {
    config.set_flag("calving_front_stress_boundary_condition", true);
    config.set_flag("part_grid", true);
    config.set_flag("part_redist", true);
    config.set_flag("kill_icebergs", true);
  }

  
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

  return 0;
}

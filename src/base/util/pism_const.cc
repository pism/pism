// Copyright (C) 2007--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petsc.h>
#include <petscfix.h>
#include "NCTool.hh"
#include "pism_const.hh"
#include <sstream>
#include <ctime>
#include <algorithm>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

//! \brief PISM verbosity level; determines how much gets printed to the
//! standard out.
static PetscInt verbosityLevel;

//! \brief Set the PISM verbosity level.
PetscErrorCode setVerbosityLevel(PetscInt level) {
  if ((level < 0) || (level > 5)) {
     SETERRQ(1,"verbosity level invalid");
  }
  verbosityLevel = level;
  return 0;  
}

//! \brief Get the verbosity level.
PetscInt getVerbosityLevel() {
  return verbosityLevel;
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
  PetscTruth     verbose, levelSet;
  
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


//! Print messages to standard out according to verbosity threshhold.
/*!
Verbosity level version of PetscPrintf.  We print according to whether 
(thresh <= verbosityLevel), in which case print, or (thresh > verbosityLevel)
in which case no print.

verbosityLevelFromOptions() actually reads the threshold.

The range 1 <= thresh <= 5  is enforced.

Use this method for messages and warnings which should
- go to stdout and
- appear only once (regardless of number of processors).

Should not be used for reporting fatal errors.
 */
PetscErrorCode verbPrintf(const int thresh, 
                          MPI_Comm comm,const char format[],...)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  size_t         len;
  char           *buffer,*sub1,*sub2;
  const char     *nformat;
  PetscReal      value;

  extern FILE *petsc_history;

  if ((thresh < 1) || (thresh > 5)) { SETERRQ(1,"invalid threshold in verbPrintf()"); }

  PetscFunctionBegin;
  if (!comm) comm = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  if (!rank && ((verbosityLevel >= thresh) || petsc_history) ) {
    va_list Argp;
    va_start(Argp,format);

    ierr = PetscStrstr(format,"%A",&sub1);CHKERRQ(ierr);
    if (sub1) {
      ierr = PetscStrstr(format,"%",&sub2);CHKERRQ(ierr);
      if (sub1 != sub2) SETERRQ(PETSC_ERR_ARG_WRONG,"%%A format must be first in format string");
      ierr    = PetscStrlen(format,&len);CHKERRQ(ierr);
      ierr    = PetscMalloc((len+16)*sizeof(char),&buffer);CHKERRQ(ierr);
      ierr    = PetscStrcpy(buffer,format);CHKERRQ(ierr);
      ierr    = PetscStrstr(buffer,"%",&sub2);CHKERRQ(ierr);
      sub2[0] = 0;
      value   = (double)va_arg(Argp,double);
      if (PetscAbsReal(value) < 1.e-12) {
        ierr    = PetscStrcat(buffer,"< 1.e-12");CHKERRQ(ierr);
      } else {
        ierr    = PetscStrcat(buffer,"%g");CHKERRQ(ierr);
        va_end(Argp);
        va_start(Argp,format);
      }
      ierr    = PetscStrcat(buffer,sub1+2);CHKERRQ(ierr);
      nformat = buffer;
    } else {
      nformat = format;
    }
    if (verbosityLevel >= thresh) {
      ierr = PetscVFPrintf(PETSC_STDOUT,nformat,Argp);CHKERRQ(ierr);
    }
    if (petsc_history) { // always print to history
      ierr = PetscVFPrintf(petsc_history,nformat,Argp);CHKERRQ(ierr);
    }
    va_end(Argp);
    if (sub1) {ierr = PetscFree(buffer);CHKERRQ(ierr);}
  }
  PetscFunctionReturn(0);
}


//! Prints rank (in process group PETSC_COMM_WORLD) to stderr.  Then attempts to end all processes.
/*!
Avoid using this if possible.  SETERRQ() should be used for all procedures that
return PetscErrorCode.  Generally needed for errors in constructors.

Printing the rank seems to be redundant because it will appear in brackets at start
of each call to PetscErrorPrintf().  But emphasis is useful, perhaps ...

Calls "MPI_Abort(PETSC_COMM_WORLD,3155)" to attempt to end all processes.
If this works main() will return value 3155.  The problem with PetscEnd() for
this purpose is that it is collective (presumably over PETSC_COMM_WORLD).
 */
void endPrintRank() {
  PetscMPIInt rank;
  if (!MPI_Comm_rank(PETSC_COMM_WORLD, &rank)) {
    PetscErrorPrintf("\n\n    rank %d process called endPrintRank()\n"
                         "    ending ...  \n\n",rank);
  } else {
    PetscErrorPrintf("\n\n    process with undeterminable rank called endPrintRank()\n"
                         "    ending ...  \n\n");
  }
  MPI_Abort(PETSC_COMM_WORLD,3155);
}

//! Print a warning telling the user that an option was ignored.
PetscErrorCode ignore_option(MPI_Comm com, const char name[]) {
  PetscErrorCode ierr;
  PetscTruth option_is_set;

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
  PetscTruth option_is_set;

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
  PetscTruth option_is_set;

  char tmp[1]; // dummy string
  ierr = PetscOptionsGetString(PETSC_NULL, name, tmp, 1, &option_is_set); CHKERRQ(ierr);

  if (option_is_set) {
    ierr = PetscPrintf(com, "PISM ERROR: command-line option '%s' is not allowed.\n",
		       name); CHKERRQ(ierr);
    PISMEnd();
  }

  return 0;
}

//! Returns true if \c str ends with \c suffix and false otherwise.
bool ends_with(string str, string suffix) {
  if (str.rfind(suffix) + suffix.size() == str.size())
    return true;

  return false;
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


//! Checks if a vector of doubles is strictly increasing.
bool is_increasing(const vector<double> &a) {
  int len = (int)a.size();
  for (PetscInt k = 0; k < len-1; k++) {
    if (a[k] >= a[k+1])  return false;
  }
  return true;
}

//! Creates a time-stamp used for the history NetCDF attribute.
string pism_timestamp() {
  time_t now;
  tm tm_now;
  char date_str[50];
  now = time(NULL);
  localtime_r(&now, &tm_now);
  // Format specifiers for strftime():
  //   %F = ISO date format,  %T = Full 24 hour time,  %Z = Time Zone name
  strftime(date_str, sizeof(date_str), "%F %T %Z", &tm_now);

  return string(date_str);
}

//! Creates a string with the user name, hostname and the time-stamp (for history strings).
string pism_username_prefix() {
  PetscErrorCode ierr;

  char username[50];
  ierr = PetscGetUserName(username, sizeof(username));
  if (ierr != 0)
    username[0] = '\0';
  char hostname[100];
  ierr = PetscGetHostName(hostname, sizeof(hostname));
  if (ierr != 0)
    hostname[0] = '\0';
  
  ostringstream message;
  message << username << "@" << hostname << " " << pism_timestamp() << ": ";

  return message.str();
}

//! \brief Uses argc and argv to create the string with current PISM
//! command-line arguments.
string pism_args_string() {
  PetscInt argc;
  char **argv;
  PetscGetArgs(&argc, &argv);

  string cmdstr, argument;
  for (int j = 0; j < argc; j++) {
    argument = argv[j];

    // enclose arguments containing spaces with double quotes:
    if (argument.find(" ") != string::npos) argument = "\"" + argument + "\"";

    cmdstr += string(" ") + argument;
  }
  cmdstr += "\n";

  return cmdstr;
}

//! Makes the process on rank \c rank_to_stop wait for the debugger to be attached.
/*
  Once the debugger is attached, giving it the "set var i = 1" (for gdb)
  command will resume the run.
 */
PetscErrorCode pism_wait_for_gdb(MPI_Comm com, PetscMPIInt rank_to_stop) {
  PetscErrorCode ierr;
  PetscMPIInt my_rank;

  ierr = MPI_Comm_rank(com, &my_rank); CHKERRQ(ierr);

  if (my_rank == rank_to_stop) {
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PISM PID %d on %s is waiting...\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
      sleep(5);
  }
  
  ierr = MPI_Barrier(com); CHKERRQ(ierr);

  return 0;
}

//! Initializes the config parameter and flag database.
/*!
  Processes -config and -config_override command line options.
 */
PetscErrorCode init_config(MPI_Comm com, PetscMPIInt rank,
			   NCConfigVariable &config, NCConfigVariable &overrides) {
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
  config.print();

  return 0;
}

/* 
   note on pass-by-reference for options: For the last argument "flag" to
   PetscOptionsXXXX(....,&flag), the flag always indicates whether the option
   has been set. Therefore "flag" is altered by this function call. For other
   arguments "value" to PetscOptionsXXXX(....,&value,&flag), the value of
   "value" is only set if the user specified the option. Therefore "flag"
   should always be given a local PetscTruth variable if we want to preserve
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
  PetscTruth opt_set = PETSC_FALSE;

  if (choices.empty()) {
    SETERRQ(1, "PISMOptionsList: empty choices argument");
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
  PetscTruth flag;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "",
			    result.c_str(), tmp,
			    TEMPORARY_STRING_LENGTH, &flag); CHKERRQ(ierr);

  is_set = (flag == PETSC_TRUE);

  if (is_set) {
    if (strlen(tmp) == 0) {
      if (allow_empty_arg)
        result = "";
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
  PetscTruth opt_set = PETSC_FALSE;

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
  PetscTruth flag;
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
  PetscTruth flag;
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
  PetscTruth flag;

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
  PetscTruth flag;

  ierr = PetscOptionsGetString(PETSC_NULL, option.c_str(), tmp, 1, &flag); CHKERRQ(ierr);

  result = (flag == PETSC_TRUE);

  return 0;
}

//! A version of PISMOptionsIsSet that prints a -help message.
PetscErrorCode PISMOptionsIsSet(string option, string text,
				bool &result) {
  PetscErrorCode ierr;
  char tmp[1];
  PetscTruth flag;

  ierr = PetscOptionsString(option.c_str(), text.c_str(), "",
			    "", tmp, 1, &flag); CHKERRQ(ierr);

  result = (flag == PETSC_TRUE);

  return 0;
}

//! \brief Adds a suffix to a filename.
/*!
 * Returns filename + separator + suffix + .nc if the original filename had the
 * .nc suffix, otherwise filename + separator. If the old filename had the form
 * "name + separator + more stuff + .nc", then removes the string after the
 * separator.
 */
string pism_filename_add_suffix(string filename, string separator, string suffix) {
  string basename = filename, result;

  // find where the separator begins:
  string::size_type j = basename.rfind(separator);
  if (j == string::npos) {
    j = basename.rfind(".nc");
  }

  // if the separator was not found, find the .nc suffix:
  if (j == string::npos) {
    j = basename.size();
  }

  // cut off everything starting from the separator (or the .nc suffix):
  basename.resize(static_cast<int>(j));

  result = basename + separator + suffix;

  if (ends_with(filename, ".nc"))
    result += ".nc";

  return result;
}

//! \brief Finalizes PETSc and MPI. Replaces PetscEnd().
/*!
 * The main reason for having this is pismebm, an executable running 2 MPI
 * communicators, only one of which runs PETSc. Using PetscEnd() in this case
 * leaves the process in the communicator \b not running PETsc hanging waiting
 * for a MPI_Finalize() call. (PetscFinalize() only calls MPI_Finalize() if
 * PetscInitialize() called MPI_Init().)
 */
void PISMEnd() {
  int flag;
  PetscFinalize();

  MPI_Finalized(&flag);
  if ( ! flag )
    MPI_Finalize();

  exit(0);
}

void PISMEndQuiet() {
  PetscOptionsSetValue("-options_left","no");
  PISMEnd();
}

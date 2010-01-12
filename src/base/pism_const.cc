// Copyright (C) 2007--2009 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "nc_util.hh"
#include "pism_const.hh"
#include <sstream>
#include <ctime>
#include <algorithm>

static PetscInt verbosityLevel;

PetscErrorCode setVerbosityLevel(PetscInt level) {
  if ((level < 0) || (level > 5)) {
     SETERRQ(1,"verbosity level invalid");
  }
  verbosityLevel = level;
  return 0;  
}

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
PetscErrorCode check_option(string name, PetscTruth &flag) {
  PetscErrorCode ierr;
  char tmp[1];

  ierr = PetscOptionsGetString(PETSC_NULL, name.c_str(), tmp, 1, &flag); CHKERRQ(ierr);

  return 0;
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
    PetscEnd();
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
    PetscEnd();
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
PetscErrorCode parse_range(MPI_Comm com, string str, double *a, double *delta, double *b) {
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

  // Convert each number from a string to double:
  for (int j = 0; j < 3; ++j) {
    double d;
    char *endptr;
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

//! Parses a time specification.
/*!
  If it is a MATLAB-style range, then calls parse_range and computes all the points.

  If it is a comma-separated list, converts to double (with error-checking).
 */
PetscErrorCode parse_times(MPI_Comm com, string str, vector<double> &result) {
  PetscErrorCode ierr;
  int N;

  if (str.find(':') != string::npos) { // it's a range specification
    
    double a, delta, b;
    ierr = parse_range(com, str, &a, &delta, &b);
    if (ierr != 0) return 1;

    if (a >= b) {
      ierr = PetscPrintf(com, "PISM ERROR: a >= b in the range specification %s.\n",
			 str.c_str()); CHKERRQ(ierr);
      return 1;
    }

    if (delta <= 0) {
      ierr = PetscPrintf(com, "PISM ERROR: delta <= 0 in the range specification %s.\n",
			 str.c_str()); CHKERRQ(ierr);
      return 1;
    }

    N = (int)floor((b - a)/delta) + 1; // number of points in the range
    result.resize(N);

    for (int j = 0; j < N; ++j)
      result[j] = a + delta*j;
    
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
	result.push_back(d);
    }
    sort(result.begin(), result.end());
  }

  return 0;
}


PetscErrorCode just_show_usage(
    MPI_Comm com, const char execname[], const char usage[]) {
  PetscErrorCode ierr;
  ierr = verbPrintf(1,com,
      "%s is a PISM (http://pism-docs.org) executable.", execname); CHKERRQ(ierr);
  ierr = verbPrintf(1,com,"\nOptions cheat-sheet:\n\n");
      CHKERRQ(ierr);
  ierr = verbPrintf(1,com,usage); CHKERRQ(ierr);
  ierr = verbPrintf(1,com,
      "\nTo run in parallel using N processors (most MPI cases):\n"
      "  mpiexec -n N %s ...\n"
      "\nFor more help on %s and PISM,\n"
      "  1. Download User's Manual for PISM: http://www.pism-docs.org/pdfs/manual0.2.pdf\n"
      "  2. Read browser for technical details on PISM: http://www.pism-docs.org/doxy/html/index.html\n"
      "  3. Search bugs and tasks at PISM source host: https://gna.org/projects/pism\n"
      "  4. Run with '-help | grep foo' to see PETSc options which relate to 'foo'.\n"
      "  5. Email for help:  help AT pism-docs.org\n", 
      execname, execname, execname);
      CHKERRQ(ierr);
  return 0;
}


//! Show provided usage message and quit.  (Consider using show_usage_check_req_opts() in preference to this one.)
PetscErrorCode show_usage_and_quit(
    MPI_Comm com, const char execname[], const char usage[]) {
  PetscErrorCode ierr;
  ierr = just_show_usage(com, execname, usage); CHKERRQ(ierr);
  PetscEnd();
  return 0;
}


//! In a single call a driver program can provide a usage string to the user and check if required options are given, and if not, end.
PetscErrorCode show_usage_check_req_opts(
    MPI_Comm com, const char execname[], vector<string> required_options,
    const char usage[]) {
  PetscErrorCode ierr;
  
  PetscTruth usageSet;
  ierr = check_option("-usage", usageSet); CHKERRQ(ierr);
  if (usageSet == PETSC_TRUE) {
    ierr = show_usage_and_quit(com, execname, usage); CHKERRQ(ierr);
  }

  // go through list of required options, and if not given, fail
  bool req_absent = false;
  for (size_t ii=0; ii < required_options.size(); ii++) {
    PetscTruth set;
    ierr = check_option(required_options[ii], set); CHKERRQ(ierr);
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
  PetscTruth helpSet;
  ierr = check_option("-help", helpSet); CHKERRQ(ierr);
  if (helpSet == PETSC_TRUE) {
    ierr = just_show_usage(com, execname, usage); CHKERRQ(ierr);
  }

  return 0;
}


//! Checks if an array of \c len doubles is strictly increasing.
bool is_increasing(int len, double *a) {
  for (PetscInt k = 0; k < len-1; k++) {
    if (a[k] >= a[k+1])  return false;
  }
  return true;
}

//! Creates a time-stamp used for the history NetCDF attribute.
string timestamp() {
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
string username_prefix() {
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
  message << username << "@" << hostname << " " << timestamp() << ": ";

  return message.str();
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

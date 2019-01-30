/* Copyright (C) 2016, 2017, 2018 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "pism_utilities.hh"

#include <sstream>              // istringstream, ostringstream
#include <cstdio>               // vsnprintf

#include <mpi.h>                // MPI_Get_library_version
#include <fftw3.h>              // fftw_version
#include <gsl/gsl_version.h>

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>             // nc_inq_libvers

#if (PISM_USE_PROJ4==1)
#include "pism/util/Proj.hh"    // pj_release
#endif

#ifdef PISM_USE_JANSSON
#include <jansson.h>            // JANSSON_VERSION
#endif

#include <petsctime.h>          // PetscTime

#include "error_handling.hh"

namespace pism {

const char *PISM_DefaultConfigFile = PISM_DEFAULT_CONFIG_FILE;

const char *PISM_Revision = PISM_REVISION;

//! Returns true if `str` ends with `suffix` and false otherwise.
bool ends_with(const std::string &str, const std::string &suffix) {
  if (suffix.size() > str.size()) {
    return false;
  }

  if (str.rfind(suffix) + suffix.size() == str.size()) {
    return true;
  }

  return false;
}

template <class T>
std::string join_impl(const T& input, const std::string& separator) {
  auto j = input.begin();
  std::string result = *j;
  ++j;
  while (j != input.end()) {
    result += separator + *j;
    ++j;
  }
  return result;
}

//! Concatenate `strings`, inserting `separator` between elements.
std::string join(const std::vector<std::string> &strings, const std::string &separator) {
  return join_impl(strings, separator);
}

std::string set_join(const std::set<std::string> &input, const std::string& separator) {
  return join_impl(input, separator);
}

//! Transform a `separator`-separated list (a string) into a vector of strings.
std::vector<std::string> split(const std::string &input, char separator) {
  std::istringstream input_list(input);
  std::string token;
  std::vector<std::string> result;

  while (getline(input_list, token, separator)) {
    if (not token.empty()) {
      result.push_back(token);
    }
  }
  return result;
}

//! Transform a `separator`-separated list (a string) into a set of strings.
std::set<std::string> set_split(const std::string &input, char separator) {
  std::istringstream input_list(input);
  std::string token;
  std::set<std::string> result;

  while (getline(input_list, token, separator)) {
    if (not token.empty()) {
      result.insert(token);
    }
  }
  return result;
}

//! Checks if a vector of doubles is strictly increasing.
bool is_increasing(const std::vector<double> &a) {
  int len = (int)a.size();
  for (int k = 0; k < len-1; k++) {
    if (a[k] >= a[k+1]) {
      return false;
    }
  }
  return true;
}

bool member(const std::string &string, const std::set<std::string> &set) {
  return (set.find(string) != set.end());
}

void GlobalReduce(MPI_Comm comm, double *local, double *result, int count, MPI_Op op) {
  int err = MPI_Allreduce(local, result, count, MPI_DOUBLE, op, comm);
  PISM_C_CHK(err, 0, "MPI_Allreduce");
}

void GlobalMin(MPI_Comm comm, double *local, double *result, int count) {
  GlobalReduce(comm, local, result, count, MPI_MIN);
}

void GlobalMax(MPI_Comm comm, double *local, double *result, int count) {
  GlobalReduce(comm, local, result, count, MPI_MAX);
}

void GlobalSum(MPI_Comm comm, double *local, double *result, int count) {
  GlobalReduce(comm, local, result, count, MPI_SUM);
}

unsigned int GlobalSum(MPI_Comm comm, unsigned int input) {
  unsigned int result;
  int err = MPI_Allreduce(&input, &result, 1, MPI_UNSIGNED, MPI_SUM, comm);
  PISM_C_CHK(err, 0, "MPI_Allreduce");
  return result;
}

int GlobalSum(MPI_Comm comm, int input) {
  int result;
  int err = MPI_Allreduce(&input, &result, 1, MPI_INT, MPI_SUM, comm);
  PISM_C_CHK(err, 0, "MPI_Allreduce");
  return result;
}

double GlobalMin(MPI_Comm comm, double local) {
  double result;
  GlobalMin(comm, &local, &result, 1);
  return result;
}

double GlobalMax(MPI_Comm comm, double local) {
  double result;
  GlobalMax(comm, &local, &result, 1);
  return result;
}

double GlobalSum(MPI_Comm comm, double local) {
  double result;
  GlobalSum(comm, &local, &result, 1);
  return result;
}

static const int TEMPORARY_STRING_LENGTH = 32768;

std::string version() {
  char buffer[TEMPORARY_STRING_LENGTH];
  std::string result;

  snprintf(buffer, sizeof(buffer), "PISM (%s)\n", PISM_Revision);
  result += buffer;

  snprintf(buffer, sizeof(buffer), "CMake %s.\n", PISM_CMAKE_VERSION);
  result += buffer;

  PetscGetVersion(buffer, TEMPORARY_STRING_LENGTH);
  result += buffer;
  result += "\n";

#ifdef PISM_PETSC_CONFIGURE_FLAGS
  snprintf(buffer, sizeof(buffer), "PETSc configure: %s\n",
           PISM_PETSC_CONFIGURE_FLAGS);
  result += buffer;
#endif

  // OpenMPI added MPI_Get_library_version in version 1.7 (relatively recently).
#ifdef OPEN_MPI
  snprintf(buffer, TEMPORARY_STRING_LENGTH, "OpenMPI %d.%d.%d\n",
           OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION, OMPI_RELEASE_VERSION);
#else
  // Assume that other MPI libraries implement this part of the MPI-3 standard...
  int string_length = TEMPORARY_STRING_LENGTH;
  MPI_Get_library_version(buffer, &string_length);
#endif
  result += buffer;

  snprintf(buffer, sizeof(buffer), "NetCDF %s.\n", nc_inq_libvers());
  result += buffer;

  snprintf(buffer, sizeof(buffer), "FFTW %s.\n", fftw_version);
  result += buffer;

  snprintf(buffer, sizeof(buffer), "GSL %s.\n", GSL_VERSION);
  result += buffer;

#if (PISM_USE_PROJ4==1)
  snprintf(buffer, sizeof(buffer), "PROJ.4 %s.\n", pj_release);
  result += buffer;
#endif

#ifdef PISM_USE_JANSSON
  snprintf(buffer, sizeof(buffer), "Jansson %s.\n", JANSSON_VERSION);
  result += buffer;
#endif

#ifdef PISM_SWIG_VERSION
  snprintf(buffer, sizeof(buffer), "SWIG %s.\n", PISM_SWIG_VERSION);
  result += buffer;
#endif

#ifdef PISM_PETSC4PY_VERSION
  snprintf(buffer, sizeof(buffer), "petsc4py %s.\n", PISM_PETSC4PY_VERSION);
  result += buffer;
#endif

  return result;
}


//! Return time since the beginning of the run, in hours.
double wall_clock_hours(MPI_Comm com, double start_time) {
  int rank = 0;
  double result = 0.0;

  MPI_Comm_rank(com, &rank);

  ParallelSection rank0(com);
  try {
    if (rank == 0) {
      result = (get_time() - start_time) / 3600.0;
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  MPI_Bcast(&result, 1, MPI_DOUBLE, 0, com);

  return result;
}

//! Creates a time-stamp used for the history NetCDF attribute.
std::string timestamp(MPI_Comm com) {
  time_t now;
  tm tm_now;
  char date_str[50];
  now = time(NULL);
  localtime_r(&now, &tm_now);
  // Format specifiers for strftime():
  //   %F = ISO date format,  %T = Full 24 hour time,  %Z = Time Zone name
  strftime(date_str, sizeof(date_str), "%F %T %Z", &tm_now);

  MPI_Bcast(date_str, 50, MPI_CHAR, 0, com);

  return std::string(date_str);
}

//! Creates a string with the user name, hostname and the time-stamp (for history strings).
std::string username_prefix(MPI_Comm com) {
  PetscErrorCode ierr;

  char username[50];
  ierr = PetscGetUserName(username, sizeof(username));
  PISM_CHK(ierr, "PetscGetUserName");
  if (ierr != 0) {
    username[0] = '\0';
  }
  char hostname[100];
  ierr = PetscGetHostName(hostname, sizeof(hostname));
  PISM_CHK(ierr, "PetscGetHostName");
  if (ierr != 0) {
    hostname[0] = '\0';
  }

  std::ostringstream message;
  message << username << "@" << hostname << " " << timestamp(com) << ": ";

  std::string result = message.str();
  unsigned int length = result.size();
  MPI_Bcast(&length, 1, MPI_UNSIGNED, 0, com);

  result.resize(length);
  MPI_Bcast(&result[0], length, MPI_CHAR, 0, com);

  return result;
}

//! \brief Uses argc and argv to create the string with current PISM
//! command-line arguments.
std::string args_string() {
  int argc;
  char **argv;
  PetscErrorCode ierr = PetscGetArgs(&argc, &argv);
  PISM_CHK(ierr, "PetscGetArgs");

  std::string cmdstr, argument;
  for (int j = 0; j < argc; j++) {
    argument = argv[j];

    // enclose arguments containing spaces with double quotes:
    if (argument.find(" ") != std::string::npos) {
      argument = "\"" + argument + "\"";
    }

    cmdstr += std::string(" ") + argument;
  }
  cmdstr += "\n";

  return cmdstr;
}

//! \brief Adds a suffix to a filename.
/*!
 * Returns filename + separator + suffix + .nc if the original filename had the
 * .nc suffix, otherwise filename + separator. If the old filename had the form
 * "name + separator + more stuff + .nc", then removes the string after the
 * separator.
 */
std::string filename_add_suffix(const std::string &filename,
                                     const std::string &separator,
                                     const std::string &suffix) {
  std::string basename = filename, result;

  // find where the separator begins:
  std::string::size_type j = basename.rfind(separator);
  if (j == std::string::npos) {
    j = basename.rfind(".nc");
  }

  // if the separator was not found, find the .nc suffix:
  if (j == std::string::npos) {
    j = basename.size();
  }

  // cut off everything starting from the separator (or the .nc suffix):
  basename.resize(static_cast<int>(j));

  result = basename + separator + suffix;

  if (ends_with(filename, ".nc")) {
    result += ".nc";
  }

  return result;
}

double get_time() {
  PetscLogDouble result;
  PetscErrorCode ierr = PetscTime(&result); PISM_CHK(ierr, "PetscTime");
  return result;
}

std::string printf(const char *format, ...) {
  std::string result(1024, ' ');
  va_list arglist;
  size_t length;

  va_start(arglist, format);
  if((length = vsnprintf(&result[0], result.size(), format, arglist)) > result.size()) {
    result.reserve(length);
    vsnprintf(&result[0], result.size(), format, arglist);
  }
  va_end(arglist);
  return result.substr(0, length);
}

} // end of namespace pism

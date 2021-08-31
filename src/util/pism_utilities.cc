/* Copyright (C) 2016, 2017, 2018, 2019, 2020, 2021 PISM Authors
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

#include <cstdarg>              // va_list, va_start(), va_end()
#include <sstream>              // istringstream, ostringstream
#include <cstdio>               // vsnprintf
#include <cassert>              // assert

#include <mpi.h>                // MPI_Get_library_version
#include <fftw3.h>              // fftw_version
#include <gsl/gsl_version.h>    // GSL_VERSION

#include "pism/pism_config.hh"  // Pism_USE_XXX, version info

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>             // nc_inq_libvers

#if (Pism_USE_PROJ==1)
#include "pism/util/Proj.hh"    // pj_release
#endif

#if (Pism_USE_JANSSON==1)
#include <jansson.h>            // JANSSON_VERSION
#endif

#include <petsctime.h>          // PetscTime

#include <cstdlib>              // strtol(), strtod()

#include "error_handling.hh"

namespace pism {

std::string string_strip(const std::string &input) {
  if (input.empty()) {
    return "";
  }

  std::string tmp = input;

  // strip leading spaces
  tmp.erase(0, tmp.find_first_not_of(" \t"));

  // strip trailing spaces
  tmp.substr(tmp.find_last_not_of(" \t"));

  return tmp;
}

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
      result.emplace_back(token);
    }
  }
  return result;
}

//! Transform a `separator`-separated list (a string) into a set of strings.
std::set<std::string> set_split(const std::string &input, char separator) {
  std::set<std::string> result;
  for (const auto &token : split(input, separator)) {
    result.insert(token);
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
  assert(local != result);
  int err = MPI_Allreduce(local, result, count, MPI_DOUBLE, op, comm);
  PISM_C_CHK(err, 0, "MPI_Allreduce");
}

void GlobalReduce(MPI_Comm comm, int *local, int *result, int count, MPI_Op op) {
  assert(local != result);
  int err = MPI_Allreduce(local, result, count, MPI_INT, op, comm);
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

void GlobalSum(MPI_Comm comm, int *local, int *result, int count) {
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

  result += pism::printf("PISM (%s)\n", pism::revision);
  result += pism::printf("CMake %s.\n", pism::cmake_version);

  PetscGetVersion(buffer, TEMPORARY_STRING_LENGTH);
  result += buffer;
  result += "\n";
  result += pism::printf("PETSc configure: %s\n", pism::petsc_configure_flags);

  // OpenMPI added MPI_Get_library_version in version 1.7 (relatively recently).
#ifdef OPEN_MPI
  result += pism::printf("OpenMPI %d.%d.%d\n",
                         OMPI_MAJOR_VERSION,
                         OMPI_MINOR_VERSION,
                         OMPI_RELEASE_VERSION);
#else
  // Assume that other MPI libraries implement this part of the MPI-3 standard...
  int string_length = TEMPORARY_STRING_LENGTH;
  MPI_Get_library_version(buffer, &string_length);
  result += buffer;
#endif

  result += pism::printf("NetCDF %s.\n", nc_inq_libvers());
  result += pism::printf("FFTW %s.\n", fftw_version);
  result += pism::printf("GSL %s.\n", GSL_VERSION);

#if (Pism_USE_PROJ==1)
  result += pism::printf("PROJ %s.\n", pj_release);
#endif

#if (Pism_USE_JANSSON==1)
  result += pism::printf("Jansson %s.\n", JANSSON_VERSION);
#endif

#if (Pism_BUILD_PYTHON_BINDINGS==1)
  result += pism::printf("SWIG %s.\n", pism::swig_version);
  result += pism::printf("petsc4py %s.\n", pism::petsc4py_version);
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

/*!
 * Validate a format string. In this application a format string should contain `%` exactly
 * once, followed by `s` (i.e. `%s`).
 *
 * Throws RuntimeError if the provided string is invalid.
 */
void validate_format_string(const std::string &format) {
  if (format.find("%s") == std::string::npos) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "format string %s does not contain %%s",
                                  format.c_str());
  }

  if (format.find("%") != format.rfind("%")) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "format string %s contains more than one %%",
                                  format.c_str());
  }
}

double vector_min(const std::vector<double> &input) {
  double my_min = input[0];
  for (auto x : input) {
    my_min = std::min(x, my_min);
  }
  return my_min;
}

double vector_max(const std::vector<double> &input) {
  double my_max = input[0];
  for (auto x : input) {
    my_max = std::max(x, my_max);
  }
  return my_max;
}


/*!
 * Fletcher's checksum
 *
 * See https://en.wikipedia.org/wiki/Fletcher%27s_checksum#Optimizations
 */
uint64_t fletcher64(const uint32_t *data, size_t length) {
  // Accumulating a sum of block_size unsigned 32-bit integers in an unsigned 64-bit
  // integer will not lead to an overflow.
  //
  // This constant is found by solving n * (n + 1) / 2 * (2^32 - 1) < (2^64 - 1).
  static const size_t block_size = 92681;

  uint64_t c0 = 0, c1 = 0;
  while (length != 0) {
    size_t block = std::min(block_size, length);

    for (size_t i = 0; i < block; ++i) {
      c0 = c0 + *data++;
      c1 = c1 + c0;
    }

    c0 = c0 % UINT32_MAX;
    c1 = c1 % UINT32_MAX;

    length = length > block_size ? length - block_size : 0;
  }
  return (c1 << 32 | c0);
}

/*!
 * Compute water column pressure vertically-averaged over the height of an ice cliff at a
 * margin.
 */
double average_water_column_pressure(double ice_thickness, double bed,
                                     double floatation_level, double rho_ice,
                                     double rho_water, double g) {

  double
    ice_bottom = std::max(bed, floatation_level - rho_ice / rho_water * ice_thickness),
    water_column_height = std::max(floatation_level - ice_bottom, 0.0);

  if (ice_thickness > 0.0) {
    return 0.5 * rho_water * g * pow(water_column_height, 2.0) / ice_thickness;
  }
  return 0.0;
}

double parse_number(const std::string &input) {
  char *endptr = NULL;
  double result = strtod(input.c_str(), &endptr);
  if (*endptr != '\0') {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Can't parse %s (expected a floating point number)",
                                  input.c_str());
  }
  return result;
}

long int parse_integer(const std::string &input) {
  char *endptr = NULL;
  long int result = strtol(input.c_str(), &endptr, 10);
  if (*endptr != '\0') {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Can't parse %s (expected an integer)",
                                  input.c_str());
  }
  return result;
}

} // end of namespace pism

/* Copyright (C) 2016 PISM Authors
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

#include <mpi.h>
#include <netcdf.h>
#include <fftw3.h>
#include <gsl/gsl_version.h>
#if (PISM_USE_PROJ4==1)
#include "base/util/Proj.hh"
#endif

#if (PISM_USE_JANSSON==1)
#include <jansson.h>
#endif

#include "pism_utilities.hh"
#include "pism_const.hh"
#include "error_handling.hh"

#include <sstream>

namespace pism {

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
  typename T::const_iterator j = input.begin();
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
    result.push_back(token);
  }
  return result;
}

//! Transform a `separator`-separated list (a string) into a vector of strings.
std::set<std::string> set_split(const std::string &input, char separator) {
  std::istringstream input_list(input);
  std::string token;
  std::set<std::string> result;

  while (getline(input_list, token, separator)) {
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

bool set_contains(const std::set<std::string> &S, const std::string &name) {
  return (S.find(name) != S.end());
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

  snprintf(buffer, sizeof(buffer), "PETSc configure: %s\n",
           PISM_PETSC_CONFIGURE_FLAGS);
  result += buffer;

  int string_length = TEMPORARY_STRING_LENGTH;
  MPI_Get_library_version(buffer, &string_length);
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

#if (PISM_USE_JANSSON==1)
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

} // end of namespace pism

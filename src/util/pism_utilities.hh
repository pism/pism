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

#ifndef PISM_UTILITIES_H
#define PISM_UTILITIES_H

#include <cstdint>              // uint32_t

#include <algorithm>            // std::min, std::max
#include <set>
#include <string>
#include <vector>

#include <mpi.h>

namespace pism {

/*!
 * Compute vertically-integrated water column pressure.
 */
double average_water_column_pressure(double ice_thickness, double bed,
                                     double floatation_level,
                                     double rho_ice, double rho_water, double g);

// Utilities that do not expose PETSc's or PISM's API.

#ifndef __GNUC__
#  define  __attribute__(x)  /* nothing */
#endif

double get_time();
std::string timestamp(MPI_Comm com);
std::string username_prefix(MPI_Comm com);
std::string args_string();
std::string filename_add_suffix(const std::string &filename,
                                const std::string &separator,
                                const std::string &suffix);

double wall_clock_hours(MPI_Comm com, double start_time);


// array
bool is_increasing(const std::vector<double> &a);

// string
bool ends_with(const std::string &str, const std::string &suffix);

// remove leading and trailing whitespace
std::string string_strip(const std::string &input);

std::string join(const std::vector<std::string> &strings, const std::string &separator);

std::vector<std::string> split(const std::string &input, char separator);

std::set<std::string> set_split(const std::string &input, char separator);

std::string set_join(const std::set<std::string> &input, const std::string& separator);

// set
bool member(const std::string &string, const std::set<std::string> &set);

/*! Helper template function for computing set unions.
 * Ensures that elements of a take precedence. For example, if
 *
 * a = {{1, 2}, {3, 4}}
 * b = {{1, 4}, {5, 6}}
 *
 * combine(a, b) will use the pair {1, 2} from a, not {1, 4} from b.
 *
 * This behavior relies on the fact that std::map::insert({a, b}) is a no-op if a key equivalent to
 * a is already present.
 *
 * This is similar to a set union, but it is not symmetric. (I would expect set_union(a, b) to be
 * the same as set_union(b, a)).
 */
template<typename T>
T combine(const T &a, const T&b) {
  T result = a;
  for (const auto &element : b) {
    result.insert(element);
  }
  return result;
}

template<typename T>
inline T clip(T x, T a, T b) {
  return std::min(std::max(a, x), b);
}

double vector_min(const std::vector<double> &input);

double vector_max(const std::vector<double> &input);

// parallel
void GlobalReduce(MPI_Comm comm, double *local, double *result, int count, MPI_Op op);

void GlobalReduce(MPI_Comm comm, int *local, int *result, int count, MPI_Op op);

void GlobalMin(MPI_Comm comm, double *local, double *result, int count);

void GlobalMax(MPI_Comm comm, double *local, double *result, int count);

void GlobalMax(MPI_Comm comm, int *local, int *result, int count);

void GlobalSum(MPI_Comm comm, double *local, double *result, int count);

void GlobalSum(MPI_Comm comm, int *local, int *result, int count);

double GlobalMin(MPI_Comm comm, double local);

double GlobalMax(MPI_Comm comm, double local);

double GlobalSum(MPI_Comm comm, double local);

unsigned int GlobalSum(MPI_Comm comm, unsigned int input);

int GlobalSum(MPI_Comm comm, int input);

std::string version();

std::string printf(const char *format, ...) __attribute__((format(printf, 1, 2)));

void validate_format_string(const std::string &format);

uint64_t fletcher64(const uint32_t *data, size_t len);

void print_checksum(MPI_Comm com,
                    const std::vector<double> &data,
                    const char *label);

double parse_number(const std::string &input);

long int parse_integer(const std::string &input);

} // end of namespace pism


#endif /* PISM_UTILITIES_H */

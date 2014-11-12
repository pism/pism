// Copyright (C) 2009--2014 Constantine Khroulev
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

#include <algorithm>
#include <gsl/gsl_math.h>
#include <assert.h>

#include "Timeseries.hh"
#include "pism_const.hh"
#include "IceGrid.hh"
#include "PIO.hh"
#include "PISMTime.hh"
#include "PISMConfig.hh"

#include "error_handling.hh"

namespace pism {

Timeseries::Timeseries(IceGrid *g, const std::string &name, const std::string &dimension_name)
  : m_unit_system(g->get_unit_system()),
    dimension(dimension_name, dimension_name, m_unit_system),
    var(name, dimension_name, m_unit_system),
    bounds(dimension_name + "_bounds", dimension_name, m_unit_system)
{
  private_constructor(g->com, name, dimension_name);
}

Timeseries::Timeseries(MPI_Comm c, const UnitSystem &unit_system,
                       const std::string & name, const std::string & dimension_name)
  : m_unit_system(unit_system),
    dimension(dimension_name, dimension_name, m_unit_system),
    var(name, dimension_name, m_unit_system),
    bounds(dimension_name + "_bounds", dimension_name, m_unit_system)
{
  private_constructor(c, name, dimension_name);
}

void Timeseries::private_constructor(MPI_Comm c, const std::string &name, const std::string &dimension_name) {
  com = c;
  dimension.set_string("bounds", dimension_name + "_bounds");

  short_name = name;
  use_bounds = true;
}

//! Ensure that time bounds have the same units as the dimension.
void Timeseries::set_bounds_units() {
  bounds.set_units(dimension.get_string("units"));
  bounds.set_glaciological_units(dimension.get_string("glaciological_units"));
}


//! Read timeseries data from a NetCDF file `filename`.
PetscErrorCode Timeseries::read(const PIO &nc, Time *time_manager) {
  PetscErrorCode ierr;

  bool exists, found_by_standard_name;
  std::vector<std::string> dims;
  std::string time_name, standard_name = var.get_string("standard_name"),
    name_found;

  nc.inq_var(short_name, standard_name,
             exists, name_found, found_by_standard_name);

  if (!exists) {
    throw RuntimeError::formatted("Can't find '%s' ('%s') in '%s'.\n",
                                  short_name.c_str(), standard_name.c_str(),
                                  nc.inq_filename().c_str());
  }

  dims = nc.inq_vardims(name_found);

  if (dims.size() != 1) {
    throw RuntimeError::formatted("Variable '%s' in '%s' depends on %d dimensions,\n"
                                  "but a time-series variable can only depend on 1 dimension.",
                                  short_name.c_str(),
                                  nc.inq_filename().c_str(),
                                  dims.size());
  }

  time_name = dims[0];

  NCTimeseries tmp_dim = dimension;
  tmp_dim.set_name(time_name);

  nc.read_timeseries(tmp_dim, time_manager, time);
  bool is_increasing = true;
  for (unsigned int j = 1; j < time.size(); ++j) {
    if (time[j] - time[j-1] < 1e-16) {
      is_increasing = false;
      break;
    }
  }
  if (!is_increasing) {
    throw RuntimeError::formatted("dimension '%s' has to be strictly increasing (read from '%s').",
                                  tmp_dim.get_name().c_str(), nc.inq_filename().c_str());
  }

  std::string time_bounds_name = nc.get_att_text(time_name, "bounds");

  if (!time_bounds_name.empty()) {
    use_bounds = true;

    set_bounds_units();
    NCTimeBounds tmp_bounds = bounds;
    tmp_bounds.set_name(time_bounds_name);

    ierr = tmp_bounds.set_units(tmp_dim.get_string("units")); CHKERRQ(ierr);

    nc.read_time_bounds(tmp_bounds, time_manager, time_bounds);
  } else {
    use_bounds = false;
  }

  nc.read_timeseries(var, time_manager, values);

  if (time.size() != values.size()) {
    throw RuntimeError::formatted("variables %s and %s in %s have different numbers of values.",
                                  dimension.get_name().c_str(),
                                  var.get_name().c_str(),
                                  nc.inq_filename().c_str());
  }

  ierr = report_range(); CHKERRQ(ierr);

  return 0;
}

//! \brief Report the range of a time-series stored in `values`.
PetscErrorCode Timeseries::report_range() {
  PetscErrorCode ierr;
  double min, max;

  // min_element and max_element return iterators; "*" is used to get
  // the value corresponding to this iterator
  min = *std::min_element(values.begin(), values.end());
  max = *std::max_element(values.begin(), values.end());

  assert(UnitConverter::are_convertible(var.get_units(), var.get_glaciological_units()) == true);
  UnitConverter c(var.get_units(), var.get_glaciological_units());
  min = c(min);
  max = c(max);

  std::string spacer(var.get_name().size(), ' ');

  ierr = verbPrintf(2, com,
                    "  FOUND  %s / %-60s\n"
                    "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
                    var.get_name().c_str(),
                    var.get_string("long_name").c_str(), spacer.c_str(), min, max,
                    var.get_string("glaciological_units").c_str()); CHKERRQ(ierr);

  return 0;
}

//! Write timeseries data to a NetCDF file `filename`.
PetscErrorCode Timeseries::write(const PIO &nc) {
  // write the dimensional variable; this call should go first
  nc.write_timeseries(dimension, 0, time);
  nc.write_timeseries(var, 0, values);

  if (use_bounds) {
    set_bounds_units();
    nc.write_time_bounds(bounds, 0, time_bounds);
  }

  return 0;
}

/** Scale all values stored in this instance by `scaling_factor`.
 *
 * This is used to convert mass balance offsets from [m s-1] to [kg m-2 s-1].
 *
 * @param[in] scaling_factor multiplicative scaling factor
 */
void Timeseries::scale(double scaling_factor) {
  for (unsigned int i = 0; i < values.size(); ++i) {
    values[i] *= scaling_factor;
  }
}

//! Get a value of timeseries at time `t`.
/*! Returns the first value or the last value if t is out of range on the left
  and right, respectively.

  Uses time bounds if present (interpreting data as piecewise-constant) and
  uses linear interpolation otherwise.
 */
double Timeseries::operator()(double t) {

  // piecewise-constant case:
  if (use_bounds) {
    std::vector<double>::iterator j;

    j = lower_bound(time_bounds.begin(), time_bounds.end(), t); // binary search

    if (j == time_bounds.end()) {
      return values.back(); // out of range (on the right)
    }

    int i = (int)(j - time_bounds.begin());

    if (i == 0) {
      return values[0];         // out of range (on the left)
    }

    if (i % 2 == 0) {
      throw RuntimeError::formatted("time bounds array in %s does not represent continguous time intervals.\n"
                                    "(PISM was trying to compute %s at time %3.3f seconds.)",
                                    bounds.get_name().c_str(), short_name.c_str(), t);
    }

    return values[(i-1)/2];
  }

  // piecewise-linear case:
  std::vector<double>::iterator end = time.end(), j;
  
  j = lower_bound(time.begin(), end, t); // binary search

  if (j == end) {
    return values.back(); // out of range (on the right)
  }

  int i = (int)(j - time.begin());

  if (i == 0) {
    return values[0];   // out of range (on the left)
  }

  double dt = time[i] - time[i - 1];
  double dv = values[i] - values[i - 1];
  
  return values[i - 1] + (t - time[i - 1]) / dt * dv;
}

//! Get a value of timeseries by index.
/*!
  Stops if the index is out of range.
 */
double Timeseries::operator[](unsigned int j) const {

#if (PISM_DEBUG==1)
  if (j >= values.size()) {
    throw RuntimeError::formatted("Timeseries %s: operator[]: invalid argument: size=%d, index=%d",
                                  var.get_name().c_str(), values.size(), j);
  }
#endif

  return values[j];
}

//! \brief Compute an average of a time-series over interval (t,t+dt) using
//! trapezoidal rule with N sub-intervals.
double Timeseries::average(double t, double dt, unsigned int N) {
  std::vector<double> V(N+1);
 
  for (unsigned int i = 0; i < N+1; ++i) {
    double t_i = t + (dt / N) * i;
    V[i] = (*this)(t_i);
  }

  double sum = 0;
  for (unsigned int i = 0; i < N; ++i) {
    sum += V[i] + V[i+1];
  }

  return sum / (2.0*N);
}

//! Append a pair (t,v) to the timeseries.
PetscErrorCode Timeseries::append(double v, double a, double b) {
  time.push_back(b);
  values.push_back(v);
  time_bounds.push_back(a);
  time_bounds.push_back(b);
  return 0;
}

NCTimeseries& Timeseries::get_metadata() {
  return var;
}

NCTimeseries& Timeseries::get_dimension_metadata() {
  return dimension;
}

//! Returns the length of the time-series stored.
/*!
  This length is changed by read() and append().
 */
int Timeseries::length() {
  return (int)values.size();
}

//----- DiagnosticTimeseries

DiagnosticTimeseries::DiagnosticTimeseries(IceGrid *g, const std::string &name, const std::string &dimension_name)
  : Timeseries(g, name, dimension_name) {

  buffer_size = (size_t)g->config.get("timeseries_buffer_size");
  start = 0;
  rate_of_change = false;
  dimension.set_string("calendar", g->time->calendar());
  dimension.set_string("long_name", "time");
  dimension.set_string("axis", "T");
}

//! Destructor; makes sure that everything is written to a file.
DiagnosticTimeseries::~DiagnosticTimeseries() {
  flush();
}

//! Adds the (t,v) pair to the interpolation buffer.
/*! The interpolation buffer holds 2 values only (for linear interpolation).
 *
 * If this DiagnosticTimeseries object is reporting a "rate of change",
 * append() has to be called with the "cumulative" quantity as the V argument.
 */
PetscErrorCode DiagnosticTimeseries::append(double V, double /*a*/, double b) {

  if (rate_of_change && v.empty()) {
    v_previous = V;
  }

  // append to the interpolation buffer
  t.push_back(b);
  v.push_back(V);

  if (t.size() == 3) {
    t.pop_front();
    v.pop_front();
  }

  return 0;
}

//! \brief Use linear interpolation to find the value of a scalar diagnostic
//! quantity at time `T` and store the obtained pair (T, value).
PetscErrorCode DiagnosticTimeseries::interp(double a, double b) {
  PetscErrorCode ierr;

  if (t.empty()) {
    throw RuntimeError("DiagnosticTimeseries::interp(...): interpolation buffer is empty");
  }

  if (t.size() == 1) {
    time.push_back(b);
    values.push_back(GSL_NAN);
    time_bounds.push_back(a);
    time_bounds.push_back(b);
    return 0;
  }

  if ((b < t[0]) || (b > t[1])) {
    throw RuntimeError::formatted("DiagnosticTimeseries::interp(...): requested time %f is not within the last time-step!",
                                  b);
  }

  double 
    // compute the "cumulative" quantity using linear interpolation
    v_current = v[0] + (b - t[0]) / (t[1] - t[0]) * (v[1] - v[0]),
    // the value to report
    value = v_current;

  if (rate_of_change) {
    // use backward-in-time finite difference to compute the rate of change:
    value = (v_current - v_previous) / (b - a);

    // remember the value of the "cumulative" quantity for differencing during
    // the next call:
    v_previous = v_current;
  }

  // use the right endpoint as the 'time' record (the midpoint is also an option)
  time.push_back(b);
  values.push_back(value);

  // save the time bounds
  time_bounds.push_back(a);
  time_bounds.push_back(b);

  if (time.size() == buffer_size) {
    ierr = flush(); CHKERRQ(ierr);
  }

  return 0;
}
PetscErrorCode DiagnosticTimeseries::init(const std::string &filename) {
  PIO nc(com, "netcdf3", m_unit_system); // OK to use netcdf3
  unsigned int len = 0;

  // Get the number of records in the file (for appending):
  bool file_exists = nc.check_if_exists(com, filename);

  if (file_exists == true) {
    nc.open(filename, PISM_READONLY);
    len = nc.inq_dimlen(dimension.get_name());
    if (len > 0) {
      // read the last value and initialize v_previous and v[0]
      std::vector<double> tmp;
      bool var_exists = nc.inq_var(short_name);

      if (var_exists) {
        nc.get_1d_var(short_name, len - 1, 1, tmp);
        // NOTE: this is WRONG if rate_of_change == true!
        v.push_back(tmp[0]);
        v_previous = tmp[0];
      }
    }
    nc.close();
  }

  output_filename = filename;
  start = len;

  return 0;
}


  //! Writes data to a file.
PetscErrorCode DiagnosticTimeseries::flush() {
  PIO nc(com, "netcdf3", m_unit_system); // OK to use netcdf3
  unsigned int len = 0;

  // return cleanly if this DiagnosticTimeseries object was created but never
  // used:
  if (output_filename.empty()) {
    return 0;
  }

  if (time.empty()) {
    return 0;
  }

  nc.open(output_filename, PISM_READWRITE);
  len = nc.inq_dimlen(dimension.get_name());

  if (len > 0) {
    double last_time;
    nc.inq_dim_limits(dimension.get_dimension_name(),
                      NULL, &last_time);
    if (last_time < time.front()) {
      start = len;
    }
  }

  if (len == (unsigned int)start) {
    nc.write_timeseries(dimension, start, time);  

    set_bounds_units();
    nc.write_time_bounds(bounds, start, time_bounds);
  }
  nc.write_timeseries(var, start, values);

  start += time.size();

  time.clear();
  values.clear();
  time_bounds.clear();

  nc.close();

  return 0;
}

void DiagnosticTimeseries::reset() {
  time.clear();
  values.clear();
  time_bounds.clear();
  start = 0;
  t.clear();
  v.clear();
}


} // end of namespace pism

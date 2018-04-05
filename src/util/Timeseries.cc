// Copyright (C) 2009--2018 Constantine Khroulev
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

#include <algorithm>            // min_element and max_element

#include "Timeseries.hh"
#include "pism_utilities.hh"
#include "IceGrid.hh"
#include "pism/util/io/PIO.hh"
#include "Time.hh"
#include "ConfigInterface.hh"

#include "error_handling.hh"
#include "io/io_helpers.hh"
#include "pism/util/Logger.hh"

namespace pism {

Timeseries::Timeseries(const IceGrid &g, const std::string &name, const std::string &dimension_name)
  : m_dimension(dimension_name, dimension_name, g.ctx()->unit_system()),
    m_variable(name, dimension_name, g.ctx()->unit_system()),
    m_bounds(dimension_name + "_bounds", dimension_name, g.ctx()->unit_system())
{
  private_constructor(g.com, dimension_name);
}

Timeseries::Timeseries(MPI_Comm c, units::System::Ptr unit_system,
                       const std::string & name, const std::string & dimension_name)
  : m_dimension(dimension_name, dimension_name, unit_system),
    m_variable(name, dimension_name, unit_system),
    m_bounds(dimension_name + "_bounds", dimension_name, unit_system)
{
  private_constructor(c, dimension_name);
}

void Timeseries::private_constructor(MPI_Comm c, const std::string &dimension_name) {
  m_com = c;
  m_dimension.set_string("bounds", dimension_name + "_bounds");

  m_use_bounds = true;
}

//! Ensure that time bounds have the same units as the dimension.
void Timeseries::set_bounds_units() {
  m_bounds.set_string("units", m_dimension.get_string("units"));
  m_bounds.set_string("glaciological_units", m_dimension.get_string("glaciological_units"));
}

bool Timeseries::get_use_bounds() const {
  return m_use_bounds;
}

void Timeseries::set_use_bounds(bool flag) {
  m_use_bounds = flag;
}


//! Read timeseries data from a NetCDF file `filename`.
void Timeseries::read(const PIO &nc, const Time &time_manager, const Logger &log) {

  bool exists, found_by_standard_name;
  std::vector<std::string> dims;
  std::string time_name, standard_name = m_variable.get_string("standard_name"),
    name_found;

  nc.inq_var(name(), standard_name,
             exists, name_found, found_by_standard_name);

  if (!exists) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' ('%s') in '%s'.\n",
                                  name().c_str(), standard_name.c_str(),
                                  nc.inq_filename().c_str());
  }

  dims = nc.inq_vardims(name_found);

  if (dims.size() != 1) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Variable '%s' in '%s' depends on %d dimensions,\n"
                                  "but a time-series variable can only depend on 1 dimension.",
                                  name().c_str(),
                                  nc.inq_filename().c_str(),
                                  (int)dims.size());
  }

  time_name = dims[0];

  TimeseriesMetadata tmp_dim = m_dimension;
  tmp_dim.set_name(time_name);

  io::read_timeseries(nc, tmp_dim, time_manager, log, m_time);
  bool is_increasing = true;
  for (unsigned int j = 1; j < m_time.size(); ++j) {
    if (m_time[j] - m_time[j-1] < 1e-16) {
      is_increasing = false;
      break;
    }
  }
  if (!is_increasing) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "dimension '%s' has to be strictly increasing (read from '%s').",
                                  tmp_dim.get_name().c_str(), nc.inq_filename().c_str());
  }

  std::string time_bounds_name = nc.get_att_text(time_name, "bounds");

  if (!time_bounds_name.empty()) {
    m_use_bounds = true;

    set_bounds_units();
    TimeBoundsMetadata tmp_bounds = m_bounds;
    tmp_bounds.set_name(time_bounds_name);

    tmp_bounds.set_string("units", tmp_dim.get_string("units"));

    io::read_time_bounds(nc, tmp_bounds, time_manager, log, m_time_bounds);
  } else {
    m_use_bounds = false;
  }

  io::read_timeseries(nc, m_variable, time_manager, log, m_values);

  if (m_time.size() != m_values.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "variables %s and %s in %s have different numbers of values.",
                                  m_dimension.get_name().c_str(),
                                  m_variable.get_name().c_str(),
                                  nc.inq_filename().c_str());
  }

  report_range(log);
}

//! \brief Report the range of a time-series stored in `values`.
void Timeseries::report_range(const Logger &log) {
  double min, max;

  // min_element and max_element return iterators; "*" is used to get
  // the value corresponding to this iterator
  min = *std::min_element(m_values.begin(), m_values.end());
  max = *std::max_element(m_values.begin(), m_values.end());

  units::Converter c(m_variable.unit_system(),
                  m_variable.get_string("units"),
                  m_variable.get_string("glaciological_units"));
  min = c(min);
  max = c(max);

  std::string spacer(m_variable.get_name().size(), ' ');

  log.message(2,
             "  FOUND  %s / %-60s\n"
             "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
             m_variable.get_name().c_str(),
             m_variable.get_string("long_name").c_str(), spacer.c_str(), min, max,
             m_variable.get_string("glaciological_units").c_str());
}

//! Write timeseries data to a NetCDF file `filename`.
void Timeseries::write(const PIO &nc) const {
  // write the dimensional variable; this call should go first
  io::write_timeseries(nc, m_dimension, 0, m_time);
  io::write_timeseries(nc, m_variable, 0, m_values);

  if (m_use_bounds) {
    io::write_time_bounds(nc, m_bounds, 0, m_time_bounds);
  }
}

//! Clear storage.
void Timeseries::reset() {
  m_time.clear();
  m_values.clear();
  m_time_bounds.clear();
}

/** Scale all values stored in this instance by `scaling_factor`.
 *
 * This is used to convert mass balance offsets from [m s-1] to [kg m-2 s-1].
 *
 * @param[in] scaling_factor multiplicative scaling factor
 */
void Timeseries::scale(double scaling_factor) {
  for (unsigned int i = 0; i < m_values.size(); ++i) {
    m_values[i] *= scaling_factor;
  }
}

//! Get a value of timeseries at time `t`.
/*! Returns the first value or the last value if t is out of range on the left
  and right, respectively.

  Uses time bounds if present (interpreting data as piecewise-constant) and
  uses linear interpolation otherwise.
 */
double Timeseries::operator()(double t) const {

  // piecewise-constant case:
  if (m_use_bounds) {
    auto j = lower_bound(m_time_bounds.begin(), m_time_bounds.end(), t); // binary search

    if (j == m_time_bounds.end()) {
      return m_values.back(); // out of range (on the right)
    }

    int i = (int)(j - m_time_bounds.begin());

    if (i == 0) {
      return m_values[0];         // out of range (on the left)
    }

    if (i % 2 == 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "time bounds array in %s does not represent continguous time intervals.\n"
                                    "(PISM was trying to compute %s at time %3.3f seconds.)",
                                    m_bounds.get_name().c_str(), name().c_str(), t);
    }

    return m_values[(i-1)/2];
  }

  // piecewise-linear case:
  auto end = m_time.end();

  auto j = lower_bound(m_time.begin(), end, t); // binary search

  if (j == end) {
    return m_values.back(); // out of range (on the right)
  }

  int i = (int)(j - m_time.begin());

  if (i == 0) {
    return m_values[0];   // out of range (on the left)
  }

  double dt = m_time[i] - m_time[i - 1];
  double dv = m_values[i] - m_values[i - 1];

  return m_values[i - 1] + (t - m_time[i - 1]) / dt * dv;
}

//! Get a value of timeseries by index.
/*!
  Stops if the index is out of range.
 */
double Timeseries::operator[](unsigned int j) const {

#if (PISM_DEBUG==1)
  if (j >= m_values.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Timeseries %s: operator[]: invalid argument: size=%d, index=%d",
                                  m_variable.get_name().c_str(), (int)m_values.size(), j);
  }
#endif

  return m_values[j];
}

//! \brief Compute an average of a time-series over interval (t,t+dt) using
//! trapezoidal rule with N sub-intervals.
double Timeseries::average(double t, double dt, unsigned int N) const {
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
void Timeseries::append(double v, double t0, double t1) {
  m_time.push_back(t1);
  m_values.push_back(v);
  m_time_bounds.push_back(t0);
  m_time_bounds.push_back(t1);
}

std::string Timeseries::name() const {
  return m_variable.get_name();
}

TimeseriesMetadata& Timeseries::variable() {
  return m_variable;
}

TimeseriesMetadata& Timeseries::dimension() {
  return m_dimension;
}

TimeBoundsMetadata& Timeseries::bounds() {
  return m_bounds;
}

const TimeseriesMetadata& Timeseries::variable() const {
  return m_variable;
}

const TimeseriesMetadata& Timeseries::dimension() const {
  return m_dimension;
}

const TimeBoundsMetadata& Timeseries::bounds() const {
  return m_bounds;
}

const std::vector<double> & Timeseries::times() const {
  return m_time;
}

const std::vector<double> & Timeseries::time_bounds() const {
  return m_time_bounds;
}

const std::vector<double> & Timeseries::values() const {
  return m_values;
}

} // end of namespace pism

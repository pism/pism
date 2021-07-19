// Copyright (C) 2009--2021 Constantine Khroulev
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
#include <gsl/gsl_interp.h>     // gsl_interp_bsearch

#include "Timeseries.hh"
#include "pism_utilities.hh"
#include "IceGrid.hh"
#include "pism/util/io/File.hh"

#include "error_handling.hh"
#include "io/io_helpers.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Context.hh"

namespace pism {

Timeseries::Timeseries(const IceGrid &grid,
                       const std::string &name)
  : m_com(grid.ctx()->com()),
    m_use_bounds(true),
    m_unit_system(grid.ctx()->unit_system()),
    m_variable(name, m_unit_system)
{
  // empty
}

Timeseries::Timeseries(MPI_Comm com,
                       units::System::Ptr unit_system,
                       const std::string &name)
  : m_com(com),
    m_use_bounds(true),
    m_unit_system(unit_system),
    m_variable(name, m_unit_system)
{
  // empty
}

//! Read timeseries data from a NetCDF file `filename`.
void Timeseries::read(const File &file, const std::string &time_units,
                      const Logger &log) {

  std::string standard_name = m_variable.get_string("standard_name");

  auto var = file.find_variable(name(), standard_name);

  if (not var.exists) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' ('%s') in '%s'.\n",
                                  name().c_str(), standard_name.c_str(),
                                  file.filename().c_str());
  }

  auto dims = file.dimensions(var.name);

  if (dims.size() != 1) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Variable '%s' in '%s' depends on %d dimensions,\n"
                                  "but a time-series variable can only depend on 1 dimension.",
                                  name().c_str(),
                                  file.filename().c_str(),
                                  (int)dims.size());
  }

  auto time_name = dims[0];

  VariableMetadata time_dimension(time_name, m_unit_system);
  time_dimension.set_string("units", time_units);

  io::read_timeseries(file, time_dimension, log, m_time);

  if (not is_increasing(m_time)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "dimension '%s' has to be strictly increasing (read from '%s').",
                                  time_dimension.get_name().c_str(), file.filename().c_str());
  }

  std::string time_bounds_name = file.read_text_attribute(time_name, "bounds");

  if (not time_bounds_name.empty()) {
    m_use_bounds = true;

    VariableMetadata time_bounds(time_bounds_name, m_unit_system);
    time_bounds.set_string("units", time_dimension.get_string("units"));
    time_bounds.set_string("glaciological_units", time_dimension.get_string("glaciological_units"));

    io::read_time_bounds(file, time_bounds, log, m_time_bounds);

    // Time bounds override the time dimension read from the input file.
    //
    // Note: we use the right end point of each interval.
    for (unsigned int k = 0; k < m_time.size(); ++k) {
      m_time[k] = m_time_bounds[2 * k + 1];
    }
  } else {
    m_use_bounds = false;
  }

  io::read_timeseries(file, m_variable, log, m_values);

  if (m_time.size() != m_values.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "variables %s and %s in %s have different numbers of values.",
                                  time_name.c_str(),
                                  m_variable.get_name().c_str(),
                                  file.filename().c_str());
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

  units::Converter c(m_unit_system,
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

//! Get a value of timeseries at time `t`.
/*! Returns the first value or the last value if t is out of range on the left
  and right, respectively.

  Uses time bounds if present (interpreting data as piecewise-constant) and
  uses linear interpolation otherwise.
 */
double Timeseries::operator()(double t) const {

  if (m_use_bounds) {
    // piecewise-constant case

    size_t k = 0;
    if (t < m_time[0]) {
      k = 0;
    } else if (t >= m_time.back()) {
      k = m_time.size() - 1;
    } else {
      k = gsl_interp_bsearch(m_time.data(), t, 0, m_time.size()) + 1;
    }

    return m_values[k];
  } else {
    // piecewise-linear case

    // extrapolation on the left
    if (t < m_time[0]) {
      return m_values[0];
    }

    size_t k = gsl_interp_bsearch(m_time.data(), t, 0, m_time.size());

    // extrapolation on the right
    if (k + 1 >= m_time.size()) {
      return m_values.back();
    }

    double alpha = (t - m_time[k]) / (m_time[k + 1] - m_time[k]);

    return m_values[k] * (1.0 - alpha) + m_values[k + 1] * alpha;
  }
}

//! Get a value of timeseries by index.
double Timeseries::operator[](unsigned int j) const {

#if (Pism_DEBUG==1)
  if (j >= m_values.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Timeseries %s: operator[]: invalid argument: size=%d, index=%d",
                                  name().c_str(), (int)m_values.size(), j);
  }
#endif

  return m_values[j];
}

//! @brief Compute an average of a time-series over interval (t,t+dt) using trapezoidal
//! rule with N sub-intervals.
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

std::string Timeseries::name() const {
  return m_variable.get_name();
}

VariableMetadata& Timeseries::variable() {
  return m_variable;
}

} // end of namespace pism

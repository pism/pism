/* Copyright (C) 2015, 2016, 2017, 2019, 2020, 2021, 2022, 2023 PISM Authors
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
#include <cmath>

#include "pism/util/Diagnostic.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Logger.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Context.hh"

namespace pism {

Diagnostic::Diagnostic(std::shared_ptr<const Grid> grid)
  : m_grid(grid),
    m_sys(grid->ctx()->unit_system()),
    m_config(grid->ctx()->config()),
    m_fill_value(m_config->get_number("output.fill_value")) {
  // empty
}

void Diagnostic::update(double dt) {
  this->update_impl(dt);
}

void Diagnostic::update_impl(double dt) {
  (void) dt;
  // empty
}

void Diagnostic::reset() {
  this->reset_impl();
}

void Diagnostic::reset_impl() {
  // empty
}

/*!
 * Convert from external (output) units to internal units.
 */
double Diagnostic::to_internal(double x) const {
  std::string
    out = m_vars.at(0).get_string("output_units");
  std::string
    in  = m_vars.at(0).get_string("units");
  return convert(m_sys, x, out, in);
}

/*!
 * Convert from internal to external (output) units.
 */
double Diagnostic::to_external(double x) const {
  std::string
    out = m_vars.at(0).get_string("output_units"),
    in  = m_vars.at(0).get_string("units");
  return convert(m_sys, x, in, out);
}

//! Get the number of NetCDF variables corresponding to a diagnostic quantity.
unsigned int Diagnostic::n_variables() const {
  return m_vars.size();
}

void Diagnostic::init(const File &input, unsigned int time) {
  this->init_impl(input, time);
}

void Diagnostic::define_state(const File &output) const {
  this->define_state_impl(output);
}

void Diagnostic::write_state(const File &output) const {
  this->write_state_impl(output);
}

void Diagnostic::init_impl(const File &input, unsigned int time) {
  (void) input;
  (void) time;
  // empty
}

void Diagnostic::define_state_impl(const File &output) const {
  (void) output;
  // empty
}

void Diagnostic::write_state_impl(const File &output) const {
  (void) output;
  // empty
}

//! Get a metadata object corresponding to variable number N.
SpatialVariableMetadata& Diagnostic::metadata(unsigned int N) {
  if (N >= m_vars.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "variable metadata index %d is out of bounds",
                                  N);
  }

  return m_vars[N];
}

void Diagnostic::define(const File &file, io::Type default_type) const {
  this->define_impl(file, default_type);
}

//! Define NetCDF variables corresponding to a diagnostic quantity.
void Diagnostic::define_impl(const File &file, io::Type default_type) const {
  for (const auto &v : m_vars) {
    io::define_spatial_variable(v, *m_grid, file, default_type);
  }
}

//! \brief A method for setting common variable attributes.
void Diagnostic::set_attrs(const std::string &long_name,
                           const std::string &standard_name,
                           const std::string &units,
                           const std::string &output_units,
                           unsigned int N) {
  if (N >= m_vars.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "N (%d) >= m_dof (%d)", N,
                                  static_cast<int>(m_vars.size()));
  }

  m_vars[N]["long_name"] = long_name;

  m_vars[N]["standard_name"] = standard_name;

  if (units == output_units) {
    // No unit conversion needed to write data, so we assume that we don't need to check
    // if units are valid. This is needed to be able to set units like "Pa s^(1/3)", which
    // are correct (ice hardness with the Glen exponent n=3) but are not supported by
    // UDUNITS.
    //
    // Also note that this automatically sets output_units.
    m_vars[N].set_units_without_validation(units);
  } else {
    m_vars[N]["units"] = units;

    if (not (m_config->get_flag("output.use_MKS") or output_units.empty())) {
      m_vars[N]["output_units"] = output_units;
    }
  }
}

array::Array::Ptr Diagnostic::compute() const {
  // use the name of the first variable
  std::vector<std::string> names;
  for (const auto &v : m_vars) {
    names.push_back(v.get_name());
  }
  std::string all_names = join(names, ",");

  m_grid->ctx()->log()->message(3, "-  Computing %s...\n", all_names.c_str());
  array::Array::Ptr result = this->compute_impl();
  m_grid->ctx()->log()->message(3, "-  Done computing %s.\n", all_names.c_str());

  return result;
}

TSDiagnostic::TSDiagnostic(std::shared_ptr<const Grid> grid, const std::string &name)
  : m_grid(grid),
    m_config(grid->ctx()->config()),
    m_sys(grid->ctx()->unit_system()),
    m_time_name(grid->ctx()->config()->get_string("time.dimension_name")),
    m_variable(name, m_sys),
    m_dimension(m_time_name, m_sys),
    m_time_bounds(m_time_name + "_bounds", m_sys) {

  m_current_time = 0;
  m_start        = 0;

  m_buffer_size = static_cast<size_t>(m_config->get_number("output.timeseries.buffer_size"));

  m_variable["ancillary_variables"] = name + "_aux";

  m_dimension["calendar"] = m_grid->ctx()->time()->calendar();
  m_dimension["long_name"] = "time";
  m_dimension["axis"] = "T";
  m_dimension["units"] = m_grid->ctx()->time()->units_string();
}

TSDiagnostic::~TSDiagnostic() {
  flush();
}

void TSDiagnostic::set_units(const std::string &units,
                             const std::string &output_units) {
  m_variable["units"] = units;

  if (not m_config->get_flag("output.use_MKS")) {
    m_variable["output_units"] = output_units;
  }
}

TSSnapshotDiagnostic::TSSnapshotDiagnostic(std::shared_ptr<const Grid> grid, const std::string &name)
  : TSDiagnostic(grid, name) {
  // empty
}

TSRateDiagnostic::TSRateDiagnostic(std::shared_ptr<const Grid> grid, const std::string &name)
  : TSDiagnostic(grid, name), m_accumulator(0.0), m_v_previous(0.0), m_v_previous_set(false) {
  // empty
}

TSFluxDiagnostic::TSFluxDiagnostic(std::shared_ptr<const Grid> grid, const std::string &name)
  : TSRateDiagnostic(grid, name) {
  // empty
}

void TSSnapshotDiagnostic::evaluate(double t0, double t1, double v) {

  // skip times before the beginning of this time step
  while (m_current_time < m_requested_times->size() and (*m_requested_times)[m_current_time] < t0) {
    m_current_time += 1;
  }

  while (m_current_time < m_requested_times->size() and (*m_requested_times)[m_current_time] <= t1) {
    const unsigned int k = m_current_time;
    m_current_time += 1;

    // skip the first time: it defines the beginning of a reporting interval
    if (k == 0) {
      continue;
    }

    const double t_s = (*m_requested_times)[k - 1];
    const double t_e = (*m_requested_times)[k];

    // store computed data in the buffer
    {
      m_time.push_back(t_e);
      m_values.push_back(v);
      m_bounds.push_back(t_s);
      m_bounds.push_back(t_e);
    }
  }
}

void TSRateDiagnostic::evaluate(double t0, double t1, double change) {
  static const double epsilon = 1e-4; // seconds
  assert(t1 > t0);

  // skip times before and including the beginning of this time step
  while (m_current_time < m_requested_times->size() and (*m_requested_times)[m_current_time] <= t0) {
    m_current_time += 1;
  }

  // number of requested times in this time step
  unsigned int N = 0;

  // loop through requested times that are within this time step
  while (m_current_time < m_requested_times->size() and (*m_requested_times)[m_current_time] <= t1) {
    const unsigned int k = m_current_time;
    m_current_time += 1;

    N += 1;

    // skip the first time: it defines the beginning of a reporting interval
    if (k == 0) {
      continue;
    }

    const double t_s = (*m_requested_times)[k - 1];
    const double t_e = (*m_requested_times)[k];

    double rate = 0.0;
    if (N == 1) {
      // it is the right end-point of the first reporting interval in this time step: count the
      // contribution from the last time step plus the one from the beginning of this time step
      const double
        total_change  = m_accumulator + change * (t_e - t0) / (t1 - t0);
      const double dt = t_e - t_s;

      rate = total_change / dt;

    } else {
      // this reporting interval is completely contained within the time step, so the rate of change
      // does not depend on its length
      rate = change / (t1 - t0);
    }

    // store computed data in the buffer
    {
      m_time.push_back(t_e);
      m_values.push_back(rate);
      m_bounds.push_back(t_s);
      m_bounds.push_back(t_e);
    }

    m_accumulator = 0.0;
  }

  if (N == 0) {
    // if this time step contained no requested times we need to add the whole change to the
    // accumulator
    m_accumulator += change;
  } else {
    // if this time step contained some requested times we need to add the change since the last one
    // to the accumulator
    const double dt = t1 - (*m_requested_times)[m_current_time - 1];
    if (dt > epsilon) {
      m_accumulator += change * (dt / (t1 - t0));
    }
  }
}

void TSDiagnostic::update(double t0, double t1) {
  this->update_impl(t0, t1);
}

void TSSnapshotDiagnostic::update_impl(double t0, double t1) {
  static const double epsilon = 1e-4; // seconds

  if (fabs(t1 - t0) < epsilon) {
    return;
  }

  assert(t1 > t0);

  evaluate(t0, t1, this->compute());
}

void TSRateDiagnostic::update_impl(double t0, double t1) {
  const double v = this->compute();

  if (m_v_previous_set) {
    assert(t1 > t0);
    evaluate(t0, t1, v - m_v_previous);
  }

  m_v_previous = v;
  m_v_previous_set = true;
}

void TSFluxDiagnostic::update_impl(double t0, double t1) {
  static const double epsilon = 1e-4; // seconds

  if (fabs(t1 - t0) < epsilon) {
    return;
  }

  assert(t1 > t0);

  evaluate(t0, t1, this->compute());
}

void TSDiagnostic::define(const File &file) const {
  auto time_name = m_config->get_string("time.dimension_name");
  io::define_timeseries(m_variable, time_name, file, io::PISM_DOUBLE);
  io::define_time_bounds(m_time_bounds, time_name, "nv", file, io::PISM_DOUBLE);
}

void TSDiagnostic::flush() {

  if (m_time.empty()) {
    return;
  }

  std::string dimension_name = m_dimension.get_name();

  File file(m_grid->com, m_output_filename, io::PISM_NETCDF3, io::PISM_READWRITE); // OK to use netcdf3

  unsigned int len = file.dimension_length(dimension_name);

  if (len > 0) {
    double last_time = vector_max(file.read_dimension(dimension_name));

    if (last_time < m_time.front()) {
      m_start = len;
    }
  }

  auto time_name = m_config->get_string("time.dimension_name");

  if (len == m_start) {
    if (not file.find_variable(m_dimension.get_name())) {
      io::define_timeseries(m_dimension, time_name, file, io::PISM_DOUBLE);
    }
    io::write_timeseries(file, m_dimension, m_start, m_time);

    if (not file.find_variable(m_time_bounds.get_name())) {
      io::define_time_bounds(m_time_bounds, time_name, "nv", file, io::PISM_DOUBLE);
    }
    io::write_time_bounds(file, m_time_bounds, m_start, m_bounds);
  }

  if (not file.find_variable(m_variable.get_name())) {
    io::define_timeseries(m_variable, time_name, file, io::PISM_DOUBLE);
  }
  io::write_timeseries(file, m_variable, m_start, m_values);

  m_start += m_time.size();

  {
    m_time.clear();
    m_bounds.clear();
    m_values.clear();
  }
}

void TSDiagnostic::init(const File &output_file,
                        std::shared_ptr<std::vector<double>> requested_times) {
  m_output_filename = output_file.filename();

  m_requested_times = std::move(requested_times);

  // Get the number of records in the file (for appending):
  m_start = output_file.dimension_length(m_dimension.get_name());
}

const VariableMetadata &TSDiagnostic::metadata() const {
  return m_variable;
}

} // end of namespace pism

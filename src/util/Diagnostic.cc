/* Copyright (C) 2015, 2016, 2017, 2019, 2020, 2021, 2022, 2023, 2024, 2025 PISM Authors
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
#include <memory>

#include "io/IO_Flags.hh"
#include "io/OutputWriter.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Logger.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Context.hh"
#include "pism/util/projection.hh"

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
    out = m_vars.at(0)["output_units"],
    in  = m_vars.at(0)["units"];
  return convert(m_sys, x, out, in);
}

/*!
 * Convert from internal to external (output) units.
 */
double Diagnostic::to_external(double x) const {
  std::string
    out = m_vars.at(0)["output_units"],
    in  = m_vars.at(0)["units"];
  return convert(m_sys, x, in, out);
}

//! Get the number of NetCDF variables corresponding to a diagnostic quantity.
unsigned int Diagnostic::n_variables() const {
  return m_vars.size();
}

void Diagnostic::init(const File &input, unsigned int time) {
  this->init_impl(input, time);
}

void Diagnostic::define_state(const OutputFile &output) const {
  this->define_state_impl(output);
}

void Diagnostic::write_state(const OutputFile &output) const {
  this->write_state_impl(output);
}

void Diagnostic::init_impl(const File &input, unsigned int time) {
  (void) input;
  (void) time;
  // empty
}

void Diagnostic::define_state_impl(const OutputFile &output) const {
  (void) output;
  // empty
}

void Diagnostic::write_state_impl(const OutputFile &output) const {
  (void) output;
  // empty
}

//! Get a metadata object corresponding to variable number N.
VariableMetadata& Diagnostic::metadata(unsigned int N) {
  if (N >= m_vars.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "variable metadata index %d is out of bounds",
                                  N);
  }

  return m_vars[N];
}

const grid::DistributedGridInfo &Diagnostic::grid_info() const {
  return m_grid->info();
}

std::shared_ptr<array::Array> Diagnostic::compute() const {
  std::vector<std::string> names;
  for (const auto &v : m_vars) {
    names.push_back(v.get_name());
  }
  std::string all_names = join(names, ",");

  m_grid->ctx()->log()->message(3, "-  Computing %s...\n", all_names.c_str());
  auto result = this->compute_impl();
  m_grid->ctx()->log()->message(3, "-  Done computing %s.\n", all_names.c_str());

  return result;
}

TSDiagnostic::TSDiagnostic(std::shared_ptr<const Grid> grid, const std::string &name)
  : m_grid(grid),
    m_config(grid->ctx()->config()),
    m_sys(grid->ctx()->unit_system()),
    m_variable(name, m_sys),
    m_time_dimension(grid->ctx()->time()->metadata()) {

  m_current_time = 0;
  m_start        = 0;

  m_buffer_size = static_cast<size_t>(m_config->get_number("output.timeseries.buffer_size"));

  m_variable["ancillary_variables"] = name + "_aux";
  m_variable.set_time_dependent(true);
}

TSDiagnostic::~TSDiagnostic() {
  flush();
}

void TSDiagnostic::set_units(const std::string &units,
                             const std::string &output_units) {
  m_variable.units(units);

  if (not m_config->get_flag("output.use_MKS")) {
    m_variable.output_units(output_units);
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

void TSDiagnostic::flush() {

  if (m_time.empty()) {
    return;
  }

  auto &file = *m_output_file;
  auto len = file.time_dimension_length();

  if (len > 0) {
    // Note: does not perform unit conversion of the time read from the file. This should
    // be OK because this file was written by PISM.
    double last_time = file.last_time_value();

    if (last_time < m_time.front()) {
      m_start = len;
    }
  }

  if (len == m_start) {
    bool with_bounds = true;
    io::define_time(file, m_time_dimension, with_bounds);

    // write requested times
    file.write_array(m_time_dimension, { m_start }, { (unsigned int)m_time.size() }, m_time);
    // write time bounds
    auto bounds = m_time_dimension.get_name() + "_bounds";
    file.write_array({ bounds, m_sys }, { m_start, 0 }, { (unsigned int)m_bounds.size() / 2, 2 },
                     m_bounds);
  }

  file.define_variable(m_variable);
  // write values of a diagnostic
  file.write_timeseries_variable(m_variable.get_name(), { m_start },
                                 { (unsigned int)m_values.size() }, m_values);

  file.sync();

  m_start += m_time.size();

  {
    m_time.clear();
    m_bounds.clear();
    m_values.clear();
  }
}

void TSDiagnostic::init(std::shared_ptr<OutputFile> output_file,
                        std::shared_ptr<std::vector<double>> requested_times) {
  m_output_file = output_file;

  m_requested_times = std::move(requested_times);

  // Get the number of records in the file (for appending):
  m_start = output_file->time_dimension_length();
}

const VariableMetadata &TSDiagnostic::metadata() const {
  return m_variable;
}

} // end of namespace pism

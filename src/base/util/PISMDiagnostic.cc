/* Copyright (C) 2015, 2016, 2017 PISM Authors
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

#include "PISMDiagnostic.hh"
#include "error_handling.hh"
#include "io/io_helpers.hh"
#include "base/util/Logger.hh"
#include "base/util/pism_utilities.hh"

namespace pism {

Diagnostic::Ptr Diagnostic::wrap(const IceModelVec2S &input) {
  return Ptr(new DiagWithDedicatedStorage<IceModelVec2S>(input));
}

Diagnostic::Ptr Diagnostic::wrap(const IceModelVec2V &input) {
  return Ptr(new DiagWithDedicatedStorage<IceModelVec2V>(input));
}

Diagnostic::Diagnostic(IceGrid::ConstPtr g)
  : m_grid(g), m_sys(g->ctx()->unit_system()), m_config(g->ctx()->config()) {
  m_dof = 1;
  m_fill_value = m_config->get_double("output.fill_value");
}

Diagnostic::~Diagnostic() {
  // empty
}

//! \brief Update a cumulative quantity needed to compute a rate of change.
//! So far we there is only one such quantity: the rate of change of the ice
//! thickness.
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

//! Get the number of NetCDF variables corresponding to a diagnostic quantity.
unsigned int Diagnostic::n_variables() const {
  return m_dof;
}

void Diagnostic::init(const PIO &input, unsigned int time) {
  this->init_impl(input, time);
}

void Diagnostic::define_state(const PIO &output) const {
  this->define_state_impl(output);
}

void Diagnostic::write_state(const PIO &output) const {
  this->write_state_impl(output);
}

void Diagnostic::init_impl(const PIO &input, unsigned int time) {
  (void) input;
  (void) time;
  // empty
}

void Diagnostic::define_state_impl(const PIO &output) const {
  (void) output;
  // empty
}

void Diagnostic::write_state_impl(const PIO &output) const {
  (void) output;
  // empty
}

//! Get a metadata object corresponding to variable number N.
SpatialVariableMetadata& Diagnostic::metadata(unsigned int N) {
  if (N >= m_dof) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "variable metadata index %d is out of bounds",
                                  N);
  }

  return m_vars[N];
}

void Diagnostic::define(const PIO &file, IO_Type default_type) const {
  this->define_impl(file, default_type);
}

//! Define NetCDF variables corresponding to a diagnostic quantity.
void Diagnostic::define_impl(const PIO &file, IO_Type default_type) const {
  for (unsigned int j = 0; j < m_dof; ++j) {
    io::define_spatial_variable(m_vars[j], *m_grid, file,
                                default_type,
                                m_grid->ctx()->config()->get_string("output.variable_order"),
                                true);
  }
}

//! \brief A method for setting common variable attributes.
void Diagnostic::set_attrs(const std::string &long_name,
                           const std::string &standard_name,
                           const std::string &units,
                           const std::string &glaciological_units,
                           unsigned int N) {
  if (N >= m_dof) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "N (%d) >= m_dof (%d)", N, m_dof);
  }

  m_vars[N].set_string("pism_intent", "diagnostic");

  m_vars[N].set_string("long_name", long_name);

  m_vars[N].set_string("standard_name", standard_name);

  m_vars[N].set_string("units", units);

  if (not glaciological_units.empty()) {
    m_vars[N].set_string("glaciological_units", glaciological_units);
  }
}

IceModelVec::Ptr Diagnostic::compute() const {
  // use the name of the first variable
  std::vector<std::string> names;
  for (auto &v : m_vars) {
    names.push_back(v.get_name());
  }
  std::string all_names = join(names, ",");

  m_grid->ctx()->log()->message(3, "-  Computing %s...\n", all_names.c_str());
  IceModelVec::Ptr result = this->compute_impl();
  m_grid->ctx()->log()->message(3, "-  Done computing %s.\n", all_names.c_str());
  result->write_in_glaciological_units = true;
  return result;
}

TSDiagnostic::TSDiagnostic(IceGrid::ConstPtr g, const std::string &name)
  : m_grid(g),
    m_config(g->ctx()->config()),
    m_sys(g->ctx()->unit_system()),
    m_ts(new DiagnosticTimeseries(*g, name, g->ctx()->config()->get_string("time.dimension_name"))) {

  m_current_time = 0;

  m_ts->dimension_metadata().set_string("units", g->ctx()->time()->CF_units_string());
}

TSDiagnostic::~TSDiagnostic() {
}

void TSDiagnostic::update(double t0, double t1) {
  this->update_impl(t0, t1);
}


void TSDiagnostic::save(double a, double b) {
  m_ts->interp(a, b);
}

void TSDiagnostic::flush() {
  m_ts->flush();
}

void TSDiagnostic::init(const std::string &output_filename,
                        std::shared_ptr<std::vector<double>> requested_times) {
  m_ts->init(output_filename);
  m_times = requested_times;
}

const VariableMetadata &TSDiagnostic::metadata() const {
  return m_ts->metadata();
}

} // end of namespace pism

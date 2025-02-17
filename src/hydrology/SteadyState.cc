/* Copyright (C) 2019, 2020, 2021, 2023, 2024, 2025 PISM Authors
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

#include "pism/hydrology/SteadyState.hh"

#include <gsl/gsl_interp.h>     // gsl_interp_bsearch

#include "pism/hydrology/EmptyingProblem.hh"

#include "pism/util/Time.hh"    // time().current()
#include "pism/util/Profiling.hh"
#include "pism/util/Context.hh"
#include "pism/util/MaxTimestep.hh"

/* FIXMEs
 *
 * - At some later date I might want to support water input rates coming from the model.
 *   In that case I will have to accumulate input at every time step and then use the
 *   average (in time) when it's time to update the model.
 *
 * - IceModel::step() updates ice geometry before calling IceModel::hydrology_step(). This
 *   means that a hydrology model has to be able to provide its outputs before the first
 *   update() call. We save subglacial water flux to ensure that this is true for
 *   re-started runs, but in runs started using bootstrapping the first time step will see
 *   "bootstrapped" hydrology outputs (in this case: zero flux).
 *
 * - In this context (i.e. computing water flux using piecewise-in-time forcing data) it
 *    makes sense to update the flux at the beginning of an "update interval". In the case
 *    described above (water input coming from the model) we would want to update at the
 *    *end* of an interval.
 *
 */

namespace pism {
namespace hydrology {

void SteadyState::initialization_message() const {
  m_log->message(2,
                 "* Initializing the \"steady state\" subglacial hydrology model ...\n");
}

SteadyState::SteadyState(std::shared_ptr<const Grid> grid)
  : NullTransport(grid) {

  m_time_name = m_config->get_string("time.dimension_name") + "_hydrology_steady";
  m_t_last = time().current();
  m_update_interval = m_config->get_number("hydrology.steady.flux_update_interval", "seconds");
  m_t_eps = 1.0;
  m_bootstrap = false;

  m_emptying_problem.reset(new EmptyingProblem(grid));

  if (m_config->get_flag("hydrology.add_water_input_to_till_storage")) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "'steady' hydrology requires hydrology.add_water_input_to_till_storage == false");
  }

  if (m_config->get_string("hydrology.surface_input.file").empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "'steady' hydrology requires hydrology.surface_input.file");
  }
}

void SteadyState::update_impl(double t, double dt, const Inputs& inputs) {
  NullTransport::update_impl(t, dt, inputs);

  double t_next = m_t_last + max_timestep(m_t_last).value();

  if (t >= t_next or std::abs(t_next - t) < m_t_eps or
      m_bootstrap) {

    m_log->message(3, " Updating the steady-state subglacial water flux...\n");

    profiling().begin("steady_emptying");

    m_emptying_problem->update(*inputs.geometry,
                               inputs.no_model_mask,
                               m_surface_input_rate);

    profiling().end("steady_emptying");
    m_Q.copy_from(m_emptying_problem->flux());

    m_t_last = t;
    m_bootstrap = false;
  }
}

std::map<std::string, Diagnostic::Ptr> SteadyState::diagnostics_impl() const {
  auto hydro_diagnostics = NullTransport::diagnostics_impl();

  return combine(m_emptying_problem->diagnostics(), hydro_diagnostics);
}

MaxTimestep SteadyState::max_timestep_impl(double t) const {

  // compute the maximum time step coming from the forcing (water input rate)
  double dt_forcing = 0.0;
  if (not m_time.empty()) {

    // the right end point of the last time interval in the forcing file
    double t_last = m_time_bounds.back();

    double t_next = 0.0;
    if (t < m_time[0]) {
      // Allow stepping until the left end point of the first interval.
      t_next = m_time[0];
    } else if (t >= t_last) {
      // Went past the right end point of the last forcing intervals: no time step
      // restriction from forcing.
      t_next = t + m_update_interval;
    } else {
      // find the index k such that m_time[k] <= t < m_time[k + 1]
      size_t k = gsl_interp_bsearch(m_time.data(), t, 0, m_time.size());

      assert(m_time[k] <= t);
      assert(k + 1 == m_time.size() or t < m_time[k + 1]);

      t_next = m_time_bounds[2 * k + 1];

      if (std::abs(t_next - t) < m_t_eps) {
        // reached t_next; use the right end point of the next interval
        if (k + 1 < m_time.size()) {
          t_next = m_time_bounds[2 * (k + 1) + 1];
        } else {
          // No time step restriction from forcing. We set dt_forcing to m_update_interval
          // because dt_interval below will not exceed this, effectively selecting
          // dt_interval.
          t_next = t + m_update_interval;
        }
      }
    }

    dt_forcing = t_next - t;

    assert(dt_forcing > 0.0);
  } else {
    dt_forcing = m_update_interval;
  }

  // compute the maximum time step using the update interval
  double dt_interval = 0.0;
  {
    if (t < m_t_last) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "time %f is less than the previous time %f",
                                    t, m_t_last);
    }

    // Find the smallest time of the form m_t_last + k * m_update_interval that is greater
    // than t
    double k = ceil((t - m_t_last) / m_update_interval);

    double
      t_next = m_t_last + k * m_update_interval;

    dt_interval = t_next - t;

    if (dt_interval < m_t_eps) {
      dt_interval = m_update_interval;
    }
  }

  double dt = std::min(dt_forcing, dt_interval);

  auto dt_null = NullTransport::max_timestep_impl(t);
  if (dt_null.finite()) {
    dt = std::min(dt, dt_null.value());
  }

  return MaxTimestep(dt, "hydrology 'steady'");
}

void SteadyState::define_model_state_impl(const File& output) const {
  NullTransport::define_model_state_impl(output);

  if (not output.variable_exists(m_time_name)) {
    output.define_variable(m_time_name, io::PISM_DOUBLE, {});

    output.write_attribute(m_time_name, "long_name",
                        "time of the last update of the steady state subglacial water flux");
    output.write_attribute(m_time_name, "calendar", time().calendar());
    output.write_attribute(m_time_name, "units", time().units_string());
  }

  m_Q.define(output, io::PISM_DOUBLE);
}

void SteadyState::write_model_state_impl(const File& output) const {
  NullTransport::write_model_state_impl(output);

  output.write_variable(m_time_name, {0}, {1}, &m_t_last);
  m_Q.write(output);
}

void SteadyState::restart_impl(const File &input_file, int record) {
  NullTransport::restart_impl(input_file, record);

  init_time(m_config->get_string("hydrology.surface_input.file"));

  // Read m_t_last
  {
    if (input_file.variable_exists(m_time_name)) {
      input_file.read_variable(m_time_name, {0}, {1}, &m_t_last);
    } else {
      m_t_last = time().current();
    }
  }

  m_Q.read(input_file, record);

  regrid("hydrology 'steady'", m_Q, REGRID_WITHOUT_REGRID_VARS);
}

void SteadyState::bootstrap_impl(const File &input_file,
                                 const array::Scalar &ice_thickness) {
  NullTransport::bootstrap_impl(input_file, ice_thickness);

  init_time(m_config->get_string("hydrology.surface_input.file"));

  // Read m_t_last
  {
    if (input_file.variable_exists(m_time_name)) {
      input_file.read_variable(m_time_name, {0}, {1}, &m_t_last);
    } else {
      m_t_last = time().current();
    }
  }

  // Read water flux
  if (input_file.variable_exists(m_Q.metadata().get_name())) {
    // Regrid from the input file.
    m_Q.regrid(input_file, io::Default::Nil());

    // Allow regridding from a different file.
    regrid("hydrology 'steady'", m_Q, REGRID_WITHOUT_REGRID_VARS);
  } else {
    // No water flux in the input file; try regridding from a different file.
    auto n = m_Q.state_counter();

    regrid("hydrology 'steady'", m_Q, REGRID_WITHOUT_REGRID_VARS);

    if (n == m_Q.state_counter()) {
      // Regridding did not fill m_Q: we need to bootstrap during the first update_impl()
      // call.
      m_bootstrap = true;
    }
  }
}

void SteadyState::init_impl(const array::Scalar &W_till,
                            const array::Scalar &W,
                            const array::Scalar &P) {
  NullTransport::init_impl(W_till, W, P);

  m_Q.set(0.0);

  m_bootstrap = true;
}


/*!
 * Read time bounds corresponding to the water input rate in the forcing file.
 *
 * These times are used to compute the maximum time step the model can take while still
 * capturing temporal variability of the forcing.
 */
void SteadyState::init_time(const std::string &input_file) {

  std::string variable_name = "water_input_rate";

  File file(m_grid->com, input_file, io::PISM_GUESS, io::PISM_READONLY);

  auto time_name = io::time_dimension(m_grid->ctx()->unit_system(),
                                      file, variable_name);

  if (time_name.empty()) {
    // Water input rate is time-independent. m_time and m_time_bounds are left empty.
    return;
  }

  std::string bounds_name = file.read_text_attribute(time_name, "bounds");

  if (bounds_name.empty()) {
    // no time bounds attribute
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Variable '%s' does not have the time_bounds attribute.\n"
                                  "Cannot use time-dependent water input rate without time bounds.",
                                  time_name.c_str());
  }

  // read time bounds data from a file
  m_time_bounds =
      io::read_bounds(file, bounds_name, time().units_string(), m_grid->ctx()->unit_system());

  // time bounds data overrides the time variable: we make t[j] be the
  // left end-point of the j-th interval
  m_time.resize(m_time_bounds.size() / 2);
  for (unsigned int k = 0; k < m_time.size(); ++k) {
    m_time[k] = m_time_bounds[2*k];
  }

  if (not is_increasing(m_time)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "time bounds in %s are invalid", input_file.c_str());
  }
}

} // end of namespace hydrology
} // end of namespace pism

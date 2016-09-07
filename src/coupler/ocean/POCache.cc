/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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

#include <algorithm>
#include <cassert>

#include "POCache.hh"
#include "base/util/PISMTime.hh"
#include "base/util/pism_options.hh"
#include "base/util/IceGrid.hh"

#include "base/util/error_handling.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace ocean {

Cache::Cache(IceGrid::ConstPtr g, OceanModel* in)
  : OceanModifier(g, in) {

  m_next_update_time = m_grid->ctx()->time()->current();
  m_update_interval_years = 10;

  m_shelf_base_mass_flux.create(m_grid, "shelfbmassflux", WITHOUT_GHOSTS);
  m_shelf_base_mass_flux.set_attrs("climate_state",
                                   "ice mass flux from ice shelf base"
                                   " (positive flux is loss from ice shelf)",
                                   "kg m-2 s-1", "");
  m_shelf_base_mass_flux.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_shelf_base_mass_flux.write_in_glaciological_units = true;

  m_shelf_base_temperature.create(m_grid, "shelfbtemp", WITHOUT_GHOSTS);
  m_shelf_base_temperature.set_attrs("climate_state",
                                     "absolute temperature at ice shelf base",
                                     "K", "");
  m_melange_back_pressure_fraction.create(m_grid,"melange_back_pressure_fraction",
                                          WITHOUT_GHOSTS);
  m_melange_back_pressure_fraction.set_attrs("climate_state",
                                             "melange back pressure fraction",
                                             "1", "");
}

Cache::~Cache() {
  // empty
}


void Cache::init_impl() {
  int update_interval = m_update_interval_years;

  m_input_model->init();

  m_log->message(2,
             "* Initializing the 'caching' ocean model modifier...\n");

  update_interval = options::Integer("-ocean_cache_update_interval",
                                     "Interval (in years) between ocean model updates",
                                     update_interval);

  if (update_interval <= 0) {
    throw RuntimeError("-ocean_cache_update_interval has to be strictly positive.");
  }

  m_update_interval_years = update_interval;
  m_next_update_time = m_grid->ctx()->time()->current();
}

void Cache::update_impl(double my_t, double my_dt) {
  // ignore my_dt and always use 1 year long time-steps when updating
  // an input model
  (void) my_dt;

  if (my_t >= m_next_update_time ||
      fabs(my_t - m_next_update_time) < 1.0) {

    double
      one_year_from_now = m_grid->ctx()->time()->increment_date(my_t, 1.0),
      update_dt         = one_year_from_now - my_t;

    assert(update_dt > 0.0);

    m_input_model->update(my_t, update_dt);

    m_next_update_time = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                   m_update_interval_years);

    m_sea_level = m_input_model->sea_level_elevation();
    m_input_model->shelf_base_temperature(m_shelf_base_temperature);
    m_input_model->shelf_base_mass_flux(m_shelf_base_mass_flux);
    m_input_model->melange_back_pressure_fraction(m_melange_back_pressure_fraction);
  }
}


void Cache::sea_level_elevation_impl(double &result) {
  result = m_sea_level;
}

void Cache::shelf_base_temperature_impl(IceModelVec2S &result) {
  result.copy_from(m_shelf_base_temperature);
}

void Cache::shelf_base_mass_flux_impl(IceModelVec2S &result) {
  result.copy_from(m_shelf_base_mass_flux);
}

void Cache::melange_back_pressure_fraction_impl(IceModelVec2S &result) {
  result.copy_from(m_melange_back_pressure_fraction);
}


void Cache::define_variables_impl(const std::set<std::string> &vars_input, const PIO &nc,
                                         IO_Type nctype) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, m_shelf_base_mass_flux.metadata().get_string("short_name"))) {
    m_shelf_base_mass_flux.define(nc, nctype);
    vars.erase(m_shelf_base_mass_flux.metadata().get_string("short_name"));
  }

  if (set_contains(vars, m_shelf_base_temperature.metadata().get_string("short_name"))) {
    m_shelf_base_temperature.define(nc, nctype);
    vars.erase(m_shelf_base_temperature.metadata().get_string("short_name"));
  }

  if (set_contains(vars, m_melange_back_pressure_fraction.metadata().get_string("short_name"))) {
    m_melange_back_pressure_fraction.define(nc, nctype);
    vars.erase(m_melange_back_pressure_fraction.metadata().get_string("short_name"));
  }

  m_input_model->define_variables(vars, nc, nctype);
}

void Cache::write_variables_impl(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, m_shelf_base_mass_flux.metadata().get_string("short_name"))) {
    m_shelf_base_mass_flux.write(nc);
    vars.erase(m_shelf_base_mass_flux.metadata().get_string("short_name"));
  }

  if (set_contains(vars, m_shelf_base_temperature.metadata().get_string("short_name"))) {
    m_shelf_base_temperature.write(nc);
    vars.erase(m_shelf_base_temperature.metadata().get_string("short_name"));
  }

  if (set_contains(vars, m_melange_back_pressure_fraction.metadata().get_string("short_name"))) {
    m_melange_back_pressure_fraction.write(nc);
    vars.erase(m_melange_back_pressure_fraction.metadata().get_string("short_name"));
  }

  m_input_model->write_variables(vars, nc);
}

MaxTimestep Cache::max_timestep_impl(double t) {
  double dt = m_next_update_time - t;

  // if we got very close to the next update time, set time step
  // length to the interval between updates
  if (dt < 1.0) {
    double update_time_after_next = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                                m_update_interval_years);

    dt = update_time_after_next - m_next_update_time;
    assert(dt > 0.0);
  }

  MaxTimestep input_max_timestep = m_input_model->max_timestep(t);
  if (input_max_timestep.is_finite()) {
    return std::min(input_max_timestep, MaxTimestep(dt));
  } else {
    return MaxTimestep(dt);
  }
}

} // end of namespace ocean
} // end of namespace pism

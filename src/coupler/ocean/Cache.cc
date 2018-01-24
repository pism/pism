/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "Cache.hh"
#include "pism/util/Time.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/IceGrid.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace ocean {

Cache::Cache(IceGrid::ConstPtr g, OceanModel* in)
  : OceanModifier(g, in) {

  m_next_update_time = m_grid->ctx()->time()->current();
  m_update_interval_years = 10;

  m_shelf_base_mass_flux.create(m_grid, "effective_shelf_base_mass_flux", WITHOUT_GHOSTS);
  m_shelf_base_mass_flux.set_attrs("climate_state",
                                   "ice mass flux from ice shelf base"
                                   " (positive flux is loss from ice shelf)",
                                   "kg m-2 s-1", "");
  m_shelf_base_mass_flux.metadata().set_string("glaciological_units", "kg m-2 year-1");
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
    throw RuntimeError(PISM_ERROR_LOCATION, "-ocean_cache_update_interval has to be strictly positive.");
  }

  m_update_interval_years = update_interval;
  m_next_update_time = m_grid->ctx()->time()->current();
}

void Cache::update_impl(double t, double dt) {
  // ignore dt and always use 1 year long time-steps when updating
  // an input model
  (void) dt;

  if (t >= m_next_update_time ||
      fabs(t - m_next_update_time) < 1.0) {

    double
      one_year_from_now = m_grid->ctx()->time()->increment_date(t, 1.0),
      update_dt         = one_year_from_now - t;

    assert(update_dt > 0.0);

    m_input_model->update(t, update_dt);

    m_next_update_time = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                   m_update_interval_years);

    m_sea_level = m_input_model->sea_level_elevation();

    m_melange_back_pressure_fraction.copy_from(m_input_model->melange_back_pressure_fraction());

    m_shelf_base_temperature.copy_from(m_input_model->shelf_base_temperature());

    m_shelf_base_mass_flux.copy_from(m_input_model->shelf_base_mass_flux());
  }
}

MaxTimestep Cache::max_timestep_impl(double t) const {
  double dt = m_next_update_time - t;

  // if we got very close to the next update time, set time step
  // length to the interval between updates
  if (dt < 1.0) {
    double update_time_after_next = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                                m_update_interval_years);

    dt = update_time_after_next - m_next_update_time;
    assert(dt > 0.0);
  }

  MaxTimestep cache_dt(dt, "ocean cache");

  MaxTimestep input_max_timestep = m_input_model->max_timestep(t);
  if (input_max_timestep.finite()) {
    return std::min(input_max_timestep, cache_dt);
  } else {
    return cache_dt;
  }
}

} // end of namespace ocean
} // end of namespace pism

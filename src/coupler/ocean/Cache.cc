/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2020, 2021 PISM Authors
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

#include <algorithm>            // std::min
#include <cassert>
#include <cmath>

#include "Cache.hh"
#include "pism/util/Time.hh"
#include "pism/util/Grid.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace ocean {

Cache::Cache(std::shared_ptr<const Grid> g, std::shared_ptr<OceanModel> in)
  : OceanModel(g, in) {

  m_next_update_time = time().current();
  m_update_interval_years = m_config->get_number("ocean.cache.update_interval", "seconds");

  // use the current year length (according to the selected calendar) to convert update
  // interval length into years
  m_update_interval_years = time().convert_time_interval(m_update_interval_years, "years");

  if (m_update_interval_years <= 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "ocean.cache.update_interval has to be strictly positive (got %f)",
                                  m_update_interval_years);
  }

  {
    m_shelf_base_temperature = allocate_shelf_base_temperature(g);
    m_shelf_base_mass_flux   = allocate_shelf_base_mass_flux(g);
    m_water_column_pressure  = allocate_water_column_pressure(g);
  }
}

void Cache::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2,
                 "* Initializing the 'caching' ocean model modifier...\n");

  m_next_update_time = time().current();
}

void Cache::update_impl(const Geometry &geometry, double t, double dt) {
  // ignore dt and always use 1 year long time-steps when updating
  // an input model
  (void) dt;

  double time_resolution = m_config->get_number("time_stepping.resolution", "seconds");

  if (t >= m_next_update_time or
      fabs(t - m_next_update_time) < time_resolution) {

    double
      one_year_from_now = time().increment_date(t, 1.0),
      update_dt         = one_year_from_now - t;

    assert(update_dt > 0.0);

    m_input_model->update(geometry, t, update_dt);

    m_next_update_time = time().increment_date(m_next_update_time,
                                                               m_update_interval_years);

    m_water_column_pressure->copy_from(m_input_model->average_water_column_pressure());

    m_shelf_base_temperature->copy_from(m_input_model->shelf_base_temperature());

    m_shelf_base_mass_flux->copy_from(m_input_model->shelf_base_mass_flux());
  }
}

MaxTimestep Cache::max_timestep_impl(double t) const {
  double dt = m_next_update_time - t;

  double time_resolution = m_config->get_number("time_stepping.resolution", "seconds");

  // if we got very close to the next update time, set time step
  // length to the interval between updates
  if (dt < time_resolution) {
    double update_time_after_next = time().increment_date(m_next_update_time,
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

const array::Scalar& Cache::shelf_base_temperature_impl() const {
  return *m_shelf_base_temperature;
}

const array::Scalar& Cache::shelf_base_mass_flux_impl() const {
  return *m_shelf_base_mass_flux;
}

const array::Scalar& Cache::average_water_column_pressure_impl() const {
  return *m_water_column_pressure;
}


} // end of namespace ocean
} // end of namespace pism

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

#include <cassert>
#include <algorithm>            // for std::min

#include "Cache.hh"
#include "pism/util/Time.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace surface {

Cache::Cache(IceGrid::ConstPtr grid, std::shared_ptr<SurfaceModel> in)
  : SurfaceModel(grid, in) {

  m_next_update_time = m_grid->ctx()->time()->current();
  m_update_interval_years = m_config->get_double("surface.cache.update_interval", "years");

  if (m_update_interval_years <= 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "surface.cache.update_interval has to be strictly positive.");
  }

  {
    m_mass_flux             = allocate_mass_flux(grid);
    m_temperature           = allocate_temperature(grid);
    m_liquid_water_fraction = allocate_liquid_water_fraction(grid);
    m_layer_mass            = allocate_layer_mass(grid);
    m_layer_thickness       = allocate_layer_thickness(grid);
  }
}

Cache::~Cache() {
  // empty
}

void Cache::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2, "* Initializing the 'caching' surface model modifier...\n");

  m_next_update_time = m_grid->ctx()->time()->current();
}

void Cache::update_impl(const Geometry &geometry, double t, double dt) {
  // ignore dt and always use 1 year long time-steps when updating
  // the input model
  (void) dt;

  if (t >= m_next_update_time or fabs(t - m_next_update_time) < 1.0) {

    double
      one_year_from_now = m_grid->ctx()->time()->increment_date(t, 1.0),
      update_dt         = one_year_from_now - t;

    assert(update_dt > 0.0);

    m_input_model->update(geometry, t, update_dt);

    m_next_update_time = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                               m_update_interval_years);

    // store outputs of the input model
    m_mass_flux->copy_from(m_input_model->mass_flux());
    m_temperature->copy_from(m_input_model->temperature());
    m_liquid_water_fraction->copy_from(m_input_model->liquid_water_fraction());
    m_layer_mass->copy_from(m_input_model->layer_mass());
    m_layer_thickness->copy_from(m_input_model->layer_thickness());
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

  assert(m_input_model != NULL);

  MaxTimestep cache_dt(dt, "surface cache");

  MaxTimestep input_max_timestep = m_input_model->max_timestep(t);
  if (input_max_timestep.finite()) {
    return std::min(input_max_timestep, cache_dt);
  } else {
    return cache_dt;
  }
}

const IceModelVec2S &Cache::layer_thickness_impl() const {
  return *m_layer_thickness;
}

const IceModelVec2S &Cache::mass_flux_impl() const {
  return *m_mass_flux;
}

const IceModelVec2S &Cache::temperature_impl() const {
  return *m_temperature;
}

const IceModelVec2S &Cache::liquid_water_fraction_impl() const {
  return *m_liquid_water_fraction;
}

const IceModelVec2S &Cache::layer_mass_impl() const {
  return *m_layer_mass;
}

} // end of namespace surface
} // end of namespace pism

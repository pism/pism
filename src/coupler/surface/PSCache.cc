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

#include <cassert>
#include <algorithm>            // for std::min

#include "PSCache.hh"
#include "base/util/PISMTime.hh"
#include "base/util/pism_options.hh"
#include "base/util/IceGrid.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace surface {

Cache::Cache(IceGrid::ConstPtr g, SurfaceModel* in)
  : SurfaceModifier(g, in) {

  m_next_update_time = m_grid->ctx()->time()->current();
  m_update_interval_years = 10;

  m_mass_flux.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  m_mass_flux.set_attrs("climate_state",
                        "surface mass balance (accumulation/ablation) rate",
                        "kg m-2 s-1",
                        "land_ice_surface_specific_mass_balance_flux");
  m_mass_flux.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_mass_flux.write_in_glaciological_units = true;

  m_temperature.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
  m_temperature.set_attrs("climate_state",
                          "ice temperature at the ice surface",
                          "K", "");

  m_liquid_water_fraction.create(m_grid, "ice_surface_liquid_water_fraction", WITHOUT_GHOSTS);
  m_liquid_water_fraction.set_attrs("diagnostic",
                                    "ice surface liquid water fraction", "1", "");

  m_mass_held_in_surface_layer.create(m_grid, "mass_held_in_surface_layer", WITHOUT_GHOSTS);
  m_mass_held_in_surface_layer.set_attrs("diagnostic",
                                         "mass held in surface layer", "kg", "");

  m_surface_layer_thickness.create(m_grid, "surface_layer_thickness", WITHOUT_GHOSTS);
  m_surface_layer_thickness.set_attrs("diagnostic",
                                      "surface layer thickness", "1", "");
}

Cache::~Cache() {
  // empty
}


void Cache::init_impl() {
  int update_interval = m_update_interval_years;

  m_input_model->init();

  m_log->message(2,
             "* Initializing the 'caching' surface model modifier...\n");

  update_interval = options::Integer("-surface_cache_update_interval",
                                     "Interval (in years) between surface model updates",
                                     update_interval);

  if (update_interval <= 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "-surface_cache_update_interval has to be strictly positive.");
  }

  m_update_interval_years = update_interval;
  m_next_update_time = m_grid->ctx()->time()->current();
}

void Cache::update_impl(double t, double dt) {
  // ignore dt and always use 1 year long time-steps when updating
  // an input model
  (void) dt;

  if (t >= m_next_update_time or fabs(t - m_next_update_time) < 1.0) {

    double
      one_year_from_now = m_grid->ctx()->time()->increment_date(t, 1.0),
      update_dt         = one_year_from_now - t;

    assert(update_dt > 0.0);

    m_input_model->update(t, update_dt);

    m_next_update_time = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                               m_update_interval_years);

    m_input_model->ice_surface_mass_flux(m_mass_flux);
    m_input_model->ice_surface_temperature(m_temperature);
    m_input_model->ice_surface_liquid_water_fraction(m_liquid_water_fraction);
    m_input_model->mass_held_in_surface_layer(m_mass_held_in_surface_layer);
    m_input_model->surface_layer_thickness(m_surface_layer_thickness);
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

void Cache::ice_surface_mass_flux_impl(IceModelVec2S &result) const {
  result.copy_from(m_mass_flux);
}

void Cache::ice_surface_temperature_impl(IceModelVec2S &result) const {
  result.copy_from(m_temperature);
}

void Cache::ice_surface_liquid_water_fraction_impl(IceModelVec2S &result) const {
  result.copy_from(m_liquid_water_fraction);
}

void Cache::mass_held_in_surface_layer_impl(IceModelVec2S &result) const {
  result.copy_from(m_mass_held_in_surface_layer);
}

void Cache::surface_layer_thickness_impl(IceModelVec2S &result) const {
  result.copy_from(m_surface_layer_thickness);
}

} // end of namespace surface
} // end of namespace pism

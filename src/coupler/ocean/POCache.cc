/* Copyright (C) 2013, 2014 PISM Authors
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

#include "POCache.hh"
#include "PISMTime.hh"
#include "pism_options.hh"
#include <algorithm>
#include <cassert>

namespace pism {

POCache::POCache(IceGrid &g, const Config &conf, OceanModel* in)
  : POModifier(g, conf, in) {

  m_next_update_time = grid.time->current();
  m_update_interval_years = 10;

  PetscErrorCode ierr = allocate_POCache();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: failed to allocate storage for PSCache.\n");
    PISMEnd();
  }
}

PetscErrorCode POCache::allocate_POCache() {
  PetscErrorCode ierr;

  ierr = m_shelf_base_mass_flux.create(grid, "shelfbmassflux", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = m_shelf_base_mass_flux.set_attrs("climate_state",
                                          "ice mass flux from ice shelf base"
                                          " (positive flux is loss from ice shelf)",
                                          "kg m-2 s-1", ""); CHKERRQ(ierr);
  ierr = m_shelf_base_mass_flux.set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);
  m_shelf_base_mass_flux.write_in_glaciological_units = true;

  ierr = m_shelf_base_temperature.create(grid, "shelfbtemp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = m_shelf_base_temperature.set_attrs("climate_state",
                                            "absolute temperature at ice shelf base",
                                            "K", ""); CHKERRQ(ierr);
  ierr = m_melange_back_pressure_fraction.create(grid,"melange_back_pressure_fraction",
                                                 WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = m_melange_back_pressure_fraction.set_attrs("climate_state",
                                                    "melange back pressure fraction",
                                                    "1", ""); CHKERRQ(ierr);


  return 0;
}

POCache::~POCache() {
  // empty
}


PetscErrorCode POCache::init(Vars &vars) {
  PetscErrorCode ierr;
  int update_interval = m_update_interval_years;
  bool flag;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the 'caching' ocean model modifier...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "-ocean ...,cache options", ""); CHKERRQ(ierr);
  {
    ierr = OptionsInt("-ocean_cache_update_interval",
                          "Interval (in years) between ocean model updates",
                          update_interval, flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (update_interval <= 0) {
    PetscPrintf(grid.com,
                "PISM ERROR: -ocean_cache_update_interval has to be strictly positive.\n");
    PISMEnd();
  }

  m_update_interval_years = update_interval;
  m_next_update_time = grid.time->current();

  return 0;
}

PetscErrorCode POCache::update(double my_t, double my_dt) {
  PetscErrorCode ierr;

  // ignore my_dt and always use 1 year long time-steps when updating
  // an input model
  (void) my_dt;

  if (my_t >= m_next_update_time ||
      fabs(my_t - m_next_update_time) < 1.0) {

    double
      one_year_from_now = grid.time->increment_date(my_t, 1.0),
      update_dt         = one_year_from_now - my_t;

    assert(update_dt > 0.0);

    ierr = input_model->update(my_t, update_dt); CHKERRQ(ierr);

    m_next_update_time = grid.time->increment_date(m_next_update_time,
                                                   m_update_interval_years);

    ierr = input_model->sea_level_elevation(m_sea_level);                 CHKERRQ(ierr); 
    ierr = input_model->shelf_base_temperature(m_shelf_base_temperature); CHKERRQ(ierr); 
    ierr = input_model->shelf_base_mass_flux(m_shelf_base_mass_flux);     CHKERRQ(ierr); 
    ierr = input_model->melange_back_pressure_fraction(m_melange_back_pressure_fraction); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode POCache::sea_level_elevation(double &result) {
  result = m_sea_level;
  return 0;
}

PetscErrorCode POCache::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = m_shelf_base_temperature.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POCache::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = m_shelf_base_mass_flux.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POCache::melange_back_pressure_fraction(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = m_melange_back_pressure_fraction.copy_to(result); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode POCache::define_variables(std::set<std::string> vars, const PIO &nc,
                                         IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, m_shelf_base_mass_flux.metadata().get_string("short_name"))) {
    ierr = m_shelf_base_mass_flux.define(nc, nctype); CHKERRQ(ierr);
    vars.erase(m_shelf_base_mass_flux.metadata().get_string("short_name"));
  }

  if (set_contains(vars, m_shelf_base_temperature.metadata().get_string("short_name"))) {
    ierr = m_shelf_base_temperature.define(nc, nctype); CHKERRQ(ierr);
    vars.erase(m_shelf_base_temperature.metadata().get_string("short_name"));
  }

  if (set_contains(vars, m_melange_back_pressure_fraction.metadata().get_string("short_name"))) {
    ierr = m_melange_back_pressure_fraction.define(nc, nctype); CHKERRQ(ierr);
    vars.erase(m_melange_back_pressure_fraction.metadata().get_string("short_name"));
  }

  ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POCache::write_variables(std::set<std::string> vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, m_shelf_base_mass_flux.metadata().get_string("short_name"))) {
    ierr = m_shelf_base_mass_flux.write(nc); CHKERRQ(ierr);
    vars.erase(m_shelf_base_mass_flux.metadata().get_string("short_name"));
  }

  if (set_contains(vars, m_shelf_base_temperature.metadata().get_string("short_name"))) {
    ierr = m_shelf_base_temperature.write(nc); CHKERRQ(ierr);
    vars.erase(m_shelf_base_temperature.metadata().get_string("short_name"));
  }

  if (set_contains(vars, m_melange_back_pressure_fraction.metadata().get_string("short_name"))) {
    ierr = m_melange_back_pressure_fraction.write(nc); CHKERRQ(ierr);
    vars.erase(m_melange_back_pressure_fraction.metadata().get_string("short_name"));
  }

  ierr = input_model->write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POCache::max_timestep(double t, double &dt, bool &restrict) {
  dt       = m_next_update_time - t;
  restrict = true;

  // if we got very close to the next update time, set time step
  // length to the interval between updates
  if (dt < 1.0) {
    double update_time_after_next = grid.time->increment_date(m_next_update_time,
                                                              m_update_interval_years);

    dt = update_time_after_next - m_next_update_time;
    assert(dt > 0.0);
  }

  bool input_restrict = false;
  double input_model_dt = 0.0;
  assert(input_model != NULL);

  PetscErrorCode ierr = input_model->max_timestep(t, input_model_dt, input_restrict); CHKERRQ(ierr);
  if (input_restrict == true) {
    dt = std::min(input_model_dt, dt);
  }

  return 0;
}

} // end of namespace pism

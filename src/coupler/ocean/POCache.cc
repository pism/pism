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

POCache::POCache(IceGrid &g, const PISMConfig &conf, PISMOceanModel* in)
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
                                          "m s-1", ""); CHKERRQ(ierr);
  ierr = m_shelf_base_mass_flux.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  m_shelf_base_mass_flux.write_in_glaciological_units = true;

  ierr = m_shelf_base_temperature.create(grid, "shelfbtemp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = m_shelf_base_temperature.set_attrs("climate_state",
                                            "absolute temperature at ice shelf base",
                                            "K", ""); CHKERRQ(ierr);

  return 0;
}

POCache::~POCache() {
  // empty
}


PetscErrorCode POCache::init(PISMVars &vars) {
  PetscErrorCode ierr;
  int update_interval = m_update_interval_years;
  bool flag;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the 'caching' ocean model modifier...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "-ocean ...,cache options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsInt("-ocean_cache_update_interval",
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

PetscErrorCode POCache::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;
  if (my_t + my_dt > m_next_update_time) {
    ierr = input_model->update(my_t + 0.5*my_dt,
                               grid.convert(1.0, "year", "seconds")); CHKERRQ(ierr);

    m_next_update_time = grid.time->increment_date(m_next_update_time,
                                                   m_update_interval_years);

    ierr = input_model->sea_level_elevation(m_sea_level);                 CHKERRQ(ierr); 
    ierr = input_model->shelf_base_temperature(m_shelf_base_temperature); CHKERRQ(ierr); 
    ierr = input_model->shelf_base_mass_flux(m_shelf_base_mass_flux);     CHKERRQ(ierr); 
  }

  return 0;
}


PetscErrorCode POCache::sea_level_elevation(PetscReal &result) {
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


PetscErrorCode POCache::define_variables(std::set<std::string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, m_shelf_base_mass_flux.metadata().get_string("short_name"))) {
    ierr = m_shelf_base_mass_flux.define(nc, nctype); CHKERRQ(ierr);
    vars.erase(m_shelf_base_mass_flux.metadata().get_string("short_name"));
  }

  if (set_contains(vars, m_shelf_base_temperature.metadata().get_string("short_name"))) {
    ierr = m_shelf_base_temperature.define(nc, nctype); CHKERRQ(ierr);
    vars.erase(m_shelf_base_temperature.metadata().get_string("short_name"));
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

  ierr = input_model->write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}


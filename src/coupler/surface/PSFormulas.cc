/* Copyright (C) 2014 PISM Authors
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

#include "PSFormulas.hh"
#include "PISMAtmosphere.hh"

namespace pism {

PSFormulas::PSFormulas(IceGrid &g, const Config &conf)
  : SurfaceModel(g, conf) {
  PetscErrorCode ierr = allocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: memory allocation failed");
    PISMEnd();
  }
}

PetscErrorCode PSFormulas::allocate() {
  PetscErrorCode ierr;

  ierr = m_climatic_mass_balance.create(grid, "climatic_mass_balance", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = m_climatic_mass_balance.set_attrs("internal",
                                           "ice-equivalent surface mass balance (accumulation/ablation) rate",
                                           "kg m-2 s-1",
                                           "land_ice_surface_specific_mass_balance_flux"); CHKERRQ(ierr);
  ierr = m_climatic_mass_balance.set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);
  m_climatic_mass_balance.write_in_glaciological_units = true;
  m_climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // annual mean air temperature at "ice surface", at level below all
  // firn processes (e.g. "10 m" or ice temperatures)
  ierr = m_ice_surface_temp.create(grid, "ice_surface_temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = m_ice_surface_temp.set_attrs("internal",
                                      "annual average ice surface temperature, below firn processes",
                                      "K", ""); CHKERRQ(ierr);
  return 0;
}

PSFormulas::~PSFormulas() {
  // empty
}


void PSFormulas::attach_atmosphere_model(AtmosphereModel *input) {
  delete input;
}

PetscErrorCode PSFormulas::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = m_climatic_mass_balance.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PSFormulas::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = m_ice_surface_temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

void PSFormulas::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  (void) keyword;

  result.insert(m_climatic_mass_balance.metadata().get_name());
  result.insert(m_ice_surface_temp.metadata().get_name());
}

PetscErrorCode PSFormulas::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                             IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, m_climatic_mass_balance.metadata().get_name())) {
    ierr = m_climatic_mass_balance.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, m_ice_surface_temp.metadata().get_name())) {
    ierr = m_ice_surface_temp.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PSFormulas::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, m_climatic_mass_balance.metadata().get_name())) {
    ierr = m_climatic_mass_balance.write(nc); CHKERRQ(ierr);
  }

  if (set_contains(vars, m_ice_surface_temp.metadata().get_name())) {
    ierr = m_ice_surface_temp.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

} // end of namespace pism

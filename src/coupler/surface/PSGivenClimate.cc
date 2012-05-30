// Copyright (C) 2011, 2012 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "PSGivenClimate.hh"
#include "IceGrid.hh"

PetscErrorCode PSGivenClimate::init(PISMVars &) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the surface model reading temperature at the top of the ice\n"
                    "  and ice surface mass flux from a file...\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr);

  ierr = temp.set_attrs("climate_forcing",
                        "temperature of the ice at the ice surface but below firn processes",
                        "Kelvin", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing",
                       "ice-equivalent surface mass balance (accumulation/ablation) rate",
                       "m s-1", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);
  ierr = mass_flux.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  mass_flux.write_in_glaciological_units = true;

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  // read time-independent data right away:
  if (temp.get_n_records() == 1 && mass_flux.get_n_records() == 1) {
    ierr = update(grid.time->current(), 0); CHKERRQ(ierr); // dt is irrelevant
  }

  return 0;
}

PetscErrorCode PSGivenClimate::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = mass_flux.at_time(t); CHKERRQ(ierr);
  ierr = temp.at_time(t); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSGivenClimate::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = mass_flux.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PSGivenClimate::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

// Copyright (C) 2011 Constantine Khroulev
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

#include "PODirectForcing.hh"

PetscErrorCode PODirectForcing::init(PISMVars &) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the ocean model reading base of the shelf temperature\n"
                    "  and sub-shelf mass flux from a file...\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", ""); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr);

  ierr = temp.set_attrs("climate_forcing",
                        "absolute temperature at ice shelf base",
                        "Kelvin", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing",
                       "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                       "m s-1", ""); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
                    "    reading boundary conditions from %s ...\n",
                    filename.c_str()); CHKERRQ(ierr);

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PODirectForcing::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr = update_internal(t_years, dt_years); CHKERRQ(ierr);

  if (enable_time_averaging) {
    ierr = mass_flux.average(t, dt); CHKERRQ(ierr); 
    ierr = temp.average(t, dt); CHKERRQ(ierr); 
  } else {
    ierr = mass_flux.get_record_years(t); CHKERRQ(ierr);
    ierr = temp.get_record_years(t); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PODirectForcing::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PODirectForcing::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = mass_flux.copy_to(result); CHKERRQ(ierr);
  return 0;
}


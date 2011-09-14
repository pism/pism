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

#include "PSDirectAnomalies.hh"

PetscErrorCode PSDirectAnomalies::init(PISMVars &vars) {
  PetscErrorCode ierr;
  string input_file;
  bool regrid;
  int start;

  ierr = PSDirectForcing::init(vars); CHKERRQ(ierr);

  // create special variables
  ierr = mass_flux_0.create(grid, "mass_flux_0", false); CHKERRQ(ierr);
  ierr = mass_flux_0.set_attrs("internal", "surface mass flux at the beginning of a run",
                               "m s-1", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);

  ierr = mass_flux_in.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = mass_flux_in.set_attrs("model_state", "surface mass flux to apply anomalies to",
                                    "m s-1", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);

  ierr = find_pism_input(input_file, regrid, start); CHKERRQ(ierr);

  if (regrid) {
    ierr = mass_flux_in.regrid(input_file.c_str(), true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = mass_flux_in.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }
  
  // get the mass balance at the beginning of the run:
  ierr = PSDirectForcing::update(grid.start_year, 0); CHKERRQ(ierr);

  ierr = mass_flux.copy_to(mass_flux_0); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSDirectAnomalies::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  ierr = PSDirectForcing::update(t_years, dt_years); CHKERRQ(ierr);

  ierr = mass_flux.begin_access(); CHKERRQ(ierr);
  ierr = mass_flux_0.begin_access(); CHKERRQ(ierr);
  ierr = mass_flux_in.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      mass_flux(i, j) = mass_flux(i, j) - mass_flux_0(i, j) + mass_flux_in(i, j);
    }
  }

  ierr = mass_flux_in.end_access(); CHKERRQ(ierr);
  ierr = mass_flux_0.end_access(); CHKERRQ(ierr);
  ierr = mass_flux.end_access(); CHKERRQ(ierr);

  return 0;
}


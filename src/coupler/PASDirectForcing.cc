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

// Implementation of surface and atmosphere models reading spatially-variable
// time-dependent B.C. data from a file (-surface given and -atmosphere given).

#include "PASDirectForcing.hh"

PetscErrorCode PSDirectForcing::init(PISMVars &) {
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

  ierr = verbPrintf(2,grid.com,
                    "    reading boundary conditions from %s ...\n",
                    filename.c_str()); CHKERRQ(ierr);

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSDirectForcing::update(PetscReal t_years, PetscReal dt_years) {
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

PetscErrorCode PSDirectForcing::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = mass_flux.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PSDirectForcing::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PADirectForcing::init(PISMVars &) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the atmosphere model reading near-surface air temperature\n"
                    "  and ice-equivalent precipitation from a file...\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", ""); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr);

  ierr = temp.set_attrs("climate_forcing", "near-surface air temperature",
                        "Kelvin", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing", "ice-equivalent precipitation rate",
                       "m s-1", ""); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "    reading boundary conditions from %s ...\n",
                    filename.c_str()); CHKERRQ(ierr);

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PADirectForcing::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr = update_internal(t_years, dt_years); CHKERRQ(ierr);

  if (enable_time_averaging) {
    ierr = mass_flux.average(t, dt); CHKERRQ(ierr); 
    ierr = temp.average(t, 1.0); CHKERRQ(ierr); // compute the "mean annual" temperature
  } else {
    ierr = mass_flux.get_record_years(t); CHKERRQ(ierr);
    ierr = temp.get_record_years(t); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PADirectForcing::mean_precip(IceModelVec2S &result) {
  PetscErrorCode ierr = mass_flux.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PADirectForcing::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr = temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PADirectForcing::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr = temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PADirectForcing::begin_pointwise_access() {
  PetscErrorCode ierr = temp.begin_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PADirectForcing::end_pointwise_access() {
  PetscErrorCode ierr = temp.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PADirectForcing::temp_time_series(int i, int j, int N,
                                                 PetscReal *ts, PetscReal *values) {

  PetscReal *ptr;

  if (bc_period > 0.01) {
    // Recall that this method is called for each map-plane point during a
    // time-step. This if-condition is here to avoid calling my_mod() if the
    // user didn't ask for periodized climate.
    ts_mod.reserve(N);

    for (int k = 0; k < N; ++k)
      ts_mod[k] = my_mod(ts[k]);

    ptr = &ts_mod[0];
  } else {
    ptr = ts;
  }

  PetscErrorCode ierr = temp.interp(i, j, N, ptr, values); CHKERRQ(ierr);

  return 0;
}

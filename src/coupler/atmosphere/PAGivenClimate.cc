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

#include "PAGivenClimate.hh"
#include "IceGrid.hh"

PetscErrorCode PAGivenClimate::init(PISMVars &) {
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

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  // read time-independent data right away:
  if (temp.get_n_records() == 1 && mass_flux.get_n_records() == 1) {
    ierr = update(grid.time->current(), 0); CHKERRQ(ierr); // dt is irrelevant
  }

  return 0;
}

PetscErrorCode PAGivenClimate::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = mass_flux.at_time(t); CHKERRQ(ierr);
  ierr = temp.at_time(t); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAGivenClimate::mean_precipitation(IceModelVec2S &result) {
  PetscErrorCode ierr = mass_flux.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAGivenClimate::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr = temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAGivenClimate::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr = temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAGivenClimate::begin_pointwise_access() {
  PetscErrorCode ierr = temp.begin_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAGivenClimate::end_pointwise_access() {
  PetscErrorCode ierr = temp.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAGivenClimate::temp_time_series(int i, int j, int N,
                                                 PetscReal *ts, PetscReal *values) {

  PetscReal *ptr;

  if (bc_period > 0.01) {
    // Recall that this method is called for each map-plane point during a
    // time-step. This if-condition is here to avoid calling
    // grid.time->mod() if the user didn't ask for periodized climate.
    ts_mod.reserve(N);

    for (int k = 0; k < N; ++k)
      ts_mod[k] = grid.time->mod(ts[k] - bc_reference_time, bc_period);

    ptr = &ts_mod[0];
  } else {
    ptr = ts;
  }

  PetscErrorCode ierr = temp.interp(i, j, N, ptr, values); CHKERRQ(ierr);

  return 0;
}

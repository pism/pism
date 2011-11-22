// Copyright (C) 2011 PISM Authors
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

#include "PAAnomaly.hh"
#include "IceGrid.hh"

PetscErrorCode PAAnomaly::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (input_model != NULL) {
    ierr = input_model->init(vars); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 1, "input_model == NULL");
  }

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the -atmosphere ...,anomaly code...\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", ""); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr);

  ierr = temp.set_attrs("climate_forcing", "anomaly of the near-surface air temperature",
                        "Kelvin", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing", "anomaly of the ice-equivalent precipitation rate",
                             "m s-1", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  mass_flux.write_in_glaciological_units = true;

  ierr = verbPrintf(2, grid.com,
                    "    reading anomalies from %s ...\n",
                    filename.c_str()); CHKERRQ(ierr);

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAAnomaly::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  if (enable_time_averaging) {
    ierr = mass_flux.average(t, dt); CHKERRQ(ierr); 
    ierr = temp.average(t, 1.0); CHKERRQ(ierr); // compute the "mean annual" temperature
  } else {
    ierr = mass_flux.at_time(t); CHKERRQ(ierr);
    ierr = temp.at_time(t); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode PAAnomaly::mean_precip(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->mean_precip(result); CHKERRQ(ierr);

  return result.add(1.0, mass_flux);
}

PetscErrorCode PAAnomaly::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->mean_annual_temp(result); CHKERRQ(ierr);

  return result.add(1.0, temp);
}

PetscErrorCode PAAnomaly::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->temp_snapshot(result); CHKERRQ(ierr);

  return result.add(1.0, temp);
}


PetscErrorCode PAAnomaly::begin_pointwise_access() {
  PetscErrorCode ierr = input_model->begin_pointwise_access(); CHKERRQ(ierr);
  return temp.begin_access();
}

PetscErrorCode PAAnomaly::end_pointwise_access() {
  PetscErrorCode ierr = input_model->end_pointwise_access(); CHKERRQ(ierr);
  return temp.end_access();
}

PetscErrorCode PAAnomaly::temp_time_series(int i, int j, int N,
                                           PetscReal *ts, PetscReal *values) {
  PetscErrorCode ierr;

  // NB! the input_model uses un-periodized times.
  ierr = input_model->temp_time_series(i, j, N, ts, values); CHKERRQ(ierr);

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

  ts_values.reserve(N);
  ierr = temp.interp(i, j, N, ptr, &ts_values[0]); CHKERRQ(ierr);

  for (int k = 0; k < N; ++k)
    values[k] += ts_values[k];

  return 0;
}


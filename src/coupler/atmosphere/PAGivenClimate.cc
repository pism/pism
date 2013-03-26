// Copyright (C) 2011, 2012, 2013 PISM Authors
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

PAGivenClimate::PAGivenClimate(IceGrid &g, const NCConfigVariable &conf)
  : PGivenClimate<PAModifier,PISMAtmosphereModel>(g, conf, NULL)
{
  temp_name      = "air_temp";
  mass_flux_name = "precipitation";
  option_prefix  = "-atmosphere_given";

  // Cannot call allocate_PAGivenClimate() here, because some surface
  // models do not use atmosphere models *and* this is the default
  // atmosphere model.

}

PAGivenClimate::~PAGivenClimate() {
  // empty
}

PetscErrorCode PAGivenClimate::allocate_PAGivenClimate() {
  PetscErrorCode ierr;

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", ""); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr);

  ierr = temp.set_attrs("climate_forcing", "near-surface air temperature",
                        "Kelvin", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing", "ice-equivalent precipitation rate",
                       "m s-1", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAGivenClimate::init(PISMVars &) {
  PetscErrorCode ierr;

  t = dt = GSL_NAN;  // every re-init restarts the clock

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the atmosphere model reading near-surface air temperature\n"
                    "  and ice-equivalent precipitation from a file...\n"); CHKERRQ(ierr);

  if (temp.was_created() == false && mass_flux.was_created() == false) {
    ierr = allocate_PAGivenClimate(); CHKERRQ(ierr);
  }

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

  // compute mean precipitation
  ierr = mass_flux.average(t, dt, bc_period, bc_reference_time); CHKERRQ(ierr);

  // Average so that the mean_annual_temp() may be reported correctly (at least
  // in the "-surface pdd" case).
  ierr = temp.average(t, dt, bc_period, bc_reference_time); CHKERRQ(ierr);

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
  ierr = mass_flux.begin_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAGivenClimate::end_pointwise_access() {
  PetscErrorCode ierr = temp.end_access(); CHKERRQ(ierr);
  ierr = mass_flux.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAGivenClimate::temp_time_series(int i, int j, PetscReal *result) {
  PetscErrorCode ierr;

  ierr = temp.interp(i, j, &result[0]); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAGivenClimate::precip_time_series(int i, int j, PetscReal *result) {
  PetscErrorCode ierr;

  ierr = mass_flux.interp(i, j, &result[0]); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAGivenClimate::init_timeseries(PetscReal *ts, int N) {
  PetscErrorCode ierr;

  ierr = temp.init_interpolation(ts, N, bc_period, bc_reference_time); CHKERRQ(ierr);

  ierr = mass_flux.init_interpolation(ts, N, bc_period, bc_reference_time); CHKERRQ(ierr);

  m_ts_times.resize(N);
  
  return 0;
}


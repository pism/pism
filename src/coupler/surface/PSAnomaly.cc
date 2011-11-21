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

#include "PSAnomaly.hh"
#include "IceGrid.hh"

PetscErrorCode PSAnomaly::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (input_model != NULL) {
    ierr = input_model->init(vars); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com,
		    "* Initializing the '-surface ...,anomaly' modifier...\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", ""); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr);

  ierr = temp.set_attrs("climate_forcing",
                        "anomaly of the temperature of the ice at the ice surface but below firn processes",
                        "Kelvin", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing",
                             "anomaly of the ice-equivalent surface mass balance (accumulation/ablation) rate",
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

PetscErrorCode PSAnomaly::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  if (enable_time_averaging) {
    ierr = mass_flux.average(t, dt); CHKERRQ(ierr);
    ierr = temp.average(t, dt); CHKERRQ(ierr);
  } else {
    ierr = mass_flux.at_time(t); CHKERRQ(ierr);
    ierr = temp.at_time(t); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PSAnomaly::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->ice_surface_mass_flux(result); CHKERRQ(ierr);
  ierr = result.add(1.0, mass_flux); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSAnomaly::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->ice_surface_temperature(result); CHKERRQ(ierr);
  ierr = result.add(1.0, temp); CHKERRQ(ierr);

  return 0;
}

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

#include "PSStuffAsAnomaly.hh"

PetscErrorCode PSStuffAsAnomaly::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (input_model != NULL) {
    ierr = input_model->init(vars); CHKERRQ(ierr);
  }

  // create mean annual ice equivalent snow precipitation rate (before melt, and not including rain)
  ierr = mass_flux.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_state",
                             "ice-equivalent surface mass balance (accumulation/ablation) rate",
                             "m s-1",
                             "land_ice_surface_specific_mass_balance"); // CF standard_name
			CHKERRQ(ierr);
  ierr = mass_flux.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  mass_flux.write_in_glaciological_units = true;

  ierr = temp.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = temp.set_attrs("climate_state",
			"ice temperature at the ice surface",
			"K",
			""); CHKERRQ(ierr);

  // create special variables
  ierr = mass_flux_0.create(grid, "mass_flux_0", false); CHKERRQ(ierr);
  ierr = mass_flux_0.set_attrs("internal", "surface mass flux at the beginning of a run",
                               "m s-1", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);

  ierr = mass_flux_input.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = mass_flux_input.set_attrs("model_state", "surface mass flux to apply anomalies to",
                                    "m s-1", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);

  ierr = temp_0.create(grid, "temp_0", false); CHKERRQ(ierr);
  ierr = temp_0.set_attrs("internal", "ice-surface temperature and the beginning of a run", "K",
                          ""); CHKERRQ(ierr);

  ierr = temp_input.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = temp_input.set_attrs("model_state", "ice-surface temperature to apply anomalies to", "K",
                              ""); CHKERRQ(ierr);

  ierr = find_pism_input(input_file, regrid, start); CHKERRQ(ierr);

  if (regrid) {
    ierr = mass_flux_input.regrid(input_file.c_str(), true); CHKERRQ(ierr); // fails if not found!
    ierr = temp_input.regrid(input_file.c_str(), true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = mass_flux_input.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
    ierr = temp_input.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }

  // get the mass balance and the temperature at the beginning of the run:
  if (input_model != NULL) {
    ierr = input_model->update(grid.time->start(), 0); CHKERRQ(ierr);
    ierr = input_model->ice_surface_mass_flux(mass_flux_0); CHKERRQ(ierr);
    ierr = input_model->ice_surface_temperature(temp_0); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PSStuffAsAnomaly::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  if (input_model != NULL) {
    ierr = input_model->ice_surface_temperature(result); CHKERRQ(ierr);
  }

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = mass_flux_0.begin_access(); CHKERRQ(ierr);
  ierr = mass_flux_input.begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i, j) = result(i, j) - mass_flux_0(i, j) + mass_flux_input(i, j);
    }
  }

  ierr = mass_flux_input.end_access(); CHKERRQ(ierr);
  ierr = mass_flux_0.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSStuffAsAnomaly::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;

  if (input_model != NULL) {
    ierr = input_model->ice_surface_temperature(result); CHKERRQ(ierr);
  }

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = temp_0.begin_access(); CHKERRQ(ierr);
  ierr = temp_input.begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i, j) = result(i, j) - temp_0(i, j) + temp_input(i, j);
    }
  }

  ierr = temp_input.end_access(); CHKERRQ(ierr);
  ierr = temp_0.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

void PSStuffAsAnomaly::add_vars_to_output(string keyword, set<string> &result) {
  if (keyword != "small") {
    result.insert("acab");
    result.insert("artm");
  }
}

PetscErrorCode PSStuffAsAnomaly::define_variables(set<string> vars, const NCTool &nc, nc_type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "artm")) {
    ierr = temp.define(nc, nctype); CHKERRQ(ierr); 
  }

  if (set_contains(vars, "acab")) {
    ierr = mass_flux.define(nc, nctype); CHKERRQ(ierr);
  }

  // ensure that no one overwrites these two
  vars.erase("artm");
  vars.erase("acab");

  ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSStuffAsAnomaly::write_variables(set<string> vars, string fname) {
  PetscErrorCode ierr;

  if (set_contains(vars, "artm")) {
    ierr = temp.write(fname); CHKERRQ(ierr);
  }

  if (set_contains(vars, "acab")) {
    ierr = mass_flux.write(fname); CHKERRQ(ierr);
  }

  // ensure that no one overwrites these two
  vars.erase("artm");
  vars.erase("acab");

  ierr = input_model->write_variables(vars, fname); CHKERRQ(ierr);

  return 0;
}


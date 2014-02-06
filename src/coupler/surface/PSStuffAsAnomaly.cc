// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include "IceGrid.hh"
#include "PISMTime.hh"

PSStuffAsAnomaly::PSStuffAsAnomaly(IceGrid &g, const PISMConfig &conf, PISMSurfaceModel *input)
    : PSModifier(g, conf, input) {
  PetscErrorCode ierr = allocate_PSStuffAsAnomaly(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

PSStuffAsAnomaly::~PSStuffAsAnomaly() {
  // empty
}

PetscErrorCode PSStuffAsAnomaly::allocate_PSStuffAsAnomaly() {
  PetscErrorCode ierr;

  ierr = mass_flux.create(grid, "climatic_mass_balance", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_state",
                             "surface mass balance (accumulation/ablation) rate",
                             "kg m-2 s-1",
                             "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);
  ierr = mass_flux.set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);
  mass_flux.write_in_glaciological_units = true;

  ierr = temp.create(grid, "ice_surface_temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = temp.set_attrs("climate_state", "ice temperature at the ice surface",
			"K", ""); CHKERRQ(ierr);

  // create special variables
  ierr = mass_flux_0.create(grid, "mass_flux_0", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = mass_flux_0.set_attrs("internal", "surface mass flux at the beginning of a run",
                               "kg m-2 s-1", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);

  ierr = mass_flux_input.create(grid, "climatic_mass_balance", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = mass_flux_input.set_attrs("model_state", "surface mass flux to apply anomalies to",
                                   "kg m-2 s-1", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);

  ierr = temp_0.create(grid, "ice_surface_temp_0", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = temp_0.set_attrs("internal", "ice-surface temperature and the beginning of a run", "K",
                          ""); CHKERRQ(ierr);

  ierr = temp_input.create(grid, "ice_surface_temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = temp_input.set_attrs("model_state", "ice-surface temperature to apply anomalies to",
                              "K", ""); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSStuffAsAnomaly::init(PISMVars &vars) {
  PetscErrorCode ierr;
  std::string input_file;
  bool do_regrid = false;
  int start = 0;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  if (input_model != NULL) {
    ierr = input_model->init(vars); CHKERRQ(ierr);
  }

  ierr = find_pism_input(input_file, do_regrid, start); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
		    "* Initializing the 'turn_into_anomaly' modifier\n"
                    "  (it applies climate data as anomalies relative to 'ice_surface_temp' and 'climatic_mass_balance'\n"
                    "  read from '%s'.\n", input_file.c_str()); CHKERRQ(ierr);

  if (do_regrid) {
    ierr = mass_flux_input.regrid(input_file, CRITICAL); CHKERRQ(ierr); // fails if not found!
    ierr = temp_input.regrid(input_file, CRITICAL); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = mass_flux_input.read(input_file, start); CHKERRQ(ierr); // fails if not found!
    ierr = temp_input.read(input_file, start); CHKERRQ(ierr); // fails if not found!
  }

  return 0;
}

PetscErrorCode PSStuffAsAnomaly::update(double my_t, double my_dt) {
  PetscErrorCode ierr;

  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12))
    return 0;

  m_t  = my_t;
  m_dt = my_dt;

  if (input_model != NULL) {
    ierr = input_model->update(m_t, m_dt); CHKERRQ(ierr);
    ierr = input_model->ice_surface_temperature(temp); CHKERRQ(ierr);
    ierr = input_model->ice_surface_mass_flux(mass_flux); CHKERRQ(ierr);

    // if we are at the beginning of the run...
    if (m_t < grid.time->start() + 1) { // this is goofy, but time-steps are
                                      // usually longer than 1 second, so it
                                      // should work
      ierr = temp.copy_to(temp_0); CHKERRQ(ierr);
      ierr = mass_flux.copy_to(mass_flux_0); CHKERRQ(ierr);
    }
  }

  ierr = mass_flux.begin_access(); CHKERRQ(ierr);
  ierr = mass_flux_0.begin_access(); CHKERRQ(ierr);
  ierr = mass_flux_input.begin_access(); CHKERRQ(ierr);

  ierr = temp.begin_access(); CHKERRQ(ierr);
  ierr = temp_0.begin_access(); CHKERRQ(ierr);
  ierr = temp_input.begin_access(); CHKERRQ(ierr);

  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      mass_flux(i, j) = mass_flux(i, j) - mass_flux_0(i, j) + mass_flux_input(i, j);
      temp(i, j) = temp(i, j) - temp_0(i, j) + temp_input(i, j);
    }
  }

  ierr = temp_input.end_access(); CHKERRQ(ierr);
  ierr = temp_0.end_access(); CHKERRQ(ierr);
  ierr = temp.end_access(); CHKERRQ(ierr);

  ierr = mass_flux_input.end_access(); CHKERRQ(ierr);
  ierr = mass_flux_0.end_access(); CHKERRQ(ierr);
  ierr = mass_flux.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSStuffAsAnomaly::ice_surface_mass_flux(IceModelVec2S &result) {
  return mass_flux.copy_to(result);
}

PetscErrorCode PSStuffAsAnomaly::ice_surface_temperature(IceModelVec2S &result) {
  return temp.copy_to(result);
}

void PSStuffAsAnomaly::add_vars_to_output(std::string keyword, std::set<std::string> &result) {
  if (input_model != NULL) {
    input_model->add_vars_to_output(keyword, result);
  }

  result.insert("ice_surface_temp");
  result.insert("climatic_mass_balance");
}

PetscErrorCode PSStuffAsAnomaly::define_variables(std::set<std::string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    ierr = temp.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = mass_flux.define(nc, nctype); CHKERRQ(ierr);
  }

  // ensure that no one overwrites these two
  vars.erase("ice_surface_temp");
  vars.erase("climatic_mass_balance");

  if (input_model != NULL) {
    ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PSStuffAsAnomaly::write_variables(std::set<std::string> vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    ierr = temp.write(nc); CHKERRQ(ierr);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = mass_flux.write(nc); CHKERRQ(ierr);
  }

  // ensure that no one overwrites these two
  vars.erase("ice_surface_temp");
  vars.erase("climatic_mass_balance");

  if (input_model != NULL) {
    ierr = input_model->write_variables(vars, nc); CHKERRQ(ierr);
  }

  return 0;
}


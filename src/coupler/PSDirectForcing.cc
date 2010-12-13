// Copyright (C) 2010 Constantine Khroulev
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

#include "PSDirectForcing.hh"

inline bool set_contains(set<string> S, string name) {
  return (S.find(name) != S.end());
}

PetscReal PSDirectForcing::my_mod(PetscReal input) {
  if (bc_period < 0.01) return input;

  return input - floor(input / bc_period) * bc_period;
}


PetscErrorCode PSDirectForcing::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;
  string filename;
  bool bc_file_set, bc_period_set;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the surface model reading temperature at the top of the ice\n"
                    "  and ice surface mass flux from a file...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Direct forcing options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-bc_file", "Specifies a file with top-surface boundary conditions",
                             filename, bc_file_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-bc_period", "Specifies the length of the climate data period",
                           bc_period, bc_period_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (bc_file_set == false) {
    PetscPrintf(grid.com, "PISM ERROR: option -bc_file is required.\n");
    PetscEnd();
  }

  if (bc_period_set == false) {
    bc_period = 0;
  }

  ierr = verbPrintf(2,grid.com,
                    "    reading boundary conditions from %s ...\n",
                    filename.c_str()); CHKERRQ(ierr);

  temperature.set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
  ierr = temperature.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = temperature.set_attrs("climate_forcing",
                               "temperature of the ice at the ice surface but below firn processes",
                               "Kelvin", ""); CHKERRQ(ierr);
  ierr = temperature.init(filename); CHKERRQ(ierr);

  mass_flux.set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
  ierr = mass_flux.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing",
                             "ice-equivalent surface mass balance (accumulation/ablation) rate",
                             "m s-1", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PSDirectForcing::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  // "Periodize" the climate:
  t_years = my_mod(t_years);

  if ((fabs(t_years - t) < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  t  = t_years;
  dt = dt_years;

  ierr = temperature.update(t, dt); CHKERRQ(ierr);
  ierr = mass_flux.update(t, dt); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSDirectForcing::max_timestep(PetscReal t_years, PetscReal &dt_years) {
  PetscReal max_dt = -1;

  // "Periodize" the climate:
  t_years = my_mod(t_years);

  max_dt = temperature.max_timestep(t_years);

  if (dt_years > 0) {
    if (max_dt > 0)
      dt_years = PetscMin(max_dt, dt_years);
  }
  else dt_years = max_dt;

  max_dt = mass_flux.max_timestep(t_years);

  if (dt_years > 0) {
    if (max_dt > 0)
      dt_years = PetscMin(max_dt, dt_years);
  }
  else dt_years = max_dt;

  return 0;
}

PetscErrorCode PSDirectForcing::write_fields(set<string> vars,
                                             PetscReal t_years, PetscReal dt_years,
                                             string filename) {
  PetscErrorCode ierr;

  // "Periodize" the climate:
  t_years = my_mod(t_years);

  if (set_contains(vars, "artm")) {
    ierr = temperature.interp(t_years + 0.5*dt_years); CHKERRQ(ierr); 
    ierr = temperature.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "acab")) {
    ierr = mass_flux.interp(t_years + 0.5*dt_years); CHKERRQ(ierr); 
    ierr = mass_flux.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PSDirectForcing::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
                                                        string filename) {
  PetscErrorCode ierr;

  // "Periodize" the climate:
  t_years = my_mod(t_years);

  ierr = temperature.interp(t_years + 0.5*dt_years); CHKERRQ(ierr); 
  ierr = temperature.write(filename.c_str()); CHKERRQ(ierr);

  ierr = mass_flux.interp(t_years + 0.5*dt_years); CHKERRQ(ierr); 
  ierr = mass_flux.write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}



PetscErrorCode PSDirectForcing::ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
                                                      IceModelVec2S &result) {
  PetscErrorCode ierr;

  // "Periodize" the climate:
  t_years = my_mod(t_years);

  ierr = update(t_years, dt_years); CHKERRQ(ierr); 

  ierr = mass_flux.interp(t_years + 0.5*dt_years); CHKERRQ(ierr);

  ierr = mass_flux.copy_to(result); CHKERRQ(ierr); 

  return 0;
}


PetscErrorCode PSDirectForcing::ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
                                                        IceModelVec2S &result) {
  PetscErrorCode ierr;

  // "Periodize" the climate:
  t_years = my_mod(t_years);

  ierr = update(t_years, dt_years); CHKERRQ(ierr); 

  ierr = temperature.interp(t_years + 0.5*dt_years); CHKERRQ(ierr);

  ierr = temperature.copy_to(result); CHKERRQ(ierr); 

  return 0;
}


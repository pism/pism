// Copyright (C) 2011 Andy Aschwanden and Constantine Khroulev
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

#include "PADirectForcing.hh"

PADirectForcing::PADirectForcing(IceGrid &g, const NCConfigVariable &conf)
  : PISMAtmosphereModel(g, conf) {
}

PADirectForcing::~PADirectForcing() {
}

PetscErrorCode PADirectForcing::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;
  PetscTruth bc_file_set;
  char bc_file[PETSC_MAX_PATH_LEN];

  ierr = verbPrintf(2, grid.com,
     "* Initializing air temperature and precipitation forcing...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", 
                           "Air temperature and precipitation forcing", ""); CHKERRQ(ierr);

  ierr = PetscOptionsString("-bc_file",
			    "Specifies the air temperature and preciptation forcing file",
			    "", "",
			    bc_file, PETSC_MAX_PATH_LEN, &bc_file_set);
			    CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (bc_file_set == false) {
    PetscPrintf(grid.com, "PISM ERROR: option -bc_file is required.\n");
    PISMEnd();
  }
  
  ierr = verbPrintf(2,grid.com,
                    "    reading air temperature from %s ...\n",
                    bc_file); CHKERRQ(ierr);
  temp.set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
  ierr = temp.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = temp.set_attrs("climate_forcing",
                                 "near-surface temperature",
                                 "Kelvin", ""); CHKERRQ(ierr);
  ierr = temp.init(bc_file); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
                    "    reading ice-equivalent precipitation rate from %s ...\n",
                    bc_file); CHKERRQ(ierr);
  
  precip.set_n_records((unsigned int) config.get("climate_forcing_buffer_size")); 
  ierr = precip.create(grid, "precip", false); CHKERRQ(ierr);
  ierr = precip.set_attrs("climate_forcing",
                                   "ice-equivalent precipitation rate",
                                   "m s-1", ""); CHKERRQ(ierr);
  ierr = precip.init(bc_file); CHKERRQ(ierr);
  ierr = precip.set_glaciological_units("m year-1");
  precip.write_in_glaciological_units = true;

  return 0;
}

PetscErrorCode PADirectForcing::max_timestep(PetscReal t_years,
                                             PetscReal &dt_years) {
  PetscReal max_dt = -1;
  
  max_dt = temp.max_timestep(t_years);

  if (dt_years > 0) {
    if (max_dt > 0)
      dt_years = PetscMin(max_dt, dt_years);
  }
  else dt_years = max_dt;

  dt_years = precip.max_timestep(t_years);

  if (dt_years > 0) {
    if (max_dt > 0)
      dt_years = PetscMin(max_dt, dt_years);
  }
  else dt_years = max_dt;

  return 0;
}

void PADirectForcing::add_vars_to_output(string keyword, set<string> &result) {
  if (keyword == "big") {
    result.insert("artm");
    result.insert("precip");
  }
}

PetscErrorCode PADirectForcing::define_variables(set<string> vars, const NCTool &nc, nc_type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "artm")) {
    ierr = temp.define(nc, nctype); CHKERRQ(ierr);
  }
  if (set_contains(vars, "precip")) {
    ierr = precip.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode PADirectForcing::write_variables(set<string> vars,  string filename) {
  PetscErrorCode ierr;
  double T = t + 0.5 * dt;

  if (set_contains(vars, "artm")) {
    ierr = temp.update(T, 0); CHKERRQ(ierr);
    ierr = temp.interp(T); CHKERRQ(ierr);
    ierr = temp.write(filename.c_str()); CHKERRQ(ierr);
    vars.erase("temp");
  }

  if (set_contains(vars, "precip")) {
    ierr = precip.update(T, 0); CHKERRQ(ierr);
    ierr = precip.interp(T); CHKERRQ(ierr);
    ierr = precip.write(filename.c_str()); CHKERRQ(ierr);
    vars.erase("precip");
  }

  return 0;
}

PetscErrorCode PADirectForcing::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  if ((fabs(t_years - t) < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  t  = t_years;
  dt = dt_years;

  ierr = temp.update(t_years, dt_years); CHKERRQ(ierr);

  ierr = precip.update(t_years, dt_years); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PADirectForcing::mean_precip(IceModelVec2S &result) {
  PetscErrorCode ierr;

  string history = "added average over time-step of precipitation\n" + result.string_attr("history");

  double anomaly;
  ierr = precip.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = precip.average(i, j, t, dt, anomaly); CHKERRQ(ierr);
      result(i,j) = anomaly;
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = precip.end_access(); CHKERRQ(ierr);

  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PADirectForcing::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr;
  
  string history = "added annual average of temperature\n"
    + result.string_attr("history");
    
  // average over the year starting at t_years
  PetscScalar av_ij;
  ierr = temp.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = temp.average(i,j,t,1.0,av_ij); CHKERRQ(ierr);
      result(i,j) = av_ij;
    }
  }
  ierr = temp.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PADirectForcing::begin_pointwise_access() {
  PetscErrorCode ierr;
  ierr = temp.begin_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PADirectForcing::end_pointwise_access() {
  PetscErrorCode ierr;
  ierr = temp.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PADirectForcing::temp_time_series(int i, int j, int N,
                                                 PetscReal *ts, PetscReal *values) {
  PetscErrorCode ierr;
  PetscScalar *tmp = new PetscScalar[N];

  ierr = temp.interp(i, j, N, ts, tmp); CHKERRQ(ierr);

  for (PetscInt k = 0; k < N; k++)
    values[k] = tmp[k];

  delete tmp;

  return 0;
}

PetscErrorCode PADirectForcing::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr;
  double T = t + 0.5 * dt;

  string history = "added temperatures at midpoint of time-step\n"
    + result.string_attr("history");
  
  ierr = temp.interp(T); CHKERRQ(ierr);
  ierr = temp.copy_to(result);
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

// Copyright (C) 2008-2011 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

// Implementation of forcing using spatially-variable time-dependent air
// temperature and precipitaiton anomalies (-atmosphere ...,forcing).

#include "PISMAtmosphere.hh"

PAForcing::PAForcing(IceGrid &g, const NCConfigVariable &conf, PISMAtmosphereModel *input)
  : PAModifier(g, conf, input) {
  temp_anomaly = NULL;
  precip_anomaly = NULL;
}

PAForcing::~PAForcing() {
  delete temp_anomaly;
  delete precip_anomaly;
}

PetscErrorCode PAForcing::init(PISMVars &vars) {
  PetscErrorCode ierr;
  PetscTruth temp_anomaly_set, precip_anomaly_set;
  char temp_anomalies_file[PETSC_MAX_PATH_LEN],
    precip_anomalies_file[PETSC_MAX_PATH_LEN];

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "* Initializing air temperature and precipitation forcing...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", 
                           "Air temperature and precipitation forcing", ""); CHKERRQ(ierr);

  ierr = PetscOptionsString("-anomaly_temp",
			    "Specifies the air temperature anomalies file",
			    "", "",
			    temp_anomalies_file, PETSC_MAX_PATH_LEN, &temp_anomaly_set);
			    CHKERRQ(ierr);
  ierr = PetscOptionsString("-anomaly_precip", 
			    "Specifies the precipitation anomalies file",
			    "", "",
			    precip_anomalies_file, PETSC_MAX_PATH_LEN, &precip_anomaly_set);
			    CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (! (temp_anomaly_set || precip_anomaly_set) ) {
    ierr = verbPrintf(2, grid.com, "  NOTE: Forcing is inactive...\n"); CHKERRQ(ierr);
  }

  // check on whether we should read temperature anomalies
  if (temp_anomaly_set) {

    ierr = verbPrintf(2,grid.com,
		      "    reading air temperature anomalies from %s ...\n",
		      temp_anomalies_file); CHKERRQ(ierr);
    temp_anomaly = new IceModelVec2T;
    temp_anomaly->set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
    ierr = temp_anomaly->create(grid, "temp_anomaly", false); CHKERRQ(ierr);
    ierr = temp_anomaly->set_attrs("climate_forcing",
                                   "near-surface temperature anomalies",
				   "Kelvin", ""); CHKERRQ(ierr);
    ierr = temp_anomaly->init(temp_anomalies_file); CHKERRQ(ierr);
  }

  // check on whether we should read precipitation anomalies
  if (precip_anomaly_set) {
    ierr = verbPrintf(2,grid.com,
		      "    reading ice-equivalent precipitation rate anomalies from %s ...\n",
		      precip_anomalies_file); CHKERRQ(ierr);

    precip_anomaly = new IceModelVec2T;
    precip_anomaly->set_n_records((unsigned int) config.get("climate_forcing_buffer_size")); 
    ierr = precip_anomaly->create(grid, "precip_anomaly", false); CHKERRQ(ierr);
    ierr = precip_anomaly->set_attrs("climate_forcing",
					 "anomalies of ice-equivalent precipitation rate",
					 "m s-1", ""); CHKERRQ(ierr);
    ierr = precip_anomaly->init(precip_anomalies_file); CHKERRQ(ierr);
    ierr = precip_anomaly->set_glaciological_units("m year-1");
    precip_anomaly->write_in_glaciological_units = true;
  }

  airtemp_var.init_2d("airtemp_plus_forcing", grid);
  airtemp_var.set_string("pism_intent", "diagnostic");
  airtemp_var.set_string("long_name",
                         "snapshot of the near-surface air temperature (including forcing)");
  ierr = airtemp_var.set_units("K"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAForcing::max_timestep(PetscReal my_t,
				       PetscReal &my_dt, bool &restrict) {
  PetscErrorCode ierr;
  PetscReal max_dt = -1;
  
  ierr = input_model->max_timestep(my_t, my_dt, restrict); CHKERRQ(ierr);

  if (temp_anomaly != NULL) {
    max_dt = temp_anomaly->max_timestep(my_t);

    if (my_dt > 0) {
      if (max_dt > 0)
	my_dt = PetscMin(max_dt, my_dt);
    }
    else my_dt = max_dt;
  }

  if (precip_anomaly != NULL) {
    my_dt = precip_anomaly->max_timestep(my_t);

    if (my_dt > 0) {
      if (max_dt > 0)
	my_dt = PetscMin(max_dt, my_dt);
    }
    else my_dt = max_dt;
  }

  if (my_dt > 0)
    restrict = true;
  else
    restrict = false;

  return 0;
}

void PAForcing::add_vars_to_output(string keyword, set<string> &result) {
  if (keyword == "big") {
    result.insert("airtemp_plus_forcing");

    if (temp_anomaly != NULL)
      result.insert("temp_anomaly");

    if (precip_anomaly != NULL)
      result.insert("precip_anomaly");
  }

  input_model->add_vars_to_output(keyword, result);
}

PetscErrorCode PAForcing::define_variables(set<string> vars, const NCTool &nc, nc_type nctype) {
  PetscErrorCode ierr;
  int varid;

  ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);

  if (set_contains(vars, "airtemp_plus_forcing")) {
    ierr = airtemp_var.define(nc, varid, nctype, false); CHKERRQ(ierr);
  }

  if (temp_anomaly != NULL) {
    if (set_contains(vars, "temp_anomaly")) {
      ierr = temp_anomaly->define(nc, nctype); CHKERRQ(ierr);
    }
  }

  if (precip_anomaly != NULL) {
    if (set_contains(vars, "precip_anomaly")) {
      ierr = precip_anomaly->define(nc, nctype); CHKERRQ(ierr);
    }
  }

  return 0;
}


PetscErrorCode PAForcing::write_variables(set<string> vars,  string filename) {
  PetscErrorCode ierr;
  double T = t + 0.5 * dt;

  if (set_contains(vars, "airtemp_plus_forcing")) {
    IceModelVec2S airtemp;
    ierr = airtemp.create(grid, "airtemp_plus_forcing", false); CHKERRQ(ierr);
    ierr = airtemp.set_metadata(airtemp_var, 0); CHKERRQ(ierr);

    ierr = temp_snapshot(airtemp); CHKERRQ(ierr);

    ierr = airtemp.write(filename.c_str()); CHKERRQ(ierr);
    vars.erase("airtemp_plus_forcing");
  }

  if (temp_anomaly != NULL) {
    if (set_contains(vars, "temp_anomaly")) {
      ierr = temp_anomaly->update(T, 0); CHKERRQ(ierr);
      ierr = temp_anomaly->interp(T); CHKERRQ(ierr);
      ierr = temp_anomaly->write(filename.c_str()); CHKERRQ(ierr);
      vars.erase("temp_anomaly");
    }
  }

  if (precip_anomaly != NULL) {
    if (set_contains(vars, "precip_anomaly")) {
      ierr = precip_anomaly->update(T, 0); CHKERRQ(ierr);
      ierr = precip_anomaly->interp(T); CHKERRQ(ierr);
      ierr = precip_anomaly->write(filename.c_str()); CHKERRQ(ierr);
      vars.erase("precip_anomaly");
    }

  }

  ierr = input_model->write_variables(vars, filename); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAForcing::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  if ((fabs(my_t - t) < 1e-12) &&
      (fabs(my_dt - dt) < 1e-12))
    return 0;

  t  = my_t;
  dt = my_dt;

  ierr = input_model->update(my_t, my_dt); CHKERRQ(ierr);

  if (temp_anomaly != NULL) {
    ierr = temp_anomaly->update(my_t, my_dt); CHKERRQ(ierr);
  }

  if (precip_anomaly != NULL) {
    ierr = precip_anomaly->update(my_t, my_dt); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::mean_precip(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->mean_precip(result); CHKERRQ(ierr);

  if (precip_anomaly != NULL) {
    string history = "added average over time-step of precipitation anomalies\n" + result.string_attr("history");

    double anomaly;

    ierr = precip_anomaly->begin_access(); CHKERRQ(ierr);
    ierr = result.begin_access(); CHKERRQ(ierr);
    for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
	ierr = precip_anomaly->average(i, j, t, dt, anomaly); CHKERRQ(ierr);
	result(i,j) += anomaly;
      }
    }
    ierr = result.end_access(); CHKERRQ(ierr);
    ierr = precip_anomaly->end_access(); CHKERRQ(ierr);

    ierr = result.set_attr("history", history); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->mean_annual_temp(result); CHKERRQ(ierr);
  
  if (temp_anomaly != NULL) {
    string history = "added annual average of temperature anomalies\n"
                     + result.string_attr("history");
    // if the temperature anomaly has sub-annual cycle we need to average it
    // FIXME: if the anomaly is provided with a flag indicating it is already
    //        mean annual, then we could save some work
    
    // average over the year starting at t
    ierr = temp_anomaly->update(t, 1.0); CHKERRQ(ierr);
    PetscScalar av_ij;
    ierr = temp_anomaly->begin_access(); CHKERRQ(ierr);
    ierr = result.begin_access(); CHKERRQ(ierr);
    for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
        ierr = temp_anomaly->average(i,j,t,1.0,av_ij); CHKERRQ(ierr);
        result(i,j) += av_ij;
      }
    }
    ierr = temp_anomaly->end_access(); CHKERRQ(ierr);
    ierr = result.end_access(); CHKERRQ(ierr);

    ierr = result.set_attr("history", history); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::begin_pointwise_access() {
  PetscErrorCode ierr;

  ierr = input_model->begin_pointwise_access(); CHKERRQ(ierr);

  if (temp_anomaly != NULL) {
    ierr = temp_anomaly->begin_access(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::end_pointwise_access() {
  PetscErrorCode ierr;

  ierr = input_model->end_pointwise_access(); CHKERRQ(ierr);

  if (temp_anomaly != NULL) {
    ierr = temp_anomaly->end_access(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::temp_time_series(int i, int j, int N,
					   PetscReal *ts, PetscReal *values) {
  PetscErrorCode ierr;
  
  ierr = input_model->temp_time_series(i, j, N, ts, values); CHKERRQ(ierr);

  if (temp_anomaly != NULL) {
    PetscScalar *tmp = new PetscScalar[N];

    ierr = temp_anomaly->interp(i, j, N, ts, tmp); CHKERRQ(ierr);

    for (PetscInt k = 0; k < N; k++)
      values[k] += tmp[k];

    delete tmp;
  }

  return 0;
}

PetscErrorCode PAForcing::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->temp_snapshot(result); CHKERRQ(ierr);

  if (temp_anomaly != NULL) {
    string history = "added temperature anomalies at midpoint of time-step\n"
                     + result.string_attr("history");

    ierr = temp_anomaly->update(t, dt); CHKERRQ(ierr);
    ierr = temp_anomaly->interp(t); CHKERRQ(ierr);
    ierr = result.add(1.0, *temp_anomaly); CHKERRQ(ierr);

    ierr = result.set_attr("history", history); CHKERRQ(ierr);
  }

  return 0;
}

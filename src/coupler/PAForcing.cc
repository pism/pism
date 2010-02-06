// Copyright (C) 2008-2010 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#include "PISMAtmosphere.hh"

void PAModifier::attach_input(PISMAtmosphereModel *input) {
  if (input_model != NULL)
    delete input_model;
  input_model = input;
}

PAForcing::PAForcing(IceGrid &g, const NCConfigVariable &conf, PISMVars &vars)
  : PAModifier(g, conf, vars) {
  dTforcing = NULL;
  delta_T = NULL;
  temp_ma_anomaly = NULL;
  snowprecip_anomaly = NULL;
  paleo_precipitation_correction = false;
}

PAForcing::~PAForcing() {
  delete dTforcing;
  delete delta_T;
  delete temp_ma_anomaly;
  delete snowprecip_anomaly;
}

PetscErrorCode PAForcing::init() {
  PetscErrorCode ierr;

  ierr = input_model->init(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "* Initializing air temperature and precipitation forcing...\n"); CHKERRQ(ierr);

  // check on whether we should read mean annual temperature anomalies
  PetscTruth optSet;
  char anomalies_file[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-anomaly_temp_ma", 
			       anomalies_file, PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
  if (optSet) {
    // stop if -dTforcing is set:
    PetscTruth dTforcing_set;
    ierr = check_option("-dTforcing", dTforcing_set); CHKERRQ(ierr);
    if (dTforcing_set) {
      ierr = PetscPrintf(grid.com, "PISM ERROR: option -anomaly_temp_ma is incompatible with -dTforcing.\n");
      PetscEnd();
    }

    ierr = verbPrintf(2,grid.com,
		      "    reading air temperature anomalies from %s ...\n",
		      anomalies_file); CHKERRQ(ierr);

    temp_ma_anomaly = new IceModelVec2T;
    temp_ma_anomaly->set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
    ierr = temp_ma_anomaly->create(grid, "temp_ma_anomaly", false); CHKERRQ(ierr);
    ierr = temp_ma_anomaly->set_attrs("climate_forcing", "mean annual near-surface temperature anomalies",
				      "Kelvin", ""); CHKERRQ(ierr);
    ierr = temp_ma_anomaly->init(anomalies_file); CHKERRQ(ierr);
  }

  // check on whether we should read snow precipitation anomalies
  ierr = PetscOptionsGetString(PETSC_NULL, "-anomaly_precip", 
			       anomalies_file, PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
  if (optSet) {
    ierr = verbPrintf(2,grid.com,
		      "    reading ice-equivalent snow precipitation rate anomalies from %s ...\n",
		      anomalies_file); CHKERRQ(ierr);

    snowprecip_anomaly = new IceModelVec2T;
    snowprecip_anomaly->set_n_records((unsigned int) config.get("climate_forcing_buffer_size")); 
    ierr = snowprecip_anomaly->create(grid, "snowprecip_anomaly", false); CHKERRQ(ierr);
    ierr = snowprecip_anomaly->set_attrs("climate_forcing",
					 "anomalies of ice-equivalent snow precipitation rate",
					 "m s-1", ""); CHKERRQ(ierr);
    ierr = snowprecip_anomaly->init(anomalies_file); CHKERRQ(ierr);
    ierr = snowprecip_anomaly->set_glaciological_units("m year-1");
    snowprecip_anomaly->write_in_glaciological_units = true;
  }

  // check user option -dTforcing for a surface temperature forcing data set
  char dTfile[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-dTforcing", dTfile,
                               PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
  if (optSet == PETSC_TRUE) {
    ierr = check_option("-paleo_precip", optSet); CHKERRQ(ierr);

    if (optSet) {
      ierr = verbPrintf(2,grid.com,
			"    using the paleo-precipitation correction...\n",
			anomalies_file); CHKERRQ(ierr);
      paleo_precipitation_correction = true;
    }

    dTforcing = new Timeseries(grid.com, grid.rank, "delta_T", "t");
    ierr = dTforcing->set_units("Celsius", ""); CHKERRQ(ierr);
    ierr = dTforcing->set_dimension_units("years", ""); CHKERRQ(ierr);
    ierr = dTforcing->set_attr("long_name", "near-surface air temperature offsets"); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com, 
		      "  reading delta T data from forcing file %s...\n", dTfile);  CHKERRQ(ierr);
	 
    ierr = dTforcing->read(dTfile); CHKERRQ(ierr);

    delta_T = new DiagnosticTimeseries(grid.com, grid.rank, "delta_T", "t");
    ierr = delta_T->set_units("Celsius", ""); CHKERRQ(ierr);
    ierr = delta_T->set_dimension_units("years", ""); CHKERRQ(ierr);
    ierr = delta_T->set_attr("long_name", "near-surface air temperature offset"); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::max_timestep(PetscReal t_years,
				       PetscReal &dt_years) {
  PetscErrorCode ierr;
  PetscReal max_dt = -1;
  
  ierr = input_model->max_timestep(t_years, dt_years); CHKERRQ(ierr);

  if (temp_ma_anomaly != NULL) {
    max_dt = temp_ma_anomaly->max_timestep(t_years);

    if (dt_years > 0) {
      if (max_dt > 0)
	dt_years = PetscMin(max_dt, dt_years);
    }
    else dt_years = max_dt;
  }

  if (snowprecip_anomaly != NULL) {
    dt_years = snowprecip_anomaly->max_timestep(t_years);

    if (dt_years > 0) {
      if (max_dt > 0)
	dt_years = PetscMin(max_dt, dt_years);
    }
    else dt_years = max_dt;
  }

  return 0;
}

PetscErrorCode PAForcing::write_input_fields(PetscReal t_years, PetscReal dt_years,
					     string filename) {
  PetscErrorCode ierr;

  ierr = input_model->write_input_fields(t_years, dt_years, filename); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAForcing::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
						  string filename) {
  PetscErrorCode ierr;
  double T = t_years + 0.5 * dt_years;

  ierr = input_model->write_diagnostic_fields(t_years, dt_years, filename); CHKERRQ(ierr);

  // also write temperature and precipitation offsets (if used)

  if (temp_ma_anomaly != NULL) {
    ierr = temp_ma_anomaly->update(t_years, dt_years); CHKERRQ(ierr);
    ierr = temp_ma_anomaly->interp(T); CHKERRQ(ierr);
    ierr = temp_ma_anomaly->write(filename.c_str()); CHKERRQ(ierr);
  }

  if (snowprecip_anomaly != NULL) {
    ierr = snowprecip_anomaly->update(t_years, dt_years); CHKERRQ(ierr);
    ierr = snowprecip_anomaly->interp(T); CHKERRQ(ierr);
    ierr = snowprecip_anomaly->write(filename.c_str()); CHKERRQ(ierr);
  }

  if (dTforcing != NULL) {
    delta_T->output_filename = filename;
    ierr = delta_T->append(T, (*dTforcing)(T)); CHKERRQ(ierr);
    ierr = delta_T->interp(T); CHKERRQ(ierr);
    ierr = delta_T->flush(); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode PAForcing::write_fields(set<string> vars, PetscReal t_years,
						   PetscReal dt_years, string filename) {
  PetscErrorCode ierr;
  double T = t_years + 0.5 * dt_years;

  if (temp_ma_anomaly != NULL) {
    if (vars.find("temp_ma_anomaly") != vars.end()) {
      ierr = temp_ma_anomaly->update(t_years, dt_years); CHKERRQ(ierr);
      ierr = temp_ma_anomaly->interp(T); CHKERRQ(ierr);
      ierr = temp_ma_anomaly->write(filename.c_str()); CHKERRQ(ierr);
      vars.erase("temp_ma_anomaly");
    }
  }

  if (snowprecip_anomaly != NULL) {
    if (vars.find("snowprecip_anomaly") != vars.end()) {
      ierr = snowprecip_anomaly->update(t_years, dt_years); CHKERRQ(ierr);
      ierr = snowprecip_anomaly->interp(T); CHKERRQ(ierr);
      ierr = snowprecip_anomaly->write(filename.c_str()); CHKERRQ(ierr);
      vars.erase("snowprecip_anomaly");
    }

  }

  if (dTforcing != NULL) {
    if (vars.find("delta_T") != vars.end()) {
      delta_T->output_filename = filename;
      ierr = delta_T->append(T, (*dTforcing)(T)); CHKERRQ(ierr);
      ierr = delta_T->interp(T); CHKERRQ(ierr);
      ierr = delta_T->flush(); CHKERRQ(ierr);
      vars.erase("delta_T");
    }
  }

  return 0;
}

PetscErrorCode PAForcing::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  if ((gsl_fcmp(t_years,  t,  1e-4) == 0) &&
      (gsl_fcmp(dt_years, dt, 1e-4) == 0))
    return 0;

  t  = t_years;
  dt = dt_years;

  if (temp_ma_anomaly != NULL) {
    ierr = temp_ma_anomaly->update(t_years, dt_years); CHKERRQ(ierr);
  }

  if (snowprecip_anomaly != NULL) {
    ierr = snowprecip_anomaly->update(t_years, dt_years); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::mean_precip(PetscReal t_years, PetscReal dt_years,
						  IceModelVec2 &result) {
  PetscErrorCode ierr;

  ierr = input_model->mean_precip(t_years, dt_years, result); CHKERRQ(ierr);

  if (snowprecip_anomaly != NULL) {
    string history = "added precipitation anomalies\n" + result.string_attr("history");

    double anomaly;
    ierr = snowprecip_anomaly->update(t_years, dt_years); CHKERRQ(ierr);

    ierr = snowprecip_anomaly->begin_access(); CHKERRQ(ierr);
    ierr = result.begin_access(); CHKERRQ(ierr);
    for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
	ierr = snowprecip_anomaly->average(i, j, t_years, dt_years, anomaly); CHKERRQ(ierr);
	result(i,j) += anomaly;
      }
    }
    ierr = result.end_access(); CHKERRQ(ierr);
    ierr = snowprecip_anomaly->end_access(); CHKERRQ(ierr);

    ierr = result.set_attr("history", history); CHKERRQ(ierr);
  }

  if ((dTforcing != NULL) && paleo_precipitation_correction) {
    string history = "added the paleo-precipitation correction\n" + result.string_attr("history");

    PetscReal precipexpfactor = config.get("precip_exponential_factor_for_temperature");
    ierr = result.scale(exp( precipexpfactor * (*dTforcing)(t_years + 0.5 * dt_years) )); CHKERRQ(ierr);

    ierr = result.set_attr("history", history); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					   IceModelVec2 &result) {
  PetscErrorCode ierr;
  double T = t_years + 0.5 * dt_years;

  ierr = input_model->mean_annual_temp(t_years, dt_years, result); CHKERRQ(ierr);

  if (dTforcing != NULL) {
    string history = "added the temperature offset\n" + result.string_attr("history");

    ierr = result.shift( (*dTforcing)(T) ); CHKERRQ(ierr);

    ierr = result.set_attr("history", history); CHKERRQ(ierr);
  }
  
  if (temp_ma_anomaly != NULL) {
    string history = "added temperature anomalies\n" + result.string_attr("history");

    ierr = temp_ma_anomaly->update(t_years, dt_years); CHKERRQ(ierr);
    ierr = temp_ma_anomaly->interp(T); CHKERRQ(ierr);
    ierr = result.add(1.0, *temp_ma_anomaly); CHKERRQ(ierr);

    ierr = result.set_attr("history", history); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::begin_pointwise_access() {
  PetscErrorCode ierr;

  ierr = input_model->begin_pointwise_access(); CHKERRQ(ierr);

  if (temp_ma_anomaly != NULL) {
    ierr = temp_ma_anomaly->begin_access(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::end_pointwise_access() {
  PetscErrorCode ierr;

  ierr = input_model->end_pointwise_access(); CHKERRQ(ierr);

  if (temp_ma_anomaly != NULL) {
    ierr = temp_ma_anomaly->end_access(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAForcing::temp_time_series(int i, int j, int N,
					   PetscReal *ts, PetscReal *values) {
  PetscErrorCode ierr;
  
  ierr = input_model->temp_time_series(i, j, N, ts, values); CHKERRQ(ierr);

  if (dTforcing != NULL) {
    for (PetscInt k = 0; k < N; k++)
      values[k] += (*dTforcing)(ts[k]);
  }

  if (temp_ma_anomaly != NULL) {
    PetscScalar *tmp = new PetscScalar[N];

    ierr = temp_ma_anomaly->interp(i, j, N, ts, tmp); CHKERRQ(ierr);

    for (PetscInt k = 0; k < N; k++)
      values[k] += tmp[k];

    delete tmp;
  }

  return 0;
}

PetscErrorCode PAForcing::temp_snapshot(PetscReal t_years, PetscReal dt_years,
					IceModelVec2 &result) {
  PetscErrorCode ierr;
  double T = t_years + 0.5 * dt_years;

  ierr = input_model->temp_snapshot(t_years, dt_years, result); CHKERRQ(ierr);

  if (temp_ma_anomaly != NULL) {
    string history = "added temperature anomalies\n" + result.string_attr("history");

    ierr = temp_ma_anomaly->update(t_years, dt_years); CHKERRQ(ierr);
    ierr = temp_ma_anomaly->interp(T); CHKERRQ(ierr);
    ierr = result.add(1.0, *temp_ma_anomaly); CHKERRQ(ierr);

    ierr = result.set_attr("history", history); CHKERRQ(ierr);
  }

  return 0;
}

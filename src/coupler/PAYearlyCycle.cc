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

//! Allocates memory and reads in the snow precipitaion data.
PetscErrorCode PAYearlyCycle::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;
  LocalInterpCtx *lic = NULL;
  bool regrid = false;
  int start = -1;

  // Allocate internal IceModelVecs:
  ierr = temp_ma.create(grid, "airtemp_ma", false); CHKERRQ(ierr);
  ierr = temp_ma.set_attrs("diagnostic",
			   "mean annual near-surface air temperature",
			   "K", 
			   ""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = temp_ma.set_attr("source", reference);

  ierr = temp_mj.create(grid, "airtemp_mj", false); CHKERRQ(ierr);
  ierr = temp_mj.set_attrs("diagnostic",
			   "mean July near-surface air temperature",
			   "Kelvin",
			   ""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = temp_mj.set_attr("source", reference);

  ierr = snowprecip.create(grid, "snowprecip", false); CHKERRQ(ierr);
  ierr = snowprecip.set_attrs("climate_state", 
			      "mean annual ice-equivalent snow precipitation rate",
			      "m s-1", 
			      ""); CHKERRQ(ierr); // no CF standard_name ??
  ierr = snowprecip.set_glaciological_units("m year-1");
  snowprecip.write_in_glaciological_units = true;
  snowprecip.time_independent = true;

  ierr = find_pism_input(snowprecip_filename, lic, regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate from file
  ierr = verbPrintf(2, grid.com, 
		    "    reading mean annual ice-equivalent snow precipitation rate 'snowprecip'\n"
		    "      from %s ... \n",
		    snowprecip_filename.c_str()); CHKERRQ(ierr); 
  if (regrid) {
    ierr = snowprecip.regrid(snowprecip_filename.c_str(), *lic, true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = snowprecip.read(snowprecip_filename.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }
  string snowprecip_history = "read from " + snowprecip_filename + "\n";

  ierr = snowprecip.set_attr("history", snowprecip_history); CHKERRQ(ierr);

  delete lic;

  return 0;
}

PetscErrorCode PAYearlyCycle::write_model_state(PetscReal /*t_years*/, PetscReal /*dt_years*/,
								string filename) {
  PetscErrorCode ierr;

  ierr = snowprecip.write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAYearlyCycle::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
								     string filename) {
  PetscErrorCode ierr;

  ierr = write_model_state(t_years, dt_years, filename); CHKERRQ(ierr);

  ierr = update(t_years, dt_years); CHKERRQ(ierr);

  ierr = temp_ma.write(filename.c_str(), NC_FLOAT); CHKERRQ(ierr);
  ierr = temp_mj.write(filename.c_str(), NC_FLOAT); CHKERRQ(ierr);

  // Compute a snapshot of the instantaneous air temperature for this time-step
  // using the standard yearly cycle and write that too:
  IceModelVec2S tmp;		// will be de-allocated at the end of scope

  ierr = tmp.create(grid, "airtemp", false); CHKERRQ(ierr);
  ierr = tmp.set_attrs("diagnostic", "near-surface air temperature snapshot",
		       "K", ""); CHKERRQ(ierr);

  ierr = temp_snapshot(t_years, dt_years, tmp); CHKERRQ(ierr);
  ierr = tmp.write(filename.c_str(), NC_FLOAT); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAYearlyCycle::write_fields(set<string> vars, PetscReal t_years,
					   PetscReal /*dt_years*/, string filename) {
  PetscErrorCode ierr;

  ierr = update(t_years, 0); CHKERRQ(ierr);

  if (vars.find("airtemp") != vars.end()) {
    IceModelVec2S airtemp;
    ierr = airtemp.create(grid, "airtemp", false); CHKERRQ(ierr);
    ierr = airtemp.set_attrs("diagnostic",
			     "snapshot of the near-surface air temperature",
			     "K",
			     ""); CHKERRQ(ierr);

    ierr = temp_snapshot(t_years, 0, airtemp); CHKERRQ(ierr);

    ierr = airtemp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (vars.find("airtemp_ma") != vars.end()) {
    ierr = temp_ma.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (vars.find("airtemp_mj") != vars.end()) {
    ierr = temp_mj.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (vars.find("snowprecip") != vars.end()) {
    ierr = snowprecip.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

//! Copies the stored snow precipitation field into result.
PetscErrorCode PAYearlyCycle::mean_precip(PetscReal t_years, PetscReal dt_years,
							 IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = update(t_years, dt_years); CHKERRQ(ierr);

  string snowprecip_history = "read from " + snowprecip_filename + "\n";

  ierr = snowprecip.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", snowprecip_history); CHKERRQ(ierr);

  return 0;
}

//! Copies the stored mean annual near-surface air temperature field into result.
PetscErrorCode PAYearlyCycle::mean_annual_temp(PetscReal t_years, PetscReal dt_years,
							      IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = update(t_years, dt_years); CHKERRQ(ierr);

  ierr = temp_ma.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history",
			 "computed using " + reference + "\n"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAYearlyCycle::temp_time_series(int i, int j, int N,
							      PetscReal *ts, PetscReal *values) {
  // constants related to the standard yearly cycle
  const PetscReal
    radpersec = 2.0 * pi / secpera, // radians per second frequency for annual cycle
    sperd = 8.64e4, // exact number of seconds per day
    julydaysec = sperd * config.get("snow_temp_july_day");

  for (PetscInt k = 0; k < N; ++k) {
    double tk = ( ts[k] - floor(ts[k]) ) * secpera; // time from the beginning of a year, in seconds
    values[k] = temp_ma(i,j) + (temp_mj(i,j) - temp_ma(i,j)) * cos(radpersec * (tk - julydaysec));
  }

  return 0;
}

PetscErrorCode PAYearlyCycle::temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2S &result) {
  PetscErrorCode ierr;
  const PetscReal
    radpersec = 2.0 * pi / secpera, // radians per second frequency for annual cycle
    sperd = 8.64e4, // exact number of seconds per day
    julydaysec = sperd * config.get("snow_temp_july_day");

  ierr = update(t_years, dt_years); CHKERRQ(ierr);

  double T = t_years + 0.5 * dt_years;

  double t_sec = ( T - floor(T) ) * secpera; // time from the beginning of a year, in seconds

  ierr = temp_mj.add(-1.0, temp_ma, result); CHKERRQ(ierr); // tmp = temp_mj - temp_ma
  ierr = result.scale(cos(radpersec * (t_sec - julydaysec))); CHKERRQ(ierr);
  ierr = result.add(1.0, temp_ma); CHKERRQ(ierr);
  // result = temp_ma + (temp_mj - temp_ma) * cos(radpersec * (T - julydaysec));

  string history = "computed using " + reference + "\n";
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAYearlyCycle::begin_pointwise_access() {
  PetscErrorCode ierr;

  ierr = temp_ma.begin_access(); CHKERRQ(ierr);
  ierr = temp_mj.begin_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAYearlyCycle::end_pointwise_access() {
  PetscErrorCode ierr;

  ierr = temp_ma.end_access(); CHKERRQ(ierr);
  ierr = temp_mj.end_access(); CHKERRQ(ierr);

  return 0;
}

///// PA_SeaRISE_Greenland

PetscErrorCode PA_SeaRISE_Greenland::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
		    "* Initializing SeaRISE-Greenland atmosphere model based on the Fausto et al (2009)\n"
		    "  air temperature parameterization and using stored time-independent precipitation...\n"); CHKERRQ(ierr);
  
  reference =
    "R. S. Fausto, A. P. Ahlstrom, D. V. As, C. E. Boggild, and S. J. Johnsen, 2009. "
    "A new present-day temperature parameterization for Greenland. J. Glaciol. 55 (189), 95-105.";

  ierr = PAYearlyCycle::init(vars); CHKERRQ(ierr);

  // initialize pointers to fields the parameterization depends on:
  surfelev = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (!surfelev) SETERRQ(1, "ERROR: surface_altitude is not available");

  lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
  if (!lat) SETERRQ(1, "ERROR: latitude is not available");

  lon = dynamic_cast<IceModelVec2S*>(vars.get("longitude"));
  if (!lon) SETERRQ(1, "ERROR: longitude is not available");

  ierr = PISMOptionsIsSet("-paleo_precip", paleo_precipitation_correction); CHKERRQ(ierr);

  if (paleo_precipitation_correction) {
    PetscTruth dTforcing_set;
    char dT_file[PETSC_MAX_PATH_LEN];

    ierr = PetscOptionsString("-dTforcing", "Specifies the air temperature offsets file",
			      "", "",
			      dT_file, PETSC_MAX_PATH_LEN, &dTforcing_set); CHKERRQ(ierr);

    if (!dTforcing_set) {
      ierr = PetscPrintf(grid.com, "ERROR: option -paleo_precip requires -dTforcing.\n"); CHKERRQ(ierr);
      PetscEnd();
    }

    dTforcing = new Timeseries(grid.com, grid.rank, "delta_T", "t");
    ierr = dTforcing->set_units("Celsius", ""); CHKERRQ(ierr);
    ierr = dTforcing->set_dimension_units("years", ""); CHKERRQ(ierr);
    ierr = dTforcing->set_attr("long_name", "near-surface air temperature offsets"); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com, 
		      "  reading delta T data from forcing file %s...\n", dT_file);  CHKERRQ(ierr);
	 
    ierr = dTforcing->read(dT_file); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PA_SeaRISE_Greenland::mean_precip(PetscReal t_years, PetscReal dt_years,
						IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = PAYearlyCycle::mean_precip(t_years, dt_years, result); CHKERRQ(ierr);

  if ((dTforcing != NULL) && paleo_precipitation_correction) {
    string history = "added the paleo-precipitation correction\n" + result.string_attr("history");

    PetscReal precipexpfactor = config.get("precip_exponential_factor_for_temperature");
    ierr = result.scale(exp( precipexpfactor * (*dTforcing)(t_years + 0.5 * dt_years) )); CHKERRQ(ierr);

    ierr = result.set_attr("history", history); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Updates mean annual and mean July near-surface air temperatures.
//! Note that the snow precipitation rate is time-independent and does not need
//! to be updated.
PetscErrorCode PA_SeaRISE_Greenland::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  if ((fabs(t_years - t) < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  t  = t_years;
  dt = dt_years;

  const PetscReal 
    d_ma     = config.get("snow_temp_fausto_d_ma"),      // K
    gamma_ma = config.get("snow_temp_fausto_gamma_ma"),  // K m-1
    c_ma     = config.get("snow_temp_fausto_c_ma"),      // K (degN)-1
    kappa_ma = config.get("snow_temp_fausto_kappa_ma"),  // K (degW)-1
    d_mj     = config.get("snow_temp_fausto_d_mj"),      // SAME UNITS as for _ma ...
    gamma_mj = config.get("snow_temp_fausto_gamma_mj"),
    c_mj     = config.get("snow_temp_fausto_c_mj"),
    kappa_mj = config.get("snow_temp_fausto_kappa_mj");
  
  PetscScalar **lat_degN, **lon_degE, **h;
  ierr = surfelev->get_array(h);   CHKERRQ(ierr);
  ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  ierr = lon->get_array(lon_degE); CHKERRQ(ierr);
  ierr = temp_ma.begin_access();  CHKERRQ(ierr);
  ierr = temp_mj.begin_access();  CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      temp_ma(i,j) = d_ma + gamma_ma * h[i][j] + c_ma * lat_degN[i][j] + kappa_ma * (-lon_degE[i][j]);
      temp_mj(i,j) = d_mj + gamma_mj * h[i][j] + c_mj * lat_degN[i][j] + kappa_mj * (-lon_degE[i][j]);
    }
  }
  
  ierr = surfelev->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = lon->end_access(); CHKERRQ(ierr);
  ierr = temp_ma.end_access();  CHKERRQ(ierr);
  ierr = temp_mj.end_access();  CHKERRQ(ierr);

  return 0;
}

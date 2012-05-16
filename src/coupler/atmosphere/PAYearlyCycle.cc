// Copyright (C) 2008-2012 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

// Implementation of the atmosphere model using constant-in-time precipitation
// and a cosine yearly cycle for near-surface air temperatures.

#include "PAYearlyCycle.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"

//! Allocates memory and reads in the precipitaion data.
PetscErrorCode PAYearlyCycle::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool regrid = false;
  int start = -1;

  variables = &vars;

  snow_temp_july_day = config.get("snow_temp_july_day");

  // Allocate internal IceModelVecs:
  ierr = air_temp_mean_annual.create(grid, "air_temp_mean_annual", false); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.set_attrs("diagnostic",
			   "mean annual near-surface air temperature (without sub-year time-dependence or forcing)",
			   "K", 
			   ""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = air_temp_mean_annual.set_attr("source", reference);

  ierr = air_temp_mean_july.create(grid, "air_temp_mean_july", false); CHKERRQ(ierr);
  ierr = air_temp_mean_july.set_attrs("diagnostic",
			   "mean July near-surface air temperature (without sub-year time-dependence or forcing)",
			   "Kelvin",
			   ""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = air_temp_mean_july.set_attr("source", reference);

  ierr = precipitation.create(grid, "precipitation", false); CHKERRQ(ierr);
  ierr = precipitation.set_attrs("climate_state", 
			      "mean annual ice-equivalent precipitation rate",
			      "m s-1", 
			      ""); CHKERRQ(ierr); // no CF standard_name ??
  ierr = precipitation.set_glaciological_units("m year-1");
  precipitation.write_in_glaciological_units = true;
  precipitation.time_independent = true;

  ierr = find_pism_input(precip_filename, regrid, start); CHKERRQ(ierr);

  // read precipitation rate from file
  ierr = verbPrintf(2, grid.com, 
		    "    reading mean annual ice-equivalent precipitation rate 'precipitation'\n"
		    "      from %s ... \n",
		    precip_filename.c_str()); CHKERRQ(ierr); 
  if (regrid) {
    ierr = precipitation.regrid(precip_filename.c_str(), true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = precipitation.read(precip_filename.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }
  string precip_history = "read from " + precip_filename + "\n";

  ierr = precipitation.set_attr("history", precip_history); CHKERRQ(ierr);

  air_temp_snapshot.init_2d("air_temp_snapshot", grid);
  air_temp_snapshot.set_string("pism_intent", "diagnostic");
  air_temp_snapshot.set_string("long_name",
                         "snapshot of the near-surface air temperature");
  ierr = air_temp_snapshot.set_units("K"); CHKERRQ(ierr);

  return 0;
}

void PAYearlyCycle::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  result["precipitation"] = precipitation.get_metadata();

  if (keyword == "big") {
    result["air_temp_mean_annual"] = air_temp_mean_annual.get_metadata();
    result["air_temp_mean_july"] = air_temp_mean_july.get_metadata();
    result["air_temp_snapshot"] = air_temp_snapshot;
  }
}


PetscErrorCode PAYearlyCycle::define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp_snapshot")) {
    ierr = air_temp_snapshot.define(nc, nctype, false); CHKERRQ(ierr);
  }

  if (set_contains(vars, "air_temp_mean_annual")) {
    ierr = air_temp_mean_annual.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "air_temp_mean_july")) {
    ierr = air_temp_mean_july.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "precipitation")) {
    ierr = precipitation.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode PAYearlyCycle::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp_snapshot")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "air_temp_snapshot", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(air_temp_snapshot, 0); CHKERRQ(ierr);

    ierr = temp_snapshot(tmp); CHKERRQ(ierr);

    ierr = tmp.write(filename); CHKERRQ(ierr);
  }

  if (set_contains(vars, "air_temp_mean_annual")) {
    ierr = air_temp_mean_annual.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "air_temp_mean_july")) {
    ierr = air_temp_mean_july.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "precipitation")) {
    ierr = precipitation.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

//! Copies the stored precipitation field into result.
PetscErrorCode PAYearlyCycle::mean_precipitation(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = precipitation.copy_to(result); CHKERRQ(ierr);
  return 0;
}

//! Copies the stored mean annual near-surface air temperature field into result.
PetscErrorCode PAYearlyCycle::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = air_temp_mean_annual.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAYearlyCycle::temp_time_series(int i, int j, int N,
                                               PetscReal *ts, PetscReal *values) {
  // constants related to the standard yearly cycle
  const PetscReal
    sperd = 8.64e4, // exact number of seconds per day
    julyday_fraction = (sperd / secpera) * snow_temp_july_day;

  for (PetscInt k = 0; k < N; ++k) {
    double tk = grid.time->year_fraction(ts[k]) - julyday_fraction;
    values[k] = air_temp_mean_annual(i,j) +
      (air_temp_mean_july(i,j) - air_temp_mean_annual(i,j)) * cos(2.0 * pi * tk);
  }

  return 0;
}

PetscErrorCode PAYearlyCycle::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr;
  const PetscReal
    sperd = 8.64e4, // exact number of seconds per day
    julyday_fraction = (sperd / secpera) * snow_temp_july_day;

  double T = grid.time->year_fraction(t + 0.5 * dt) - julyday_fraction;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.begin_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_july.begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i,j) = air_temp_mean_annual(i,j) +
        (air_temp_mean_july(i,j) - air_temp_mean_annual(i,j)) * cos(2.0 * pi * T);
    }
  }

  ierr = air_temp_mean_july.end_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAYearlyCycle::begin_pointwise_access() {
  PetscErrorCode ierr;

  ierr = air_temp_mean_annual.begin_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_july.begin_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAYearlyCycle::end_pointwise_access() {
  PetscErrorCode ierr;

  ierr = air_temp_mean_annual.end_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_july.end_access(); CHKERRQ(ierr);

  return 0;
}


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

PetscErrorCode PAConstant::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;
  LocalInterpCtx *lic = NULL;
  bool regrid = false;
  int start = -1;

  ierr = verbPrintf(2, grid.com, "* Initializing the constant-in-time atmosphere model...\n"); CHKERRQ(ierr);

  // allocate IceModelVecs for storing temperature and precipitation fields:

  // create mean annual ice equivalent snow precipitation rate (before melt, and not including rain)
  ierr = snowprecip.create(grid, "snowprecip", false); CHKERRQ(ierr);
  ierr = snowprecip.set_attrs("climate_state", 
			      "mean annual ice-equivalent snow precipitation rate",
			      "m s-1", 
			      ""); CHKERRQ(ierr); // no CF standard_name ??
  ierr = snowprecip.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  snowprecip.write_in_glaciological_units = true;
  snowprecip.time_independent = true;

  ierr = temperature.create(grid, "temp_ma", false); CHKERRQ(ierr); // FIXME! choose the right name
  ierr = temperature.set_attrs(
            "climate_state",
            "mean annual near-surface (2 m) air temperature",
            "K",
	    ""); CHKERRQ(ierr);
  temperature.time_independent = true;
  
  // find PISM input file to read data from:

  ierr = find_pism_input(input_file, lic, regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate and temperatures from file
  ierr = verbPrintf(2, grid.com, 
		    "    reading mean annual ice-equivalent snow precipitation rate 'snowprecip'\n"
		    "    and mean annual near-surface air temperature 'temp_ma' from %s ... \n",
		    input_file.c_str()); CHKERRQ(ierr); 
  if (regrid) {
    ierr = snowprecip.regrid(input_file.c_str(), *lic, true); CHKERRQ(ierr); // fails if not found!
    ierr = temperature.regrid(input_file.c_str(), *lic, true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = snowprecip.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
    ierr = temperature.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }

  delete lic;

  t = grid.year;
  dt = 0;
	    
  return 0;
}

PetscErrorCode PAConstant::mean_precip(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						   IceModelVec2 &result) {
  PetscErrorCode ierr;

  string snowprecip_history = "read from " + input_file + "\n";

  ierr = snowprecip.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", snowprecip_history);
  return 0;
}

PetscErrorCode PAConstant::mean_annual_temp(PetscReal /*t_years*/, PetscReal /*dt_years*/,
							IceModelVec2 &result) {
  PetscErrorCode ierr;

  string temp_history = "read from " + input_file + "\n";

  ierr = temperature.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", temp_history);
  return 0;
}

PetscErrorCode PAConstant::begin_pointwise_access() {
  PetscErrorCode ierr;
  ierr = temperature.begin_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAConstant::end_pointwise_access() {
  PetscErrorCode ierr;
  ierr = temperature.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAConstant::temp_time_series(int i, int j, int N,
					    PetscReal */*ts*/, PetscReal *values) {
  for (PetscInt k = 0; k < N; k++)
    values[k] = temperature(i,j);
  return 0;
}

PetscErrorCode PAConstant::temp_snapshot(PetscReal t_years, PetscReal dt_years,
					 IceModelVec2 &result) {
  PetscErrorCode ierr;

  ierr = mean_annual_temp(t_years, dt_years, result); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAConstant::write_input_fields(PetscReal /*t_years*/, PetscReal /*dt_years*/,
					      string filename) {
  PetscErrorCode ierr;

  ierr = snowprecip.write(filename.c_str()); CHKERRQ(ierr);
  ierr = temperature.write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAConstant::write_diagnostic_fields(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						   string filename) {
  PetscErrorCode ierr;

  ierr = snowprecip.write(filename.c_str()); CHKERRQ(ierr);
  ierr = temperature.write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAConstant::write_fields(set<string> vars, PetscReal t_years,
					PetscReal dt_years, string filename) {
  PetscErrorCode ierr;

  if (vars.find("airtemp") != vars.end()) {
    IceModelVec2 airtemp;
    ierr = airtemp.create(grid, "airtemp", false); CHKERRQ(ierr);
    ierr = airtemp.set_attrs("diagnostic",
			     "snapshot of the near-surface air temperature",
			     "K",
			     ""); CHKERRQ(ierr);

    ierr = temp_snapshot(t_years, dt_years, airtemp); CHKERRQ(ierr);

    ierr = airtemp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (vars.find("snowprecip") != vars.end()) {
    ierr = snowprecip.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

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

// Implementation of the constant-in-time atmosphere model (-atmosphere
// constant).

#include "PISMAtmosphere.hh"
#include "IceGrid.hh"
#include "PAConstant.hh"

PetscErrorCode PAConstant::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;
  bool regrid = false;
  int start = -1;

  ierr = verbPrintf(2, grid.com, "* Initializing the constant-in-time atmosphere model...\n"); CHKERRQ(ierr);

  // allocate IceModelVecs for storing temperature and precipitation fields:

  // create mean annual ice equivalent precipitation rate (before separating rain,
  //   and before melt, etc. in PISMSurfaceModel)
  ierr = precip.create(grid, "precip", false); CHKERRQ(ierr);
  ierr = precip.set_attrs("climate_state", 
			      "mean annual ice-equivalent precipitation rate",
			      "m s-1",
			      ""); CHKERRQ(ierr); // no CF standard_name ??
  ierr = precip.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  precip.write_in_glaciological_units = true;
  precip.time_independent = true;

  ierr = temperature.create(grid, "temp_ma", false); CHKERRQ(ierr); // FIXME! choose the right name
  ierr = temperature.set_attrs(
            "climate_state",
            "mean annual near-surface (2 m) air temperature",
            "K",
	    ""); CHKERRQ(ierr);
  temperature.time_independent = true;
  
  // find PISM input file to read data from:

  ierr = find_pism_input(input_file, regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate and temperatures from file
  ierr = verbPrintf(2, grid.com, 
		    "    reading mean annual ice-equivalent precipitation rate 'precip'\n"
		    "    and mean annual near-surface air temperature 'temp_ma' from %s ... \n",
		    input_file.c_str()); CHKERRQ(ierr); 
  if (regrid) {
    ierr = precip.regrid(input_file.c_str(), true); CHKERRQ(ierr); // fails if not found!
    ierr = temperature.regrid(input_file.c_str(), true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = precip.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
    ierr = temperature.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }

  airtemp_var.init_2d("airtemp", grid);
  airtemp_var.set_string("pism_intent", "diagnostic");
  airtemp_var.set_string("long_name",
                         "snapshot of the near-surface air temperature");
  ierr = airtemp_var.set_units("K"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAConstant::mean_precip(IceModelVec2S &result) {
  PetscErrorCode ierr;

  string precip_history = "read from " + input_file + "\n";

  ierr = precip.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", precip_history);
  return 0;
}

PetscErrorCode PAConstant::mean_annual_temp(IceModelVec2S &result) {
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

PetscErrorCode PAConstant::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = mean_annual_temp(result); CHKERRQ(ierr);

  return 0;
}

void PAConstant::add_vars_to_output(string keyword, set<string> &result) {
  result.insert("precip");
  result.insert("temp_ma");
  
  if (keyword == "big") {
    result.insert("airtemp");
  }
}

PetscErrorCode PAConstant::define_variables(set<string> vars, const NCTool &nc,
                                            nc_type nctype) {
  PetscErrorCode ierr;
  int varid;

  if (set_contains(vars, "airtemp")) {
    ierr = airtemp_var.define(nc, varid, nctype, false); CHKERRQ(ierr);
  }

  if (set_contains(vars, "precip")) {
    ierr = precip.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "temp_ma")) {
    ierr = temperature.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAConstant::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "airtemp")) {
    IceModelVec2S airtemp;
    ierr = airtemp.create(grid, "airtemp", false); CHKERRQ(ierr);
    ierr = airtemp.set_metadata(airtemp_var, 0); CHKERRQ(ierr);

    ierr = temp_snapshot(airtemp); CHKERRQ(ierr);

    ierr = airtemp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "precip")) {
    ierr = precip.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "temp_ma")) {
    ierr = temperature.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}


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

#include "PISMAtmosphere.hh"

PetscErrorCode PAConstantPIK::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool regrid = false;
  int start = -1;

  ierr = verbPrintf(2, grid.com,
     "* Initializing the constant-in-time atmosphere model PAConstantPIK.\n"
     "  It reads a precipitation field directly from the file and holds it constant.\n"
     "  Near-surface air temperature is parameterized as in Martin et al. 2011, Eqn. 2.0.2.\n"); CHKERRQ(ierr);

  // allocate IceModelVecs for storing temperature and precipitation fields:

  // create mean annual ice equivalent precipitation rate (before separating
  // rain, and before melt, etc. in PISMSurfaceModel)
  ierr = precip.create(grid, "precip", false); CHKERRQ(ierr);
  ierr = precip.set_attrs("climate_state", 
                          "mean annual ice-equivalent precipitation rate",
                          "m s-1",
                          ""); CHKERRQ(ierr); // no CF standard_name ??
  ierr = precip.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  precip.write_in_glaciological_units = true;
  precip.time_independent = true;

  ierr = temperature.create(grid, "temp_ma", false); CHKERRQ(ierr); // FIXME! choose the right name
  ierr = temperature.set_attrs("climate_state",
                               "mean annual near-surface (2 m) air temperature",
                               "K",
                               ""); CHKERRQ(ierr);
  temperature.time_independent = true;
  
  // find PISM input file to read data from:

  ierr = find_pism_input(input_file, regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate and temperatures from file
  ierr = verbPrintf(2, grid.com, 
		    "    reading mean annual ice-equivalent precipitation rate 'precip'\n"
		    "    from %s ... \n",
		    input_file.c_str()); CHKERRQ(ierr); 
  if (regrid) {
    ierr = precip.regrid(input_file.c_str(), true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = precip.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }

  usurf = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (usurf == NULL) SETERRQ(1, "surface_altitude is not available");

  lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
  if (lat == NULL) SETERRQ(1, "latitude is not available");

  airtemp_var.init_2d("airtemp", grid);
  airtemp_var.set_string("pism_intent", "diagnostic");
  airtemp_var.set_string("long_name",
                         "snapshot of the near-surface air temperature");
  ierr = airtemp_var.set_units("K"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAConstantPIK::update(PetscReal, PetscReal) {
  PetscErrorCode ierr;

  // Compute near-surface air temperature using a latitude- and
  // elevation-dependent parameterization:

  ierr = temperature.begin_access();   CHKERRQ(ierr);
  ierr = usurf->begin_access();   CHKERRQ(ierr);
  ierr = lat->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      temperature(i, j) = 273.15 + 30 - 0.0075 * (*usurf)(i,j) - 0.68775 * (*lat)(i,j)*(-1.0) ;

    }
  }
  ierr = usurf->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = temperature.end_access();   CHKERRQ(ierr);

  return 0;
}


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

#include "PSConstant.hh"
#include "IceGrid.hh"

///// Constant-in-time surface model.

PetscErrorCode PSConstant::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;
  bool regrid = false;
  int start = -1;

  ierr = verbPrintf(2, grid.com,
     "* Initializing the constant-in-time surface processes model PSConstant.\n"
     "  It reads surface mass balance and ice upper-surface temperature\n"
     "  directly from the file and holds them constant.\n"
     "  Any choice of atmosphere coupler (option '-atmosphere') is ignored.\n"); CHKERRQ(ierr);

  // allocate IceModelVecs for storing temperature and surface mass balance fields

  // create mean annual ice equivalent snow precipitation rate (before melt, and not including rain)
  ierr = acab.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = acab.set_attrs("climate_state",
			"constant-in-time ice-equivalent surface mass balance (accumulation/ablation) rate",
			"m s-1",
			"land_ice_surface_specific_mass_balance"); // CF standard_name
			CHKERRQ(ierr);
  ierr = acab.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  acab.write_in_glaciological_units = true;

  ierr = artm.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = artm.set_attrs("climate_state",
			"constant-in-time ice temperature at the ice surface",
			"K",
			""); CHKERRQ(ierr);

  // find PISM input file to read data from:
  ierr = find_pism_input(input_file, regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate and temperatures from file
  ierr = verbPrintf(2, grid.com,
    "    reading ice-equivalent surface mass balance rate 'acab'\n"
    "    and ice surface temperature  'artm' from %s ... \n",
    input_file.c_str()); CHKERRQ(ierr);
  if (regrid) {
    ierr = acab.regrid(input_file.c_str(), true); CHKERRQ(ierr); // fails if not found!
    ierr = artm.regrid(input_file.c_str(), true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = acab.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
    ierr = artm.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }


  return 0;
}

PetscErrorCode PSConstant::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  string history  = "read from " + input_file + "\n";

  ierr = acab.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSConstant::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;
  string history  = "read from " + input_file + "\n";

  ierr = artm.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

void PSConstant::add_vars_to_output(string /*keyword*/, set<string> &result) {
  result.insert("acab");
  result.insert("artm");
  // does not call atmosphere->add_vars_to_output().
}

PetscErrorCode PSConstant::define_variables(set<string> vars, const NCTool &nc, nc_type nctype) {
  PetscErrorCode ierr;

  ierr = PISMSurfaceModel::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  if (set_contains(vars, "artm")) {
    ierr = artm.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "acab")) {
    ierr = acab.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PSConstant::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "artm")) {
    ierr = artm.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "acab")) {
    ierr = acab.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

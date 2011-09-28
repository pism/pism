// Copyright (C) 2008-2011 PISM Authors
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

#include "PISMSurface.hh"
#include "PISMIO.hh"


///// Constant-in-time surface model for accumulation,
///// ice surface temperature parameterized as in PISM-PIK dependent on latitude and surface elevation

PetscErrorCode PSConstantPIK::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool regrid = false;
  int start = -1;

  //variables = &vars;

  ierr = verbPrintf(2, grid.com,
     "* Initializing the constant-in-time surface processes model PSConstantPIK.\n"
     "  It reads surface mass balance directly from the file and holds them constant.\n"
     "  Ice upper-surface temperature is parameterized as in Martin et al. 2011, Eqn. 2.0.2.\n"
     "  Any choice of atmosphere coupler (option '-atmosphere') is ignored.\n"); CHKERRQ(ierr);


  usurf = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
   if (!usurf) SETERRQ(12, "ERROR: 'usurf' is not available or is wrong type in dictionary");

  lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
  if (!lat) SETERRQ(1, "ERROR: latitude is not available");



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

  // read snow precipitation rate from file
  ierr = verbPrintf(2, grid.com,
    "    reading ice-equivalent surface mass balance rate 'acab' from %s ... \n",
    input_file.c_str()); CHKERRQ(ierr);
  if (regrid) {
    ierr = acab.regrid(input_file.c_str(), true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = acab.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }


  // parameterizing the ice surface temperature 'artm'
  ierr = verbPrintf(2, grid.com,"    parameterizing the ice surface temperature 'artm' ... \n"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSConstantPIK::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  string history  = "read from " + input_file + "\n";

  ierr = acab.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSConstantPIK::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;
  string history  = "parmeterized ice surface temperature \n";

  ierr = result.begin_access();   CHKERRQ(ierr);
  ierr = artm.begin_access();   CHKERRQ(ierr);
  ierr = usurf->begin_access();   CHKERRQ(ierr);
  ierr = lat->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      result(i,j) = 273.15 + 30 - 0.0075 * (*usurf)(i,j) - 0.68775 * (*lat)(i,j)*(-1.0) ;
      artm(i,j)=result(i,j);

    }
  }
  ierr = usurf->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = result.end_access();   CHKERRQ(ierr);
  ierr = artm.end_access();   CHKERRQ(ierr);

  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

void PSConstantPIK::add_vars_to_output(string /*keyword*/, set<string> &result) {
  result.insert("acab");
  result.insert("artm");
  // does not call atmosphere->add_vars_to_output().
}

PetscErrorCode PSConstantPIK::define_variables(set<string> vars, const NCTool &nc, nc_type nctype) {
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

PetscErrorCode PSConstantPIK::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "artm")) {
    ierr = artm.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "acab")) {
    ierr = acab.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

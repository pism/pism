// Copyright (C) 2011, 2012 PISM Authors
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

#include "PSSimple.hh"
#include "IceGrid.hh"
#include "pism_const.hh"
#include "iceModelVec.hh"

///// Simple PISM surface model.

PetscErrorCode PSSimple::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (atmosphere == NULL)
    SETERRQ(grid.com, 1, "PISMSurfaceModel::init(PISMVars &vars): atmosphere == NULL");

  ierr = atmosphere->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "* Initializing the simplest PISM surface (snow) processes model PSSimple.\n"
     "  It passes atmospheric state directly to upper ice fluid surface:\n"
     "    surface mass balance          := precipitation,\n"
     "    ice upper surface temperature := 2m air temperature.\n");
     CHKERRQ(ierr);

  climatic_mass_balance.init_2d("climatic_mass_balance", grid);
  climatic_mass_balance.set_string("pism_intent", "diagnostic");
  climatic_mass_balance.set_string("long_name",
                  "ice-equivalent surface mass balance (accumulation/ablation) rate");
  climatic_mass_balance.set_string("standard_name",
                  "land_ice_surface_specific_mass_balance");
  ierr = climatic_mass_balance.set_units("m s-1"); CHKERRQ(ierr);
  ierr = climatic_mass_balance.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  ice_surface_temp.init_2d("ice_surface_temp", grid);
  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  ierr = ice_surface_temp.set_units("K"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSSimple::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_precipitation(result); CHKERRQ(ierr);

  string history = result.string_attr("history");
  history = "re-interpreted precipitation as surface mass balance (PSSimple)\n" + history;
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSSimple::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_annual_temp(result); CHKERRQ(ierr);

  string history = result.string_attr("history");
  history = "re-interpreted mean annual 2 m air temperature as instantaneous ice temperature at the ice surface (PSSimple)\n" + history;
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

void PSSimple::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  PISMSurfaceModel::add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result["ice_surface_temp"] = ice_surface_temp;
    result["climatic_mass_balance"] = climatic_mass_balance;
  }
}

PetscErrorCode PSSimple::define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    ierr = ice_surface_temp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = climatic_mass_balance.define(nc, nctype, true); CHKERRQ(ierr);
  }

  ierr = PISMSurfaceModel::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSSimple::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "ice_surface_temp", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(ice_surface_temp, 0); CHKERRQ(ierr);

    ierr = ice_surface_temperature(tmp); CHKERRQ(ierr);

    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);

    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "climatic_mass_balance", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(climatic_mass_balance, 0); CHKERRQ(ierr);

    ierr = ice_surface_mass_flux(tmp); CHKERRQ(ierr);
    tmp.write_in_glaciological_units = true;
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);

    vars.erase("climatic_mass_balance");
  }

  ierr = PISMSurfaceModel::write_variables(vars, filename); CHKERRQ(ierr);

  return 0;
}

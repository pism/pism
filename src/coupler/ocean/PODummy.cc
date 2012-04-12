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

#include "PODummy.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "pism_options.hh"

PODummy::PODummy(IceGrid &g, const NCConfigVariable &conf)
  : PISMOceanModel(g, conf) {

  shelfbmassflux.init_2d("shelfbmassflux", g);
  shelfbmassflux.set_string("pism_intent", "climate_state");
  shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  shelfbmassflux.set_units("m s-1");
  shelfbmassflux.set_glaciological_units("m year-1");

  shelfbtemp.init_2d("shelfbtemp", g);
  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");
}

PetscErrorCode PODummy::init(PISMVars &) {
  return 0;
}

PetscErrorCode PODummy::sea_level_elevation(PetscReal &result) {
  result = sea_level;
  return 0;
}

PetscErrorCode PODummy::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = result.set(273.15); CHKERRQ(ierr);

  return 0;
}

//! Computes mass flux in ice-equivalent m s-1, from assumption that basal heat flux rate converts to mass flux.
PetscErrorCode PODummy::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = result.set(0); CHKERRQ(ierr);

  return 0;
}

void PODummy::add_vars_to_output(string, set<string> &result) {
  result.insert("shelfbtemp");
  result.insert("shelfbmassflux");
}

PetscErrorCode PODummy::define_variables(set<string> vars, const PIO &nc,
                                            PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "shelfbtemp")) {
    ierr = shelfbtemp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    ierr = shelfbmassflux.define(nc, nctype, true); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PODummy::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbtemp, 0); CHKERRQ(ierr);
    ierr = shelf_base_temperature(tmp); CHKERRQ(ierr);
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbmassflux, 0); CHKERRQ(ierr);
    tmp.write_in_glaciological_units = true;
    ierr = shelf_base_mass_flux(tmp); CHKERRQ(ierr);
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

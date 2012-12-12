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

#include "PO_delta_SL.hh"

/// -ocean_delta_SL_file, ...

PO_delta_SL::PO_delta_SL(IceGrid &g, const NCConfigVariable &conf, PISMOceanModel* in)
  : PScalarForcing<PISMOceanModel,POModifier>(g, conf, in)
{
  option_prefix = "-ocean_delta_SL";
  offset_name = "delta_SL";
  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));

  offset->set_units("m", "");
  offset->set_dimension_units(grid.time->units(), "");
  offset->set_attr("long_name", "sea level elevation offsets");

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

PetscErrorCode PO_delta_SL::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "* Initializing sea level forcing...\n"); CHKERRQ(ierr);

  ierr = init_internal(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PO_delta_SL::sea_level_elevation(PetscReal &result) {
  PetscErrorCode ierr = input_model->sea_level_elevation(result); CHKERRQ(ierr);

  if (offset)
    result += (*offset)(t + 0.5*dt);

  return 0;
}

void PO_delta_SL::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  input_model->add_vars_to_output(keyword, result);

  result["shelfbtemp"] = shelfbtemp;
  result["shelfbmassflux"] = shelfbmassflux;
}

PetscErrorCode PO_delta_SL::define_variables(set<string> vars, const PIO &nc,
                                             PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "shelfbtemp")) {
    ierr = shelfbtemp.define(nc, nctype, true); CHKERRQ(ierr);
    vars.erase("shelfbtemp");
  }

  if (set_contains(vars, "shelfbmassflux")) {
    ierr = shelfbmassflux.define(nc, nctype, true); CHKERRQ(ierr);
    vars.erase("shelfbmassflux");
  }

  ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PO_delta_SL::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbtemp, 0); CHKERRQ(ierr);
    ierr = shelf_base_temperature(tmp); CHKERRQ(ierr);
    ierr = tmp.write(nc); CHKERRQ(ierr);
    vars.erase("shelfbtemp");
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbmassflux, 0); CHKERRQ(ierr);
    tmp.write_in_glaciological_units = true;
    ierr = shelf_base_mass_flux(tmp); CHKERRQ(ierr);
    ierr = tmp.write(nc); CHKERRQ(ierr);
    vars.erase("shelfbmassflux");
  }

  ierr = input_model->write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}

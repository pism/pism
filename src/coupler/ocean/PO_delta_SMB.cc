// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "PO_delta_SMB.hh"
#include "PISMConfig.hh"

namespace pism {

PO_delta_SMB::PO_delta_SMB(IceGrid &g, const Config &conf, OceanModel* in)
  : PScalarForcing<OceanModel,POModifier>(g, conf, in),
    shelfbmassflux(g.get_unit_system(), "shelfbmassflux", grid),
    shelfbtemp(g.get_unit_system(), "shelfbtemp", grid)
{
  PetscErrorCode ierr = allocate_PO_delta_SMB(); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    throw std::runtime_error("PO_delta_SMB allocation failed");
  }
}

PO_delta_SMB::~PO_delta_SMB() {
  // empty
}

PetscErrorCode PO_delta_SMB::allocate_PO_delta_SMB() {
  option_prefix = "-ocean_delta_mass_flux";
  offset_name   = "delta_mass_flux";

  offset->get_metadata().set_units("m s-1");
  offset->get_dimension_metadata().set_units(grid.time->units_string());
  offset->get_metadata().set_string("long_name",
                                    "ice-shelf-base mass flux offsets, ice equivalent thickness per time");

  shelfbmassflux.set_string("pism_intent", "climate_state");
  shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  shelfbmassflux.set_units("kg m-2 s-1");
  shelfbmassflux.set_glaciological_units("kg m-2 year-1");

  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");

  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));

  return 0;
}

PetscErrorCode PO_delta_SMB::init(Vars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing ice shelf base mass flux forcing using scalar offsets...\n"); CHKERRQ(ierr);

  ierr = init_internal(); CHKERRQ(ierr);

  // convert from [m s-1] to [kg m-2 s-1]:
  offset->scale(config.get("ice_density"));

  return 0;
}

PetscErrorCode PO_delta_SMB::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->shelf_base_mass_flux(result); CHKERRQ(ierr);
  ierr = offset_data(result); CHKERRQ(ierr);
  return 0;
}

void PO_delta_SMB::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  result.insert("shelfbtemp");
  result.insert("shelfbmassflux");
}

PetscErrorCode PO_delta_SMB::define_variables(const std::set<std::string> &vars_input, const PIO &nc,
                                              IO_Type nctype) {
  std::set<std::string> vars = vars_input;
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

PetscErrorCode PO_delta_SMB::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", WITHOUT_GHOSTS); CHKERRQ(ierr);
    }

    tmp.metadata() = shelfbtemp;
    ierr = shelf_base_temperature(tmp); CHKERRQ(ierr);
    ierr = tmp.write(nc); CHKERRQ(ierr);
    vars.erase("shelfbtemp");
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", WITHOUT_GHOSTS); CHKERRQ(ierr);
    }

    tmp.metadata() = shelfbmassflux;
    tmp.write_in_glaciological_units = true;
    ierr = shelf_base_mass_flux(tmp); CHKERRQ(ierr);
    ierr = tmp.write(nc); CHKERRQ(ierr);
    vars.erase("shelfbmassflux");
  }

  ierr = input_model->write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism

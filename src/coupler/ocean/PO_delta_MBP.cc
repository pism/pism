/* Copyright (C) 2013, 2014 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "PO_delta_MBP.hh"
#include "PISMConfig.hh"

namespace pism {

PO_delta_MBP::PO_delta_MBP(IceGrid &g, const Config &conf, OceanModel* in)
  : PScalarForcing<OceanModel,POModifier>(g, conf, in),
    shelfbmassflux(g.get_unit_system()),
    shelfbtemp(g.get_unit_system())
{
  PetscErrorCode ierr = allocate_PO_delta_MBP(); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    throw std::runtime_error("PO_delta_MBP allocation failed");
  }
}

PO_delta_MBP::~PO_delta_MBP() {
  // empty
}

PetscErrorCode PO_delta_MBP::allocate_PO_delta_MBP() {
  option_prefix = "-ocean_delta_MBP";
  offset_name   = "delta_MBP";

  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));

  offset->get_metadata().set_units("1");
  offset->get_metadata().set_string("long_name", "melange back pressure fraction");
  offset->get_dimension_metadata().set_units(grid.time->units_string());

  shelfbmassflux.init_2d("shelfbmassflux", grid);
  shelfbmassflux.set_string("pism_intent", "climate_state");
  shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  shelfbmassflux.set_units("kg m-2 s-1");
  shelfbmassflux.set_glaciological_units("kg m-2 year-1");

  shelfbtemp.init_2d("shelfbtemp", grid);
  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");

  return 0;
}

PetscErrorCode PO_delta_MBP::init(Vars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "* Initializing melange back pressure fraction forcing...\n"); CHKERRQ(ierr);

  ierr = init_internal(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PO_delta_MBP::melange_back_pressure_fraction(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->melange_back_pressure_fraction(result); CHKERRQ(ierr);

  ierr = result.shift((*offset)(m_t + 0.5*m_dt));

  return 0;
}

void PO_delta_MBP::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  result.insert("shelfbtemp");
  result.insert("shelfbmassflux");
}

PetscErrorCode PO_delta_MBP::define_variables(const std::set<std::string> &vars_input, const PIO &nc,
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

PetscErrorCode PO_delta_MBP::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
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

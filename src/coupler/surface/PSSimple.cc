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

#include "PSSimple.hh"
#include "IceGrid.hh"
#include "pism_const.hh"
#include "iceModelVec.hh"
#include "PISMConfig.hh"

#include <assert.h>

namespace pism {

///// Simple PISM surface model.
PSSimple::PSSimple(IceGrid &g, const Config &conf)
  : SurfaceModel(g, conf),
    climatic_mass_balance(g.get_unit_system()),
    ice_surface_temp(g.get_unit_system())
{
  PetscErrorCode ierr = allocate_PSSimple(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

PetscErrorCode PSSimple::allocate_PSSimple() {
  PetscErrorCode ierr;

  climatic_mass_balance.init_2d("climatic_mass_balance", grid);
  climatic_mass_balance.set_string("pism_intent", "diagnostic");
  climatic_mass_balance.set_string("long_name",
                                   "surface mass balance (accumulation/ablation) rate");
  climatic_mass_balance.set_string("standard_name",
                                   "land_ice_surface_specific_mass_balance");
  ierr = climatic_mass_balance.set_units("kg m-2 s-1"); CHKERRQ(ierr);
  ierr = climatic_mass_balance.set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);

  ice_surface_temp.init_2d("ice_surface_temp", grid);
  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  ierr = ice_surface_temp.set_units("K"); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PSSimple::init(Vars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  assert(atmosphere != NULL);
  ierr = atmosphere->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "* Initializing the simplest PISM surface (snow) processes model PSSimple.\n"
     "  It passes atmospheric state directly to upper ice fluid surface:\n"
     "    surface mass balance          := precipitation,\n"
     "    ice upper surface temperature := 2m air temperature.\n");
     CHKERRQ(ierr);
 
  return 0;
}

PetscErrorCode PSSimple::update(double my_t, double my_dt)
{
  m_t = my_t;
  m_dt = my_dt;
  if (atmosphere) {
    PetscErrorCode ierr = atmosphere->update(my_t, my_dt); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PSSimple::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = atmosphere->mean_precipitation(result); CHKERRQ(ierr);
  ierr = result.scale(config.get("ice_density")); // convert from m/s ice equivalent to kg m-2 s-1
  return 0;
}

PetscErrorCode PSSimple::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = atmosphere->mean_annual_temp(result); CHKERRQ(ierr);
  return 0;
}

void PSSimple::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  SurfaceModel::add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("ice_surface_temp");
    result.insert("climatic_mass_balance");
  }
}

PetscErrorCode PSSimple::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    ierr = ice_surface_temp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = climatic_mass_balance.define(nc, nctype, true); CHKERRQ(ierr);
  }

  ierr = SurfaceModel::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSSimple::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "ice_surface_temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
    tmp.metadata() = ice_surface_temp;

    ierr = ice_surface_temperature(tmp); CHKERRQ(ierr);

    ierr = tmp.write(nc); CHKERRQ(ierr);

    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "climatic_mass_balance", WITHOUT_GHOSTS); CHKERRQ(ierr);
    tmp.metadata() = climatic_mass_balance;

    ierr = ice_surface_mass_flux(tmp); CHKERRQ(ierr);
    tmp.write_in_glaciological_units = true;
    ierr = tmp.write(nc); CHKERRQ(ierr);

    vars.erase("climatic_mass_balance");
  }

  ierr = SurfaceModel::write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism

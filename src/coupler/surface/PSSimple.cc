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
#include <stdexcept>

namespace pism {

///// Simple PISM surface model.
PSSimple::PSSimple(IceGrid &g, const Config &conf)
  : SurfaceModel(g, conf),
    climatic_mass_balance(g.get_unit_system(), "climatic_mass_balance", grid),
    ice_surface_temp(g.get_unit_system(), "ice_surface_temp", grid) {

  climatic_mass_balance.set_string("pism_intent", "diagnostic");
  climatic_mass_balance.set_string("long_name",
                                   "surface mass balance (accumulation/ablation) rate");
  climatic_mass_balance.set_string("standard_name",
                                   "land_ice_surface_specific_mass_balance_flux");
  climatic_mass_balance.set_units("kg m-2 s-1");
  climatic_mass_balance.set_glaciological_units("kg m-2 year-1");

  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  ice_surface_temp.set_units("K");
}


void PSSimple::init(Vars &vars) {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  assert(atmosphere != NULL);
  atmosphere->init(vars);

  verbPrintf(2, grid.com,
             "* Initializing the simplest PISM surface (snow) processes model PSSimple.\n"
             "  It passes atmospheric state directly to upper ice fluid surface:\n"
             "    surface mass balance          := precipitation,\n"
             "    ice upper surface temperature := 2m air temperature.\n");
}

void PSSimple::update(double my_t, double my_dt) {
  m_t = my_t;
  m_dt = my_dt;
  if (atmosphere) {
    atmosphere->update(my_t, my_dt);
  }
}


void PSSimple::ice_surface_mass_flux(IceModelVec2S &result) {
  atmosphere->mean_precipitation(result);
  result.scale(config.get("ice_density"));
}

void PSSimple::ice_surface_temperature(IceModelVec2S &result) {
  atmosphere->mean_annual_temp(result);
}

void PSSimple::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  SurfaceModel::add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("ice_surface_temp");
    result.insert("climatic_mass_balance");
  }
}

void PSSimple::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {

  if (set_contains(vars, "ice_surface_temp")) {
    ice_surface_temp.define(nc, nctype, true);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    climatic_mass_balance.define(nc, nctype, true);
  }

  SurfaceModel::define_variables(vars, nc, nctype);
}

void PSSimple::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    tmp.create(grid, "ice_surface_temp", WITHOUT_GHOSTS);
    tmp.metadata() = ice_surface_temp;

    ice_surface_temperature(tmp);

    tmp.write(nc);

    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    IceModelVec2S tmp;
    tmp.create(grid, "climatic_mass_balance", WITHOUT_GHOSTS);
    tmp.metadata() = climatic_mass_balance;

    ice_surface_mass_flux(tmp);
    tmp.write_in_glaciological_units = true;
    tmp.write(nc);

    vars.erase("climatic_mass_balance");
  }

  SurfaceModel::write_variables(vars, nc);
}

} // end of namespace pism

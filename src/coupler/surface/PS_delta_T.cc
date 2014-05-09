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

#include "PS_delta_T.hh"
#include "PISMConfig.hh"

namespace pism {

/// -surface ...,delta_T (scalar forcing of ice surface temperatures)

PS_delta_T::PS_delta_T(IceGrid &g, const Config &conf, SurfaceModel* in)
  : PScalarForcing<SurfaceModel,PSModifier>(g, conf, in),
    climatic_mass_balance(g.get_unit_system()),
    ice_surface_temp(g.get_unit_system())
{
  PetscErrorCode ierr = allocate_PS_delta_T(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

PS_delta_T::~PS_delta_T() {
  // empty
}

PetscErrorCode PS_delta_T::allocate_PS_delta_T() {
  option_prefix = "-surface_delta_T";
  offset_name   = "delta_T";

  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));

  offset->get_metadata().set_units("Kelvin");
  offset->get_metadata().set_string("long_name", "ice-surface temperature offsets");
  offset->get_dimension_metadata().set_units(grid.time->units_string());

  climatic_mass_balance.init_2d("climatic_mass_balance", grid);
  climatic_mass_balance.set_string("pism_intent", "diagnostic");
  climatic_mass_balance.set_string("long_name",
                                   "surface mass balance (accumulation/ablation) rate");
  climatic_mass_balance.set_string("standard_name",
                                   "land_ice_surface_specific_mass_balance_flux");
  climatic_mass_balance.set_units("kg m-2 s-1");
  climatic_mass_balance.set_glaciological_units("kg m-2 year-1");

  ice_surface_temp.init_2d("ice_surface_temp", grid);
  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  ice_surface_temp.set_units("K");

  return 0;
}

PetscErrorCode PS_delta_T::init(Vars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing ice-surface temperature forcing using scalar offsets...\n"); CHKERRQ(ierr);

  ierr = init_internal(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PS_delta_T::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->ice_surface_temperature(result); CHKERRQ(ierr);
  ierr = offset_data(result); CHKERRQ(ierr);
  return 0;
}

void PS_delta_T::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("ice_surface_temp");
    result.insert("climatic_mass_balance");
  }
}

PetscErrorCode PS_delta_T::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    ierr = ice_surface_temp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = climatic_mass_balance.define(nc, nctype, true); CHKERRQ(ierr);
  }

  ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PS_delta_T::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
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

  ierr = input_model->write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism

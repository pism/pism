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

#include "PSLapseRates.hh"

PetscErrorCode PSLapseRates::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool smb_lapse_rate_set;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "  [using temperature and mass balance lapse corrections]\n"); CHKERRQ(ierr);

  ierr = init_internal(vars); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Lapse rate options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsReal("-smb_lapse_rate",
                           "Elevation lapse rate for the surface mass balance, in m/year per km",
                           smb_lapse_rate, smb_lapse_rate_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "   ice upper-surface temperature lapse rate: %3.3f K per km\n"
                    "   ice-equivalent surface mass balance lapse rate: %3.3f m/year per km\n",
                    temp_lapse_rate, smb_lapse_rate); CHKERRQ(ierr);

  temp_lapse_rate = convert(temp_lapse_rate, "K/km", "K/m");

  smb_lapse_rate = convert(smb_lapse_rate, "m/year / km", "m/s / m");

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

PetscErrorCode PSLapseRates::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = input_model->ice_surface_mass_flux(result); CHKERRQ(ierr);
  ierr = lapse_rate_correction(result, smb_lapse_rate); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PSLapseRates::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = input_model->ice_surface_temperature(result); CHKERRQ(ierr);
  ierr = lapse_rate_correction(result, temp_lapse_rate); CHKERRQ(ierr);
  return 0;
}

void PSLapseRates::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  if (keyword == "medium" || keyword == "big") {
    result["ice_surface_temp"] = ice_surface_temp;
    result["climatic_mass_balance"] = climatic_mass_balance;
  }

  input_model->add_vars_to_output(keyword, result);
}

PetscErrorCode PSLapseRates::define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype) {
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

PetscErrorCode PSLapseRates::write_variables(set<string> vars, string filename) {
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

  ierr = input_model->write_variables(vars, filename); CHKERRQ(ierr);

  return 0;
}


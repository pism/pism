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

#include "PSAnomaly.hh"
#include "IceGrid.hh"

PSAnomaly::PSAnomaly(IceGrid &g, const PISMConfig &conf, PISMSurfaceModel* in)
  : PGivenClimate<PSModifier,PISMSurfaceModel>(g, conf, in),
    climatic_mass_balance(g.get_unit_system()),
    ice_surface_temp(g.get_unit_system())
{
  PetscErrorCode ierr = allocate_PSAnomaly(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

PSAnomaly::~PSAnomaly() {
  // empty
}

PetscErrorCode PSAnomaly::allocate_PSAnomaly() {
  PetscErrorCode ierr;
  option_prefix  = "-surface_anomaly";

  // will be de-allocated by the parent's destructor
  climatic_mass_balance_anomaly = new IceModelVec2T;
  ice_surface_temp_anomaly      = new IceModelVec2T;

  m_fields["climatic_mass_balance_anomaly"] = climatic_mass_balance_anomaly;
  m_fields["ice_surface_temp_anomaly"] = ice_surface_temp_anomaly;

  ierr = process_options(); CHKERRQ(ierr);

  std::map<std::string, std::string> standard_names;
  ierr = set_vec_parameters(standard_names); CHKERRQ(ierr);

  ierr = ice_surface_temp_anomaly->create(grid, "ice_surface_temp_anomaly", false); CHKERRQ(ierr);
  ierr = climatic_mass_balance_anomaly->create(grid, "climatic_mass_balance_anomaly", false); CHKERRQ(ierr);

  ierr = ice_surface_temp_anomaly->set_attrs("climate_forcing",
                                             "anomaly of the temperature of the ice at the ice surface but below firn processes",
                                             "Kelvin", ""); CHKERRQ(ierr);
  ierr = climatic_mass_balance_anomaly->set_attrs("climate_forcing",
                                                  "anomaly of the surface mass balance (accumulation/ablation) rate",
                                                  "kg m-2 s-1", ""); CHKERRQ(ierr);
  ierr = climatic_mass_balance_anomaly->set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);
  climatic_mass_balance_anomaly->write_in_glaciological_units = true;

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

PetscErrorCode PSAnomaly::init(PISMVars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  if (input_model != NULL) {
    ierr = input_model->init(vars); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the '-surface ...,anomaly' modifier...\n"); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "    reading anomalies from %s ...\n", filename.c_str()); CHKERRQ(ierr);

  ierr = ice_surface_temp_anomaly->init(filename, bc_period, bc_reference_time); CHKERRQ(ierr);
  ierr = climatic_mass_balance_anomaly->init(filename, bc_period, bc_reference_time); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSAnomaly::update(double my_t, double my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = climatic_mass_balance_anomaly->average(m_t, m_dt); CHKERRQ(ierr);
  ierr = ice_surface_temp_anomaly->average(m_t, m_dt); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSAnomaly::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->ice_surface_mass_flux(result); CHKERRQ(ierr);
  ierr = result.add(1.0, *climatic_mass_balance_anomaly); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSAnomaly::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->ice_surface_temperature(result); CHKERRQ(ierr);
  ierr = result.add(1.0, *ice_surface_temp_anomaly); CHKERRQ(ierr);

  return 0;
}

void PSAnomaly::add_vars_to_output(std::string keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("ice_surface_temp");
    result.insert("climatic_mass_balance");
  }
}

PetscErrorCode PSAnomaly::define_variables(std::set<std::string> vars, const PIO &nc, PISM_IO_Type nctype) {
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

PetscErrorCode PSAnomaly::write_variables(std::set<std::string> vars, const PIO &nc) {
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

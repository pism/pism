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

#include "PSGivenClimate.hh"
#include "IceGrid.hh"

namespace pism {

PSGivenClimate::PSGivenClimate(IceGrid &g, const Config &conf)
  : PGivenClimate<PSModifier,SurfaceModel>(g, conf, NULL)
{
  option_prefix = "-surface_given";

  PetscErrorCode ierr = allocate_PSGivenClimate(); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    throw std::runtime_error("PSGivenClimate allocation failed");
  }

}

PSGivenClimate::~PSGivenClimate() {
  // empty
}

PetscErrorCode PSGivenClimate::allocate_PSGivenClimate() {
  PetscErrorCode ierr;

  ice_surface_temp      = new IceModelVec2T;
  climatic_mass_balance = new IceModelVec2T;

  m_fields["ice_surface_temp"]      = ice_surface_temp;
  m_fields["climatic_mass_balance"] = climatic_mass_balance;

  process_options();

  std::map<std::string, std::string> standard_names;
  standard_names["climatic_mass_balance"] = "land_ice_surface_specific_mass_balance_flux";
  set_vec_parameters(standard_names);

  ierr = ice_surface_temp->create(grid, "ice_surface_temp", false); CHKERRQ(ierr);
  ierr = climatic_mass_balance->create(grid, "climatic_mass_balance", false); CHKERRQ(ierr);

  ierr = ice_surface_temp->set_attrs("climate_forcing",
                                     "temperature of the ice at the ice surface but below firn processes",
                                     "Kelvin", ""); CHKERRQ(ierr);
  ierr = climatic_mass_balance->set_attrs("climate_forcing",
                                          "surface mass balance (accumulation/ablation) rate",
                                          "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux"); CHKERRQ(ierr);
  ierr = climatic_mass_balance->set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);
  climatic_mass_balance->write_in_glaciological_units = true;

  return 0;
}

void PSGivenClimate::attach_atmosphere_model(AtmosphereModel *input) {
  delete input;
  input = NULL;
}

void PSGivenClimate::init(Vars &) {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  verbPrintf(2, grid.com,
             "* Initializing the surface model reading temperature at the top of the ice\n"
             "  and ice surface mass flux from a file...\n");

  ice_surface_temp->init(filename, bc_period, bc_reference_time);
  climatic_mass_balance->init(filename, bc_period, bc_reference_time);

  // read time-independent data right away:
  if (ice_surface_temp->get_n_records() == 1 && climatic_mass_balance->get_n_records() == 1) {
    update(grid.time->current(), 0); // dt is irrelevant
  }
}

void PSGivenClimate::update(double my_t, double my_dt) {
  update_internal(my_t, my_dt);

  climatic_mass_balance->average(m_t, m_dt);
  ice_surface_temp->average(m_t, m_dt);
}

void PSGivenClimate::ice_surface_mass_flux(IceModelVec2S &result) {
  climatic_mass_balance->copy_to(result);
}

void PSGivenClimate::ice_surface_temperature(IceModelVec2S &result) {
  ice_surface_temp->copy_to(result);
}

} // end of namespace pism

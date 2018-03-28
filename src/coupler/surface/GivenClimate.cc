// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "GivenClimate.hh"
#include "pism/util/IceGrid.hh"

namespace pism {
namespace surface {

Given::Given(IceGrid::ConstPtr grid, std::shared_ptr<atmosphere::AtmosphereModel> input)
  : PGivenClimate<SurfaceModel>(grid, nullptr)
{
  (void) input;

  m_temperature = new IceModelVec2T;
  m_mass_flux   = new IceModelVec2T;

  m_fields["ice_surface_temp"]      = m_temperature;
  m_fields["climatic_mass_balance"] = m_mass_flux;

  process_options("-surface_given");

  set_vec_parameters(m_filename, {{"climatic_mass_balance", "land_ice_surface_specific_mass_balance_flux"}});

  m_temperature->set_attrs("climate_forcing",
                           "temperature of the ice at the ice surface but below firn processes",
                           "Kelvin", "");
  m_temperature->metadata().set_doubles("valid_range", {0.0, 323.15}); // [0C, 50C]

  const double smb_max = m_config->get_double("surface.given.smb_max", "kg m-2 second-1");

  m_mass_flux->set_attrs("climate_forcing",
                         "surface mass balance (accumulation/ablation) rate",
                         "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux");
  m_mass_flux->metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_mass_flux->metadata().set_doubles("valid_range", {-smb_max, smb_max});

  // this class uses the "modifier" version of the SurfaceModel constructor, so we need to
  // allocate these here
  {
    m_liquid_water_fraction = allocate_liquid_water_fraction(grid);
    m_layer_mass            = allocate_layer_mass(grid);
    m_layer_thickness       = allocate_layer_thickness(grid);

    // default values
    m_layer_thickness->set(0.0);
    m_layer_mass->set(0.0);
    m_liquid_water_fraction->set(0.0);
  }
}

Given::~Given() {
  // empty
}

void Given::init_impl(const Geometry &geometry) {

  m_log->message(2,
                 "* Initializing the surface model reading temperature at the top of the ice\n"
                 "  and ice surface mass flux from a file...\n");

  m_temperature->init(m_filename, m_bc_period, m_bc_reference_time);
  m_mass_flux->init(m_filename, m_bc_period, m_bc_reference_time);

  // read time-independent data right away:
  if (m_temperature->get_n_records() == 1 && m_mass_flux->get_n_records() == 1) {
    update(geometry, m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void Given::update_impl(const Geometry &geometry, double t, double dt) {
  update_internal(geometry, t, dt);

  m_mass_flux->average(t, dt);
  m_temperature->average(t, dt);
}

const IceModelVec2S &Given::mass_flux_impl() const {
  return *m_mass_flux;
}

const IceModelVec2S &Given::temperature_impl() const {
  return *m_temperature;
}

} // end of namespace surface
} // end of namespace pism

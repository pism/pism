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

#include "Runoff_SMB.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace ocean {


Runoff_SMB::Runoff_SMB(IceGrid::ConstPtr g, OceanModel* in)
  : PScalarForcing<OceanModel,OceanModifier>(g, in) {

  m_option_prefix = "-ocean_runoff_smb";
  m_offset_name = "delta_T";

  m_offset.reset(new Timeseries(*m_grid, m_offset_name, m_config->get_string("time.dimension_name")));
  m_offset->variable().set_string("units", "Kelvin");
  m_offset->variable().set_string("long_name", "air temperature offsets");
  m_offset->dimension().set_string("units", m_grid->ctx()->time()->units_string());

  m_temp_to_runoff_a       = m_config->get_double("surface.temp_to_runoff_a");
  m_runoff_to_ocean_melt_b = m_config->get_double("ocean.runoff_to_ocean_melt_b");

  m_runoff_to_ocean_melt_power_alpha = m_config->get_double("ocean.runoff_to_ocean_melt_power_alpha");
  m_runoff_to_ocean_melt_power_beta  = m_config->get_double("ocean.runoff_to_ocean_melt_power_beta");
}

Runoff_SMB::~Runoff_SMB() {
  // empty
}

void Runoff_SMB::init_impl() {

  m_input_model->init();

  m_log->message(2,
                 "* Initializing ice shelf base mass flux forcing using scalar multiplier...\n"
                 "*   derived from delta_T air temperature modifier\n");

  init_internal();
}

MaxTimestep Runoff_SMB::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("ocean runoff_SMB");
}

void Runoff_SMB::update_impl(double t, double dt) {
  super::update_impl(t, dt);

  mass_flux(m_shelf_base_mass_flux);
}

void Runoff_SMB::mass_flux(IceModelVec2S &result) const {
  // m_current_forcing is set by PScalarForcing::update_impl()
  double delta_T = m_current_forcing;

  // short-cuts, just to make the formula below easier to read
  double
    a     = m_temp_to_runoff_a,
    B     = m_runoff_to_ocean_melt_b,
    alpha = m_runoff_to_ocean_melt_power_alpha,
    beta = m_runoff_to_ocean_melt_power_beta,
    scale_factor = 1.0;

  /*
     This parameterization only works if delta_T > 0 because
     negative numbers cannot be raised to a fractional power
     so we do not scale the forcing if delta_T < 0
  */
  if (delta_T > 0.0) {
    scale_factor = 1 + B * pow(a * delta_T, alpha) * pow(delta_T, beta);
  }
    result.scale(scale_factor);
}

} // end of namespace ocean
} // end of namespace pism

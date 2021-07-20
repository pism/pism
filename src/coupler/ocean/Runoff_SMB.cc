// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021 PISM Authors
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
#include "pism/coupler/util/ScalarForcing.hh"

namespace pism {
namespace ocean {

Runoff_SMB::Runoff_SMB(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> in)
  : OceanModel(g, in) {

  m_forcing.reset(new ScalarForcing(g->ctx(),
                                    "-ocean_runoff_smb",
                                    "delta_T",
                                    "Kelvin",
                                    "Kelvin",
                                    "air temperature offsets"));

  m_temp_to_runoff_a       = m_config->get_number("surface.temp_to_runoff_a");
  m_runoff_to_ocean_melt_b = m_config->get_number("ocean.runoff_to_ocean_melt_b");

  m_runoff_to_ocean_melt_power_alpha = m_config->get_number("ocean.runoff_to_ocean_melt_power_alpha");
  m_runoff_to_ocean_melt_power_beta  = m_config->get_number("ocean.runoff_to_ocean_melt_power_beta");

  m_shelf_base_mass_flux = allocate_shelf_base_mass_flux(g);
}

Runoff_SMB::~Runoff_SMB() {
  // empty
}

void Runoff_SMB::init_impl(const Geometry &geometry) {

  m_input_model->init(geometry);

  m_log->message(2,
                 "* Initializing ice shelf base mass flux forcing using scalar multiplier...\n"
                 "*   derived from delta_T air temperature modifier\n");
}

void Runoff_SMB::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_forcing->update(t, dt);

  mass_flux(m_forcing->value(), *m_shelf_base_mass_flux);
}

void Runoff_SMB::mass_flux(double delta_T, IceModelVec2S &result) const {

  // short-cuts, just to make the formula below easier to read
  double
    a            = m_temp_to_runoff_a,
    B            = m_runoff_to_ocean_melt_b,
    alpha        = m_runoff_to_ocean_melt_power_alpha,
    beta         = m_runoff_to_ocean_melt_power_beta,
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

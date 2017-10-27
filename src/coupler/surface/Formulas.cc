/* Copyright (C) 2014, 2015, 2016, 2017 PISM Authors
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

#include "Formulas.hh"
#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace surface {

PSFormulas::PSFormulas(IceGrid::ConstPtr g)
  : SurfaceModel(g) {
  m_climatic_mass_balance.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  m_climatic_mass_balance.set_attrs("internal",
                                    "ice-equivalent surface mass balance (accumulation/ablation) rate",
                                    "kg m-2 s-1",
                                    "land_ice_surface_specific_mass_balance_flux");
  m_climatic_mass_balance.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // annual mean air temperature at "ice surface", at level below all
  // firn processes (e.g. "10 m" or ice temperatures)
  m_ice_surface_temp.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
  m_ice_surface_temp.set_attrs("internal",
                               "annual average ice surface temperature, below firn processes",
                               "K", "");
}

PSFormulas::~PSFormulas() {
  // empty
}


void PSFormulas::attach_atmosphere_model_impl(atmosphere::AtmosphereModel *input) {
  delete input;
}

void PSFormulas::mass_flux_impl(IceModelVec2S &result) const {
  result.copy_from(m_climatic_mass_balance);
}

void PSFormulas::temperature_impl(IceModelVec2S &result) const {
  result.copy_from(m_ice_surface_temp);
}

void PSFormulas::define_model_state_impl(const PIO &output) const {
  // these are *not* model state, but I want to be able to re-start from a file produced using this
  // class
  m_climatic_mass_balance.define(output);
  m_ice_surface_temp.define(output);
}

void PSFormulas::write_model_state_impl(const PIO &output) const {
  // these are *not* model state, but I want to be able to re-start from a file produced using this
  // class
  m_climatic_mass_balance.write(output);
  m_ice_surface_temp.write(output);
}

} // end of namespace surface
} // end of namespace pism

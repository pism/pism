/* Copyright (C) 2014, 2015, 2016, 2017, 2018 PISM Authors
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

PSFormulas::PSFormulas(IceGrid::ConstPtr grid)
  : SurfaceModel(grid) {

  m_mass_flux   = allocate_mass_flux(grid);
  m_temperature = allocate_temperature(grid);
}

PSFormulas::~PSFormulas() {
  // empty
}

const IceModelVec2S &PSFormulas::mass_flux_impl() const {
  return *m_mass_flux;
}

const IceModelVec2S & PSFormulas::temperature_impl() const {
  return *m_temperature;
}

void PSFormulas::define_model_state_impl(const PIO &output) const {
  // these are *not* model state, but I want to be able to re-start from a file produced using this
  // class
  m_mass_flux->define(output);
  m_temperature->define(output);
}

void PSFormulas::write_model_state_impl(const PIO &output) const {
  // these are *not* model state, but I want to be able to re-start from a file produced using this
  // class
  m_mass_flux->write(output);
  m_temperature->write(output);
}

} // end of namespace surface
} // end of namespace pism

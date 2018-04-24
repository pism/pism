/* Copyright (C) 2018 PISM Authors
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

#include "CompleteOceanModel.hh"

namespace pism {
namespace ocean {

// "modifier" constructor
CompleteOceanModel::CompleteOceanModel(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> input)
  : OceanModel(g, input) {

  m_shelf_base_temperature = allocate_shelf_base_temperature(g);
  m_shelf_base_mass_flux = allocate_shelf_base_mass_flux(g);
}

// "model" constructor
CompleteOceanModel::CompleteOceanModel(IceGrid::ConstPtr g)
  : CompleteOceanModel(g, nullptr) {
  // empty
}

CompleteOceanModel::~CompleteOceanModel() {
  // empty
}

const IceModelVec2S& CompleteOceanModel::shelf_base_temperature_impl() const {
  return *m_shelf_base_temperature;
}

const IceModelVec2S& CompleteOceanModel::shelf_base_mass_flux_impl() const {
  return *m_shelf_base_mass_flux;
}

} // end of namespace ocean
} // end of namespace pism

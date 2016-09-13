/* Copyright (C) 2015, 2016 PISM Authors
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

#include "PISMYieldStress.hh"
#include "base/util/PISMConfigInterface.hh"

namespace pism {

YieldStress::YieldStress(IceGrid::ConstPtr g)
  : Component(g) {
  m_basal_yield_stress.create(m_grid, "tauc", WITH_GHOSTS,
                              m_config->get_double("grid.max_stencil_width"));
  // PROPOSED standard_name = land_ice_basal_material_yield_stress
  m_basal_yield_stress.set_attrs("model_state",
                                 "yield stress for basal till (plastic or pseudo-plastic model)",
                                 "Pa", "");
}

YieldStress::~YieldStress() {
  // empty
}

void YieldStress::init() {
  this->init_impl();
}

void YieldStress::update() {
  this->update_impl();
}

const IceModelVec2S& YieldStress::basal_material_yield_stress() {
  return m_basal_yield_stress;
}

} // end of namespace pism

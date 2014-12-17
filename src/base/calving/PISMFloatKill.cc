/* Copyright (C) 2013, 2014 PISM Authors
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

#include "PISMFloatKill.hh"
#include "Mask.hh"
#include "iceModelVec.hh"

namespace pism {

FloatKill::FloatKill(IceGrid &g)
  : Component(g) {
  // empty
}

FloatKill::~FloatKill() {
  // empty
}

void FloatKill::init(Vars &/*vars*/) {
  verbPrintf(2, m_grid.com,
             "* Initializing the 'calving at the grounding line' mechanism...\n");
}

/**
 * Updates ice cover mask and the ice thickness using the calving rule
 * removing all floating ice.
 *
 * @param[in,out] pism_mask ice cover (cell type) mask
 * @param[in,out] ice_thickness ice thickness
 *
 * @return 0 on success
 */
void FloatKill::update(IceModelVec2Int &pism_mask, IceModelVec2S &ice_thickness) {
  MaskQuery M(pism_mask);

  IceModelVec::AccessList list;
  list.add(pism_mask);
  list.add(ice_thickness);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (M.floating_ice(i, j)) {
      ice_thickness(i, j) = 0.0;
      pism_mask(i, j)     = MASK_ICE_FREE_OCEAN;
    }
  }

  pism_mask.update_ghosts();
  ice_thickness.update_ghosts();
}

void FloatKill::add_vars_to_output(const std::string &/*keyword*/,
                                       std::set<std::string> &/*result*/) {
  // empty
}

void FloatKill::define_variables(const std::set<std::string> &/*vars*/, const PIO &/*nc*/,
                                               IO_Type /*nctype*/) {
  // empty
}

void FloatKill::write_variables(const std::set<std::string> &/*vars*/, const PIO& /*nc*/) {
  // empty
}

} // end of namespace pism

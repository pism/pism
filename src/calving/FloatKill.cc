/* Copyright (C) 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#include "FloatKill.hh"

#include "pism/util/Mask.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace calving {

FloatKill::FloatKill(IceGrid::ConstPtr g)
  : Component(g) {
  m_margin_only = m_config->get_boolean("calving.float_kill.margin_only");
  m_calve_near_grounding_line = m_config->get_boolean("calving.float_kill.calve_near_grounding_line");
}

FloatKill::~FloatKill() {
  // empty
}

void FloatKill::init() {
  m_log->message(2,
                 "* Initializing calving using the floatation criterion (float_kill)...\n");

  if (m_margin_only) {
    m_log->message(2,
                   "  [only cells at the ice margin are calved during a given time step]\n");
  }

  if (not m_calve_near_grounding_line) {
    m_log->message(2,
                   "  [keeping floating cells near the grounding line]\n");
  }
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
void FloatKill::update(IceModelVec2CellType &mask, IceModelVec2S &ice_thickness) {

  IceModelVec::AccessList list{&mask, &ice_thickness};

  const bool dont_calve_near_grounded_ice = not m_calve_near_grounding_line;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.floating_ice(i, j)) {
      if (m_margin_only and not mask.next_to_ice_free_ocean(i, j)) {
        continue;
      }

      if (dont_calve_near_grounded_ice and mask.next_to_grounded_ice(i, j)) {
        continue;
      }

      ice_thickness(i, j) = 0.0;
      mask(i, j)          = MASK_ICE_FREE_OCEAN;
    }
  }

  mask.update_ghosts();
  ice_thickness.update_ghosts();
}

} // end of namespace calving
} // end of namespace pism

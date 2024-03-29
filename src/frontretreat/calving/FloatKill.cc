/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2022, 2023 PISM Authors
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

#include "pism/frontretreat/calving/FloatKill.hh"

#include "pism/util/Mask.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/CellType.hh"

namespace pism {
namespace calving {

FloatKill::FloatKill(std::shared_ptr<const Grid> g)
  : Component(g),
    m_old_mask(m_grid, "old_mask") {
  m_margin_only = m_config->get_flag("calving.float_kill.margin_only");
  m_calve_near_grounding_line = m_config->get_flag("calving.float_kill.calve_near_grounding_line");
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
void FloatKill::update(array::Scalar &cell_type, array::Scalar &ice_thickness) {

  // this call fills ghosts of m_old_mask
  m_old_mask.copy_from(cell_type);

  array::AccessScope list{&cell_type, &m_old_mask, &ice_thickness};

  const bool dont_calve_near_grounded_ice = not m_calve_near_grounding_line;

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_old_mask.floating_ice(i, j)) {
      if (m_margin_only and not m_old_mask.next_to_ice_free_ocean(i, j)) {
        continue;
      }

      if (dont_calve_near_grounded_ice and m_old_mask.next_to_grounded_ice(i, j)) {
        continue;
      }

      ice_thickness(i, j) = 0.0;
      cell_type(i, j)     = MASK_ICE_FREE_OCEAN;
    }
  }

  cell_type.update_ghosts();
  ice_thickness.update_ghosts();
}

} // end of namespace calving
} // end of namespace pism

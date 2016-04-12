/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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
#include "base/util/Mask.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_const.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {
namespace calving {

FloatKill::FloatKill(IceGrid::ConstPtr g)
  : Component(g) {
  m_margin_only = m_config->get_boolean("float_kill_margin_only");
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

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(ice_thickness);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.floating_ice(i, j)) {
      if (m_margin_only and not mask.next_to_ice_free_ocean(i, j)) {
        continue;
      }
      ice_thickness(i, j) = 0.0;
      mask(i, j)          = MASK_ICE_FREE_OCEAN;
    }
  }

  mask.update_ghosts();
  ice_thickness.update_ghosts();
}

void FloatKill::add_vars_to_output_impl(const std::string &/*keyword*/,
                                       std::set<std::string> &/*result*/) {
  // empty
}

void FloatKill::define_variables_impl(const std::set<std::string> &/*vars*/, const PIO &/*nc*/,
                                               IO_Type /*nctype*/) {
  // empty
}

void FloatKill::write_variables_impl(const std::set<std::string> &/*vars*/, const PIO& /*nc*/) {
  // empty
}

} // end of namespace calving
} // end of namespace pism

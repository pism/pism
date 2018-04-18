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

#include "WeertmanSliding.hh"

#include "pism/rheology/FlowLawFactory.hh"
#include "pism/geometry/Geometry.hh"
#include "StressBalance.hh"

namespace pism {
namespace stressbalance {

WeertmanSliding::WeertmanSliding(IceGrid::ConstPtr grid)
  : ShallowStressBalance(grid) {
  // Use the SIA flow law.
  rheology::FlowLawFactory ice_factory("stress_balance.sia.", m_config, m_EC);
  m_flow_law = ice_factory.create();
}

WeertmanSliding::~WeertmanSliding() {
  // empty
}

void WeertmanSliding::init_impl() {
  m_log->message(2, "* Initializing Weertman-style basal sliding...\n");
}

void WeertmanSliding::update(const Inputs &inputs, bool full_update) {

  (void) full_update;

  const IceModelVec2S &bed_elevation = inputs.geometry->bed_elevation;

  IceModelVec::AccessList list{&m_velocity};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      m_velocity(i, j) = {0.0, 0.0};
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

}

} // end of namespace stressbalance
} // end of namespace pism

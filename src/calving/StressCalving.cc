/* Copyright (C) 2016, 2017, 2018 PISM Authors
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

#include "StressCalving.hh"

namespace pism {
namespace calving {

StressCalving::StressCalving(IceGrid::ConstPtr grid,
                             unsigned int stencil_width)
  // mask has to have a wider stencil to be able to call is_marginal() away from the current cell
  : CalvingFrontRetreat(grid, stencil_width + 1),
    m_stencil_width(stencil_width) {

  m_strain_rates.create(m_grid, "strain_rates", WITH_GHOSTS,
                        m_stencil_width,
                        2);     // 2 components

  m_strain_rates.metadata(0).set_name("eigen1");
  m_strain_rates.set_attrs("internal",
                           "major principal component of horizontal strain-rate",
                           "second-1", "", 0);

  m_strain_rates.metadata(1).set_name("eigen2");
  m_strain_rates.set_attrs("internal",
                           "minor principal component of horizontal strain-rate",
                           "second-1", "", 1);
}

StressCalving::~StressCalving() {
  // empty
}

} // end of namespace calving
} // end of namespace pism

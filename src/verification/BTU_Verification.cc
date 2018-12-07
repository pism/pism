/* Copyright (C) 2016, 2017 PISM Authors
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

#include "BTU_Verification.hh"
#include "pism/util/Time.hh"
#include "pism/verification/tests/exactTestK.h"
#include "pism/verification/tests/exactTestO.h"
#include "pism/util/error_handling.hh"

namespace pism {
namespace energy {

BTU_Verification::BTU_Verification(IceGrid::ConstPtr g,
                                   const BTUGrid &vertical_grid,
                                   int testname, bool bedrock_is_ice)
  : BTU_Full(g, vertical_grid) {

  m_testname       = testname;
  m_bedrock_is_ice = bedrock_is_ice;
}

BTU_Verification::~BTU_Verification() {
  // empty
}

void BTU_Verification::initialize_bottom_surface_flux() {
  // hard-wired value used in exact solutions (tests K and O)
  m_bottom_surface_flux.set(0.042);
}

void BTU_Verification::bootstrap(const IceModelVec2S &bedrock_top_temperature) {
  (void) bedrock_top_temperature;

  std::vector<double> Tbcol(m_Mbz),
    zlevels = m_temp.levels();

  double time = m_grid->ctx()->time()->current();

  // evaluate exact solution in a column; all columns are the same
  switch (m_testname) {
  default:
  case 'K':
    for (unsigned int k = 0; k < m_Mbz; k++) {
      TestKParameters P = exactK(time, zlevels[k], m_bedrock_is_ice);
      if (P.error_code != 0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "exactK() reports that level %9.7f is below B0 = -1000.0 m",
                                      zlevels[k]);
      }
      Tbcol[k] = P.T;
    }
    break;
  case 'O':
    for (unsigned int k = 0; k < m_Mbz; k++) {
      Tbcol[k] = exactO(zlevels[k]).TT;
    }
    break;
  }

  // copy column values into 3D arrays
  IceModelVec::AccessList list(m_temp);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      m_temp.set_column(p.i(), p.j(), &Tbcol[0]);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

} // end of namespace energy
} // end of namespace pism

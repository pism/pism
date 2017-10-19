/* Copyright (C) 2015, 2016, 2017 PISM Authors
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

#include "SIAFD_Regional.hh"
#include "pism/stressbalance/StressBalance.hh"


namespace pism {

namespace stressbalance {

SIAFD_Regional::SIAFD_Regional(IceGrid::ConstPtr g)
  : SIAFD(g) {
  // empty
}

SIAFD_Regional::~SIAFD_Regional() {
  // empty
}

void SIAFD_Regional::init() {

  SIAFD::init();

  m_log->message(2, "  using the regional version of the SIA solver...\n");
}

void SIAFD_Regional::compute_surface_gradient(const Inputs &inputs,
                                              IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) const {

  SIAFD::compute_surface_gradient(inputs, h_x, h_y);

  const IceModelVec2Int &nmm = *inputs.no_model_mask;
  const IceModelVec2S &hst = *inputs.no_model_surface_elevation;

  const int Mx = m_grid->Mx(), My = m_grid->My();
  const double dx = m_grid->dx(), dy = m_grid->dy();  // convenience

  IceModelVec::AccessList list{&h_x, &h_y, &nmm, &hst};

  for (PointsWithGhosts p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // x-component, i-offset
    if (nmm(i, j) > 0.5 || nmm(i + 1, j) > 0.5) {

      if (i < 0 || i + 1 > Mx - 1) {
        h_x(i, j, 0) = 0.0;
      } else {
        h_x(i, j, 0) = (hst(i + 1, j) - hst(i, j)) / dx;
      }
    }

    // x-component, j-offset
    if (nmm(i - 1, j + 1) > 0.5 || nmm(i + 1, j + 1) > 0.5 ||
        nmm(i - 1, j)     > 0.5 || nmm(i + 1, j)     > 0.5) {

      if (i - 1 < 0 || j + 1 > My - 1 || i + 1 > Mx - 1) {
        h_x(i, j, 1) = 0.0;
      } else {
        h_x(i, j, 1) = ( + hst(i + 1, j + 1) + hst(i + 1, j)
                         - hst(i - 1, j + 1) - hst(i - 1, j)) / (4.0 * dx);
      }

    }

    // y-component, i-offset
    if (nmm(i, j + 1) > 0.5 || nmm(i + 1, j + 1) > 0.5 ||
        nmm(i, j - 1) > 0.5 || nmm(i + 1, j - 1) > 0.5) {
      if (i < 0 || j + 1 > My - 1 || i + 1 > Mx - 1 || j - 1 < 0) {
        h_y(i, j, 0) = 0.0;
      } else {
        h_y(i, j, 0) = ( + hst(i + 1, j + 1) + hst(i, j + 1)
                         - hst(i + 1, j - 1) - hst(i, j - 1)) / (4.0 * dy);
      }
    }

    // y-component, j-offset
    if (nmm(i, j) > 0.5 || nmm(i, j + 1) > 0.5) {

      if (j < 0 || j + 1 > My - 1) {
        h_y(i, j, 1) = 0.0;
      } else {
        h_y(i, j, 1) = (hst(i, j + 1) - hst(i, j)) / dy;
      }
    }

  }
}

} // end of namespace stressbalance

} // end of namespace pism

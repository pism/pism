/* Copyright (C) 2015, 2016, 2017, 2020 PISM Authors
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

#include "SSAFD_Regional.hh"
#include "pism/util/Vars.hh"
#include "pism/stressbalance/StressBalance.hh"

namespace pism {

namespace stressbalance {

SSAFD_Regional::SSAFD_Regional(IceGrid::ConstPtr g)
  : SSAFD(g) {

  m_h_stored      = NULL;
  m_H_stored      = NULL;
  m_no_model_mask = NULL;
}

void SSAFD_Regional::init() {

  SSAFD::init();

  m_log->message(2, "  using the regional version of the SSA solver...\n");

  if (m_config->get_flag("stress_balance.ssa.dirichlet_bc")) {
    m_log->message(2, "  using stored SSA velocities as Dirichlet B.C. in the no_model_strip...\n");
  }
}

void SSAFD_Regional::update(const Inputs &inputs, bool full_update) {
  m_h_stored      = inputs.no_model_surface_elevation;
  m_H_stored      = inputs.no_model_ice_thickness;
  m_no_model_mask = inputs.no_model_mask;

  SSA::update(inputs, full_update);

  m_h_stored      = NULL;
  m_H_stored      = NULL;
  m_no_model_mask = NULL;
}

static int weight(int M_ij, int M_n, double h_ij, double h_n) {
  // fjord walls, nunataks, headwalls
  if ((mask::icy(M_ij) and mask::ice_free(M_n) and h_n > h_ij) or
      (mask::ice_free(M_ij) and mask::icy(M_n) and h_ij > h_n)) {
    return 0;
  }

  return 1;
}

void SSAFD_Regional::compute_driving_stress(const array::Scalar &ice_thickness,
                                            const array::Scalar &surface_elevation,
                                            const array::CellType1 &cell_type,
                                            const array::Scalar *no_model_mask,
                                            IceModelVec2V &result) const {

  SSAFD::compute_driving_stress(ice_thickness, surface_elevation, cell_type, no_model_mask, result);

  double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  int
    Mx = m_grid->Mx(),
    My = m_grid->My();

  IceModelVec::AccessList list{&result, &cell_type, no_model_mask, m_h_stored, m_H_stored};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto M = no_model_mask->star(i, j);

    if (M.ij == 0) {
      // this grid point is in the modeled area so we don't need to modify the driving
      // stress
      continue;
    }

    double pressure = m_EC->pressure((*m_H_stored)(i, j));
    if (pressure <= 0) {
      result(i, j) = 0.0;
      continue;
    }

    auto h = m_h_stored->star(i, j);
    auto CT = cell_type.star(i, j);

    // x-derivative
    double h_x = 0.0;
    {
      double
        west = M.w == 1 and i > 0,
        east = M.e == 1 and i < Mx - 1;

      // don't use differences spanning "cliffs"
      west *= weight(CT.ij, CT.w, h.ij, h.w);
      east *= weight(CT.ij, CT.e, h.ij, h.e);

      if (east + west > 0) {
        h_x = 1.0 / ((west + east) * dx) * (west * (h.ij - h.w) + east * (h.e - h.ij));
      } else {
        h_x = 0.0;
      }
    }

    // y-derivative
    double h_y = 0.0;
    {
      double
        south = M.s == 1 and j > 0,
        north = M.n == 1 and j < My - 1;

      // don't use differences spanning "cliffs"
      south *= weight(CT.ij, CT.s, h.ij, h.s);
      north *= weight(CT.ij, CT.n, h.ij, h.n);

      if (north + south > 0) {
        h_y = 1.0 / ((south + north) * dy) * (south * (h.ij - h.s) + north * (h.n - h.ij));
      } else {
        h_y = 0.0;
      }
    }

    result(i, j) = - pressure * Vector2(h_x, h_y);
  } // end of the loop over grid points
}

SSA * SSAFD_RegionalFactory(IceGrid::ConstPtr grid) {
  return new SSAFD_Regional(grid);
}

} // end of namespace stressbalance

} // end of namespace pism

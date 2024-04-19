/* Copyright (C) 2015, 2016, 2017, 2019, 2022, 2023, 2024 PISM Authors
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

#include "pism/regional/SIAFD_Regional.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {

namespace stressbalance {

SIAFD_Regional::SIAFD_Regional(std::shared_ptr<const Grid> grid)
  : SIAFD(grid),
    m_h_x_no_model(grid, "h_x_no_model"),
    m_h_y_no_model(grid, "h_y_no_model") {
  // empty
}

void SIAFD_Regional::init() {

  SIAFD::init();

  m_log->message(2, "  using the regional version of the SIA solver...\n");
}

void SIAFD_Regional::compute_surface_gradient(const Inputs &inputs,
                                              array::Staggered1 &h_x,
                                              array::Staggered1 &h_y) {

  SIAFD::compute_surface_gradient(inputs, h_x, h_y);

  if (inputs.no_model_surface_elevation == nullptr) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "no_model surface elevation was not provided to SIAFD_Regional");
  }

  // this call updates ghosts of h_x_no_model and h_y_no_model
  surface_gradient_haseloff(*inputs.no_model_surface_elevation,
                            inputs.geometry->cell_type,
                            m_h_x_no_model, m_h_y_no_model);

  const array::Scalar2 &no_model = *inputs.no_model_mask;

  const int Mx = m_grid->Mx(), My = m_grid->My();

  array::AccessScope list{&h_x, &h_y, &no_model, &m_h_x_no_model, &m_h_y_no_model};

  for (auto p = m_grid->points(1); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto M = no_model.box(i, j);

    // x-component, i-offset
    if (M.c > 0.5 or M.e > 0.5) {

      if (i < 0 or i + 1 > Mx - 1) {
        h_x(i, j, 0) = 0.0;
      } else {
        h_x(i, j, 0) = m_h_x_no_model(i, j, 0);
      }
    }

    // x-component, j-offset
    if (M.nw > 0.5 or M.ne > 0.5 or M.w  > 0.5 or M.e  > 0.5) {

      if (i - 1 < 0 or j + 1 > My - 1 or i + 1 > Mx - 1) {
        h_x(i, j, 1) = 0.0;
      } else {
        h_x(i, j, 1) = m_h_x_no_model(i, j, 1);
      }

    }

    // y-component, i-offset
    if (M.n > 0.5 or M.ne > 0.5 or M.s > 0.5 or M.se > 0.5) {

      if (i < 0 or j + 1 > My - 1 or i + 1 > Mx - 1 or j - 1 < 0) {
        h_y(i, j, 0) = 0.0;
      } else {
        h_y(i, j, 0) = m_h_y_no_model(i, j, 0);
      }

    }

    // y-component, j-offset
    if (M.c > 0.5 or M.n > 0.5) {

      if (j < 0 or j + 1 > My - 1) {
        h_y(i, j, 1) = 0.0;
      } else {
        h_y(i, j, 1) = m_h_y_no_model(i, j, 1);
      }

    }
  } // end of the loop over grid points
}

} // end of namespace stressbalance

} // end of namespace pism

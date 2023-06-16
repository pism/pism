/* Copyright (C) 2020, 2022 PISM Authors
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

#include "BlatterTestXY.hh"

#include "pism/rheology/IsothermalGlen.hh"

#include "manufactured_solutions.hh"

namespace pism {
namespace stressbalance {

BlatterTestXY::BlatterTestXY(std::shared_ptr<const Grid> grid, int Mz, int coarsening_factor)
  : Blatter(grid, Mz, coarsening_factor) {

  // use the isothermal Glen flow law
  m_flow_law.reset(new rheology::IsothermalGlen("stress_balance.blatter.", *m_config, m_EC));

  // make sure we use the same Glen exponent
  assert(m_flow_law->exponent() == 3.0);

  // store constant hardness (enthalpy and pressure values are irrelevant)
  m_B = m_flow_law->hardness(1e5, 0);
}

bool BlatterTestXY::marine_boundary(int face,
                                    const int *node_type,
                                    const double *ice_bottom,
                                    const double *sea_level) {
  (void) face;
  (void) node_type;
  (void) ice_bottom;
  (void) sea_level;

  return false;
}

bool BlatterTestXY::dirichlet_node(const DMDALocalInfo &info, const fem::Element3::GlobalIndex& I) {
  // use Dirichlet BC on the whole map-plane boundary
  return (I.i == 0 or I.i == info.mx - 1 or
          I.j == 0 or I.j == info.my - 1);
}

Vector2d BlatterTestXY::u_bc(double x, double y, double z) const {
  (void) z;

  return blatter_xy_exact(x, y);
}

void BlatterTestXY::residual_source_term(const fem::Q1Element3 &element,
                                         const double *surface,
                                         const double *bed,
                                         Vector2d *residual) {
  (void) surface;
  (void) bed;

  // compute x and y coordinates of quadrature points
  double
    *x = m_work[0],
    *y = m_work[1];
  {
    double
      *x_nodal = m_work[2],
      *y_nodal = m_work[3];

    for (int n = 0; n < element.n_chi(); ++n) {
      x_nodal[n] = element.x(n);
      y_nodal[n] = element.y(n);
    }

    element.evaluate(x_nodal, x);
    element.evaluate(y_nodal, y);
  }

  // loop over all quadrature points
  for (int q = 0; q < element.n_pts(); ++q) {
    auto W = element.weight(q) / m_scaling;

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      const auto &psi = element.chi(q, t);

      residual[t] += W * psi.val * blatter_xy_source(x[q], y[q], m_B);
    }
  }
}

} // end of namespace stressbalance
} // end of namespace pism

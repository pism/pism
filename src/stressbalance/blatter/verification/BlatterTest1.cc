/* Copyright (C) 2020 PISM Authors
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

#include "BlatterTest1.hh"

#include "pism/rheology/IsothermalGlen.hh"

#include "blatter_mms.hh"

namespace pism {
namespace stressbalance {

BlatterTest1::BlatterTest1(IceGrid::ConstPtr grid, int Mz, int n_levels, int coarsening_factor)
  : Blatter(grid, Mz, n_levels, coarsening_factor) {

  // use the isothermal Glen flow law
  m_flow_law.reset(new rheology::IsothermalGlen("stress_balance.blatter.", *m_config, m_EC));

  // make sure we use the same Glen exponent
  assert(m_flow_law->exponent() == 3.0);

  // store constant hardness (enthalpy and pressure values are irrelevant)
  m_B = m_flow_law->hardness(1e5, 0);
}

bool BlatterTest1::neumann_bc_face(int face, const int *node_type) {
  (void) face;
  (void) node_type;

  return false;
}

bool BlatterTest1::dirichlet_node(const DMDALocalInfo &info, const fem::Element3::GlobalIndex& I) {
  // use Dirichlet BC on the whole map-plane boundary
  return (I.i == 0 or I.i == info.mx - 1 or
          I.j == 0 or I.j == info.my - 1);
}

Vector2 BlatterTest1::u_bc(double x, double y, double z) {
  (void) z;

  return exact_xy(x, y);
}

void BlatterTest1::residual_source_term(const fem::Q1Element3 &element,
                                        const double *surface,
                                        Vector2 *residual) {
  (void) surface;

  // compute x and y coordinates of quadrature points
  const int Nq = element.n_pts();
  std::vector<double> x(Nq), y(Nq);
  {
    const int Nk = fem::q13d::n_chi;
    double x_nodal[Nk], y_nodal[Nk];

    for (int n = 0; n < Nk; ++n) {
      x_nodal[n] = element.x(n);
      y_nodal[n] = element.y(n);
    }

    element.evaluate(x_nodal, x.data());
    element.evaluate(y_nodal, y.data());
  }

  // loop over all quadrature points
  for (int q = 0; q < Nq; ++q) {
    auto W = element.weight(q);

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      const auto &psi = element.chi(q, t);

      residual[t] += W * psi.val * source_xy(x[q], y[q], m_B);
    }
  }
}

} // end of namespace stressbalance
} // end of namespace pism

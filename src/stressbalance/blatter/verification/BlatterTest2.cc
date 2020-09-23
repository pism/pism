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

#include "BlatterTest2.hh"

#include "pism/rheology/IsothermalGlen.hh"

#include "manufactured_solutions.hh"

namespace pism {
namespace stressbalance {

BlatterTest2::BlatterTest2(IceGrid::ConstPtr grid, int Mz, int n_levels, int coarsening_factor)
  : Blatter(grid, Mz, n_levels, coarsening_factor) {

  // use the isothermal Glen flow law
  m_flow_law.reset(new rheology::IsothermalGlen("stress_balance.blatter.", *m_config, m_EC));

  // make sure we use the same Glen exponent
  assert(m_flow_law->exponent() == 3.0);

  // make sure the grid is periodic in the Y direction
  assert(m_grid->periodicity() == Y_PERIODIC);

  // store constant hardness (enthalpy and pressure values are irrelevant)
  m_B = m_flow_law->hardness(1e5, 0);

  // surface elevation parameter (1/m)
  m_alpha = 4e-8;

  // surface elevation parameter (m)
  m_s0 = 2000.0;

  // ice density times g
  m_rhog = (m_config->get_number("constants.ice.density") *
            m_config->get_number("constants.standard_gravity"));

  // ice thickness (m)
  m_H0 = 1000.0;
}

bool BlatterTest2::neumann_bc_face(int face, const int *node_type) {
  (void) face;
  (void) node_type;

  return false;
}

bool BlatterTest2::dirichlet_node(const DMDALocalInfo &info, const fem::Element3::GlobalIndex& I) {
  // use Dirichlet BC on the whole map-plane boundary
  return (I.i == 0 or I.i == info.mx - 1);
}

Vector2 BlatterTest2::u_bc(double x, double y, double z) {
  (void) y;

  return exact_xz(x, z, m_B, m_rhog, m_s0, m_alpha, m_H0);
}

void BlatterTest2::residual_source_term(const fem::Q1Element3 &element,
                                        const double *surface,
                                        Vector2 *residual) {
  (void) surface;

  // compute x and y coordinates of quadrature points
  const int Nq = element.n_pts();
  std::vector<double> x(Nq), z(Nq);
  {
    const int Nk = fem::q13d::n_chi;
    double x_nodal[Nk], z_nodal[Nk];

    for (int n = 0; n < Nk; ++n) {
      x_nodal[n] = element.x(n);
      z_nodal[n] = element.z(n);
    }

    element.evaluate(x_nodal, x.data());
    element.evaluate(z_nodal, z.data());
  }

  // loop over all quadrature points
  for (int q = 0; q < Nq; ++q) {
    auto W = element.weight(q);

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      const auto &psi = element.chi(q, t);

      residual[t] += W * psi.val * source_xz(x[q], z[q], m_B, m_rhog, m_s0, m_alpha, m_H0);
    }
  }
}

} // end of namespace stressbalance
} // end of namespace pism

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

#include "BlatterTestHalfar.hh"

#include "pism/rheology/IsothermalGlen.hh"

#include "manufactured_solutions.hh"

namespace pism {
namespace stressbalance {

BlatterTestHalfar::BlatterTestHalfar(IceGrid::ConstPtr grid, int Mz, int n_levels, int coarsening_factor)
  : Blatter(grid, Mz, n_levels, coarsening_factor) {

  // use the isothermal Glen flow law
  m_flow_law.reset(new rheology::IsothermalGlen("stress_balance.blatter.", *m_config, m_EC));

  // make sure we use the same Glen exponent
  assert(m_flow_law->exponent() == 3.0);

  // make sure the grid is periodic in the Y direction
  assert(m_grid->periodicity() == Y_PERIODIC);

  // store constant ice softness (enthalpy and pressure values are irrelevant)
  m_B = m_flow_law->hardness(1e5, 0);

  m_R0 = 2.0 * grid->Lx();

  // ice thickness (m)
  m_H0 = 3600.0;

  m_rho = m_config->get_number("constants.ice.density");

  m_g   = m_config->get_number("constants.standard_gravity");
}

bool BlatterTestHalfar::dirichlet_node(const DMDALocalInfo &info, const fem::Element3::GlobalIndex& I) {
  (void) info;
  // use Dirichlet BC at x == 0
  return I.i == 0;
}

Vector2 BlatterTestHalfar::u_bc(double x, double y, double z) const {
  (void) y;

  return blatter_xz_halfar_exact(x, z, m_H0, m_R0, m_rho, m_g, m_B);
}

void BlatterTestHalfar::residual_source_term(const fem::Q1Element3 &element,
                                             const double *surface,
                                             const double *bed,
                                             Vector2 *residual) {

  (void) surface;
  (void) bed;

  // compute x and z coordinates of quadrature points
  double
    *x = m_work[0],
    *z = m_work[1];
  {
    double
      *x_nodal = m_work[2],
      *z_nodal = m_work[3];

    for (int n = 0; n < fem::q13d::n_chi; ++n) {
      x_nodal[n] = element.x(n);
      z_nodal[n] = element.z(n);
    }

    element.evaluate(x_nodal, x);
    element.evaluate(z_nodal, z);
  }

  // loop over all quadrature points
  for (int q = 0; q < element.n_pts(); ++q) {
    auto W = element.weight(q);

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      const auto &psi = element.chi(q, t);

      residual[t] += W * psi.val * blatter_xz_halfar_source(x[q], z[q], m_H0, m_R0, m_rho, m_g, m_B);
    }
  }
}

void BlatterTestHalfar::residual_lateral(const fem::Q1Element3 &element,
                                         const fem::Q1Element3Face &face,
                                         const double *surface_nodal,
                                         const double *z_nodal,
                                         const double *sl_nodal,
                                         Vector2 *residual) {
  (void) surface_nodal;
  (void) sl_nodal;

  // compute x and z coordinates of quadrature points
  double
    *x = m_work[0],
    *z = m_work[1];
  {
    double
      *x_nodal = m_work[2];

    for (int n = 0; n < fem::q13d::n_chi; ++n) {
      x_nodal[n] = element.x(n);
    }

    face.evaluate(x_nodal, x);
    face.evaluate(z_nodal, z);
  }

  // loop over all quadrature points
  for (int q = 0; q < face.n_pts(); ++q) {
    auto W = face.weight(q);

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      auto psi = face.chi(q, t);

      Vector2 F(0.0, 0.0);
      if (x[q] > 0.0) {
        // this lateral BC does not apply at the left (x = 0) boundary.
        F = blatter_xz_halfar_lateral(x[q], z[q], m_H0, m_R0, m_rho, m_g, m_B);
      }

      // N = (1, 0, 0) is implied ("right" boundary)
      residual[t] += W * psi * F;
    }
  }
}

} // end of namespace stressbalance
} // end of namespace pism

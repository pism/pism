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

#include "BlatterTestXZ.hh"

#include "pism/rheology/IsothermalGlen.hh"

#include "manufactured_solutions.hh"

namespace pism {
namespace stressbalance {

BlatterTestXZ::BlatterTestXZ(IceGrid::ConstPtr grid,
                             int Mz, int coarsening_factor)
  : Blatter(grid, Mz, coarsening_factor) {

  // use the isothermal Glen flow law
  m_flow_law.reset(new rheology::IsothermalGlen("stress_balance.blatter.", *m_config, m_EC));

  // make sure we use the same Glen exponent
  assert(m_flow_law->exponent() == 3.0);

  // make sure the grid is periodic in the Y direction
  assert(m_grid->periodicity() == Y_PERIODIC);

  // store constant ice softness (enthalpy and pressure values are irrelevant)
  m_A = m_flow_law->softness(1e5, 0);

  // surface elevation parameter (1/m)
  m_alpha = 4e-8;

  // surface elevation parameter (m)
  m_s0 = 2000.0;

  // ice thickness (m)
  m_H0 = 1000.0;

  // sliding parameter
  m_beta = convert(m_sys, 1.0, "kPa year m-1", "Pa s m-1");

  m_rho = m_config->get_number("constants.ice.density");
  m_g   = m_config->get_number("constants.standard_gravity");
}

bool BlatterTestXZ::marine_boundary(int face,
                                    const int *node_type,
                                    const double *ice_bottom,
                                    const double *sea_level) {
  (void) face;
  (void) node_type;
  (void) ice_bottom;
  (void) sea_level;

  return false;
}

bool BlatterTestXZ::dirichlet_node(const DMDALocalInfo &info,
                                   const fem::Element3::GlobalIndex& I) {
  // use Dirichlet BC at x == -Lx and x == Lx
  return (I.i == 0 or I.i == info.mx - 1);
}

Vector2 BlatterTestXZ::u_bc(double x, double y, double z) const {
  (void) y;

  return blatter_xz_exact(x, z, m_A, m_rho, m_g, m_s0, m_alpha, m_H0, m_beta);
}

void BlatterTestXZ::residual_source_term(const fem::Q1Element3 &element,
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

      residual[t] += W * psi.val * blatter_xz_source(x[q], z[q], m_A, m_rho, m_g,
                                                     m_s0, m_alpha, m_H0, m_beta);
    }
  }
}

/*!
 * Basal contribution to the residual.
 */
void BlatterTestXZ::residual_basal(const fem::Q1Element3 &element,
                                   const fem::Q1Element3Face &face,
                                   const double *tauc_nodal,
                                   const double *f_nodal,
                                   const Vector2 *u_nodal,
                                   Vector2 *residual) {

  // The basal sliding BC contribution:
  Blatter::residual_basal(element, face, tauc_nodal, f_nodal, u_nodal, residual);

  // Additional contribution needed to satisfy the manufactured solution:
  {
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

      face.evaluate(x_nodal, x);
      face.evaluate(z_nodal, z);
    }

    for (int q = 0; q < face.n_pts(); ++q) {
      auto W = face.weight(q);

      // loop over all test functions
      for (int t = 0; t < element.n_chi(); ++t) {
        auto psi = face.chi(q, t);

        residual[t] += - W * psi * blatter_xz_source_bed(x[q], z[q], m_A, m_rho, m_g,
                                                         m_s0, m_alpha, m_H0, m_beta);
      }
    }
  }
}

/*!
 * Top surface contribution to the residual.
 */
void BlatterTestXZ::residual_surface(const fem::Q1Element3 &element,
                                     const fem::Q1Element3Face &face,
                                     Vector2 *residual) {

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

    face.evaluate(x_nodal, x);
    face.evaluate(z_nodal, z);
  }

  for (int q = 0; q < face.n_pts(); ++q) {
    auto W = face.weight(q);

    auto F = blatter_xz_source_surface(x[q], z[q], m_A, m_rho, m_g,
                                       m_s0, m_alpha, m_H0, m_beta);

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      auto psi = face.chi(q, t);

      residual[t] += - W * psi * F;
    }
  }
}

} // end of namespace stressbalance
} // end of namespace pism

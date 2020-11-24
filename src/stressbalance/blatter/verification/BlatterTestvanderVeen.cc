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

#include "BlatterTestvanderVeen.hh"

#include "pism/rheology/IsothermalGlen.hh"

#include "manufactured_solutions.hh"

namespace pism {
namespace stressbalance {

BlatterTestvanderVeen::BlatterTestvanderVeen(IceGrid::ConstPtr grid,
                                             int Mz, int n_levels, int coarsening_factor)
  : Blatter(grid, Mz, n_levels, coarsening_factor) {

  // use the isothermal Glen flow law
  m_flow_law.reset(new rheology::IsothermalGlen("stress_balance.blatter.", *m_config, m_EC));

  // make sure we use the same Glen exponent
  assert(m_flow_law->exponent() == 3.0);

  // make sure the grid is periodic in the Y direction
  assert(m_grid->periodicity() == Y_PERIODIC);

  // store constant ice hardness (enthalpy and pressure values are irrelevant)
  m_B = m_flow_law->hardness(1e5, 0);

  // ice thickness (m)
  m_H0 = 100.0;

  // ice velocity (m / s)
  m_V0 = convert(m_sys, 100.0, "m / year", "m / s");

  m_rho_i = m_config->get_number("constants.ice.density");

  m_rho_w = m_config->get_number("constants.sea_water.density");

  m_g = m_config->get_number("constants.standard_gravity");

  {
    double alpha = 1.0 - m_rho_i / m_rho_w;
    m_C0 = pow(alpha * m_g * m_rho_i / (2.0 * m_B), 3.0);
  }
}

bool BlatterTestvanderVeen::dirichlet_node(const DMDALocalInfo &info,
                                           const fem::Element3::GlobalIndex& I) {
  (void) info;
  // use Dirichlet BC at x == 0
  return I.i == 0 or I.i == info.mx - 1; /*  FIXME: remove Dirichlet BC at the right end point */
}

Vector2 BlatterTestvanderVeen::u_bc(double x, double y, double z) const {
  (void) y;
  (void) z;

  return u_exact(x);
}

Vector2 BlatterTestvanderVeen::u_exact(double x) {
  double Q0 = m_H0 * m_V0;
  return { Q0 / H_exact(x), 0.0 };
}

double BlatterTestvanderVeen::H_exact(double x) {
  double Q0 = m_H0 * m_V0;
  return pow(4 * m_C0 * x / Q0 + pow(m_H0, -4), -1.0 / 4.0);
}

double BlatterTestvanderVeen::beta_exact(double x) {
  double
    H     = H_exact(x),
    Q0    = m_H0 * m_V0,
    alpha = 1.0 - m_rho_i / m_rho_w;

  return 2 * m_B * pow(m_C0, 4.0 / 3.0) * (1 - alpha) * pow(H, 7) / pow(Q0, 2);
}

void BlatterTestvanderVeen::residual_lateral(const fem::Q1Element3 &element,
                                             const fem::Q1Element3Face &face,
                                             const double *surface_nodal,
                                             const double *z_nodal,
                                             const double *sl_nodal,
                                             Vector2 *residual) {
  (void) surface_nodal;
  (void) sl_nodal;
  (void) z_nodal;

  // compute x coordinates of quadrature points
  double *x = m_work[0];
  {
    double *x_nodal = m_work[1];

    for (int n = 0; n < element.n_chi(); ++n) {
      x_nodal[n] = element.x(n);
    }

    face.evaluate(x_nodal, x);
  }

  // loop over all quadrature points
  for (int q = 0; q < face.n_pts(); ++q) {
    auto W = face.weight(q);
    auto N3 = face.normal(q);
    Vector2 N = {N3.x, N3.y};

    double F = m_g * m_rho_i * (1.0 - m_rho_i / m_rho_w) * H_exact(x[q]);

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      auto psi = face.chi(q, t);

      residual[t] += - W * psi * F * N;
    }
  }
}

} // end of namespace stressbalance
} // end of namespace pism

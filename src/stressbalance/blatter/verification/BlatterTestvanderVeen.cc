/* Copyright (C) 2020, 2021, 2022 PISM Authors
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
#include "pism/util/node_types.hh"

#include "pism/rheology/IsothermalGlen.hh"

#include "manufactured_solutions.hh"

namespace pism {
namespace stressbalance {

BlatterTestvanderVeen::BlatterTestvanderVeen(std::shared_ptr<const IceGrid> grid,
                                             int Mz, int coarsening_factor)
  : Blatter(grid, Mz, coarsening_factor) {

  // use the isothermal Glen flow law
  m_flow_law.reset(new rheology::IsothermalGlen("stress_balance.blatter.", *m_config, m_EC));

  // make sure we use the same Glen exponent
  assert(m_flow_law->exponent() == 3.0);

  // make sure the grid is periodic in the Y direction
  assert(m_grid->periodicity() == Y_PERIODIC);

  // store constant ice hardness (enthalpy and pressure values are irrelevant)
  m_B = m_flow_law->hardness(1e5, 0);

  // ice thickness (m)
  m_H0 = 600.0;

  // ice velocity (m / s)
  m_V0 = convert(m_sys, 300.0, "m / year", "m / s");

  m_rho_ice = m_config->get_number("constants.ice.density");

  // s = alpha * H, where s is surface elevation and H ice thickness
  double rho_w = m_config->get_number("constants.sea_water.density");
  m_alpha = 1.0 - m_rho_ice / rho_w;

  m_g = m_config->get_number("constants.standard_gravity");
}

bool BlatterTestvanderVeen::dirichlet_node(const DMDALocalInfo &info,
                                           const fem::Element3::GlobalIndex& I) {
  (void) info;
  // use Dirichlet BC at x == 0
  return I.i == 0;
}

bool BlatterTestvanderVeen::marine_boundary(int face,
                                            const int *node_type,
                                            const double *ice_bottom,
                                            const double *sea_level) {
  // Ignore sea level and ice bottom elevation:
  (void) ice_bottom;
  (void) sea_level;

  auto nodes = fem::q13d::incident_nodes[face];

  // number of nodes per face
  int N = 4;

  // exclude faces that contain at least one node that is not a part of the boundary
  for (int n = 0; n < N; ++n) {
    if (not (node_type[nodes[n]] == NODE_BOUNDARY)) {
      return false;
    }
  }

  return true;
}

Vector2d BlatterTestvanderVeen::u_bc(double x, double y, double z) const {
  (void) y;
  (void) z;

  return u_exact(x);
}

Vector2d BlatterTestvanderVeen::u_exact(double x) const {
  double Q0 = m_H0 * m_V0;
  return blatter_xz_vanderveen_exact(x, m_alpha, m_H0, Q0, m_rho_ice, m_g, m_B);
}

double BlatterTestvanderVeen::H_exact(double x) const {
  double Q0 = m_H0 * m_V0;
  return blatter_xz_vanderveen_thickness(x, m_alpha, m_H0, Q0, m_rho_ice, m_g, m_B);
}

double BlatterTestvanderVeen::b_exact(double x) const {
  return (m_alpha - 1.0) * H_exact(x);
}

double BlatterTestvanderVeen::beta_exact(double x) const {
  double Q0 = m_H0 * m_V0;
  return blatter_xz_vanderveen_beta(x, m_alpha, m_H0, Q0, m_rho_ice, m_g, m_B);
}

void BlatterTestvanderVeen::residual_lateral(const fem::Q1Element3 &element,
                                             const fem::Q1Element3Face &face,
                                             const double *surface_nodal,
                                             const double *z_nodal,
                                             const double *sl_nodal,
                                             Vector2d *residual) {
  (void) surface_nodal;
  (void) sl_nodal;
  (void) z_nodal;

  double Q0 = m_H0 * m_V0;

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
    auto W = face.weight(q) / m_scaling;

    // the normal direction (1, 0, 0) is hard-wired here
    auto F = blatter_xz_vanderveen_source_lateral(x[q], m_alpha, m_H0, Q0,
                                                  m_rho_ice, m_g, m_B);

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      auto psi = face.chi(q, t);

      residual[t] += - W * psi * F;
    }
  }
}

/*!
 * Top surface contribution to the residual.
 */
void BlatterTestvanderVeen::residual_surface(const fem::Q1Element3 &element,
                                             const fem::Q1Element3Face &face,
                                             Vector2d *residual) {
  double Q0 = m_H0 * m_V0;

  // compute x coordinates of quadrature points
  double
    *x = m_work[0];
  {
    double *x_nodal = m_work[1];

    for (int n = 0; n < fem::q13d::n_chi; ++n) {
      x_nodal[n] = element.x(n);
    }

    face.evaluate(x_nodal, x);
  }

  for (int q = 0; q < face.n_pts(); ++q) {
    auto W = face.weight(q) / m_scaling;

    auto F = blatter_xz_vanderveen_source_surface(x[q], m_alpha, m_H0, Q0,
                                                  m_rho_ice, m_g, m_B);

    // loop over all test functions
    for (int t = 0; t < element.n_chi(); ++t) {
      auto psi = face.chi(q, t);

      residual[t] += - W * psi * F;
    }
  }
}

} // end of namespace stressbalance
} // end of namespace pism

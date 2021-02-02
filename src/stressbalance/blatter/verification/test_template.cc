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

namespace pism {
namespace stressbalance {

TestTemplate::TestTemplate(IceGrid::ConstPtr grid, int Mz, int n_levels, int coarsening_factor)
  : Blatter(grid, Mz, n_levels, coarsening_factor) {
  // empty
}

bool TestTemplate::marine_boundary(int face,
                                   const int *node_type,
                                   const double *ice_bottom,
                                   const double *sea_level) {
  (void) face;
  (void) node_type;

  return false;
}

bool TestTemplate::dirichlet_node(const DMDALocalInfo &info, const fem::Element3::GlobalIndex& I) {
  (void) info;
  (void) I;
  return false;
}

Vector2 TestTemplate::u_bc(double x, double y, double z) {
  (void) x;
  (void) y;
  (void) z;

  return 0.0
}

void TestTemplate::residual_source_term(const fem::Q1Element3 &element,
                                        const double *surface,
                                        const double *bed,
                                        Vector2 *residual) {
  (void) element;
  (void) surface;
  (void) bed;
  (void) residual;
}

void TestTemplate::jacobian_f(const fem::Q1Element3 &element,
                              const Vector2 *u_nodal,
                              const double *B_nodal,
                              double K[2 * fem::q13d::n_chi][2 * fem::q13d::n_chi]) {
  // use the same system of PDEs
  Blatter::jacobian_f(element, u_nodal, B_nodal, K);
}

void TestTemplate::residual_f(const fem::Q1Element3 &element,
                              const Vector2 *u_nodal,
                              const double *B_nodal,
                              Vector2 *residual) {
  // use the same system of PDEs
  Blatter::residual_f(element, u_nodal, B_nodal, residual);
}

void TestTemplate::jacobian_basal(const fem::Q1Element3Face &face,
                                  const double *tauc_nodal,
                                  const double *f_nodal,
                                  const Vector2 *u_nodal,
                                  double K[2 * fem::q13d::n_chi][2 * fem::q13d::n_chi]) {
  (void) face;
  (void) tauc_nodal;
  (void) f_nodal;
  (void) u_nodal;
  (void) K;
  // empty: this verification case does not use Robin BC at the base
}

void TestTemplate::residual_basal(const fem::Q1Element3 &element,
                                  const fem::Q1Element3Face &face,
                                  const double *tauc_nodal,
                                  const double *f_nodal,
                                  const Vector2 *u_nodal,
                                  Vector2 *residual) {
  (void) element;
  (void) face;
  (void) tauc_nodal;
  (void) f_nodal;
  (void) u_nodal;
  (void) residual;
  // empty: this verification case does not use Robin BC at the base
}

void TestTemplate::residual_lateral(const fem::Q1Element3 &element,
                                    const fem::Q1Element3Face &face,
                                    const double *z_nodal,
                                    const double *sl_nodal,
                                    Vector2 *residual) {
  (void) element;
  (void) face;
  (void) z_nodal;
  (void) sl_nodal;
  (void) residual;
  // empty: this verification test does not use lateral BC
}

} // end of namespace stressbalance
} // end of namespace pism

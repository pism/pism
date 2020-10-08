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

#ifndef BLATTERTESTCFBC_H
#define BLATTERTESTCFBC_H

#include "pism/stressbalance/blatter/Blatter.hh"

namespace pism {
namespace stressbalance {

/*!
 * Implements Dirichlet BC for an X-Z flow line setup testing the implementation of lateral (ice margin) boundary conditions.
 */
class BlatterTestCFBC : public Blatter {
public:
  BlatterTestCFBC(IceGrid::ConstPtr grid, int Mz, int n_levels, int coarsening_factor);

private:
  void init_2d_parameters(const Inputs &inputs);

  bool dirichlet_node(const DMDALocalInfo &info, const fem::Element3::GlobalIndex& I);

  Vector2 u_bc(double x, double y, double z);

  void residual_surface(const fem::Q1Element3 &element,
                        const fem::Q1Element3Face &face,
                        Vector2 *residual);

  void residual_basal(const fem::Q1Element3 &element,
                      const fem::Q1Element3Face &face,
                      const double *tauc_nodal,
                      const double *f_nodal,
                      const Vector2 *u_nodal,
                      Vector2 *residual);

  void jacobian_basal(const fem::Q1Element3Face &face,
                      const double *tauc_nodal,
                      const double *f_nodal,
                      const Vector2 *u_nodal,
                      double K[2 * fem::q13d::n_chi][2 * fem::q13d::n_chi]);
  double m_B;
  double m_g;
  double m_rho_i;
  double m_rho_w;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* BLATTERTESTCFBC_H */

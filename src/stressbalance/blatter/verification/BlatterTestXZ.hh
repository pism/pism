/* Copyright (C) 2020, 2021 PISM Authors
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

#ifndef BLATTERTESTXZ_H
#define BLATTERTESTXZ_H

#include "pism/stressbalance/blatter/Blatter.hh"

namespace pism {
namespace stressbalance {

/*!
 * Implements Dirichlet BC and the source term for a verification test from Tezaur et al,
 * 2015 (), section 4.2.
 *
 * Domain: [-50km, 50km] * [-1, 1] * [b, s].
 *
 * Here b is the bed elevation, s is the surface elevation.
 *
 * s = b + H0, H0 = 1000m.
 *
 * Dirichlet BC are imposed at all the nodes along the lateral boundary.
 *
 * Natural boundary conditions are used on the top boundary.
 *
 * The basal boundary condition includes the sliding condition and a correction for the
 * chosen manufactured solution.
 *
 */
class BlatterTestXZ : public Blatter {
public:
  BlatterTestXZ(IceGrid::ConstPtr grid, int Mz, int coarsening_factor);

private:
  bool marine_boundary(int face,
                       const int *node_type,
                       const double *ice_bottom,
                       const double *sea_level);

  bool dirichlet_node(const DMDALocalInfo &info, const fem::Element3::GlobalIndex& I);

  Vector2 u_bc(double x, double y, double z) const;

  void residual_source_term(const fem::Q1Element3 &element,
                            const double *surface,
                            const double *bed,
                            Vector2 *residual);

  void residual_basal(const fem::Q1Element3 &element,
                      const fem::Q1Element3Face &face,
                      const double *tauc_nodal,
                      const double *f_nodal,
                      const Vector2 *u_nodal,
                      Vector2 *residual);

  void residual_surface(const fem::Q1Element3 &element,
                        const fem::Q1Element3Face &face,
                        Vector2 *residual);

  //! constant ice hardness
  double m_A;

  double m_alpha;

  double m_s0;

  double m_H0;

  double m_rho;

  double m_g;

  double m_beta;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* BLATTERTESTXZ_H */

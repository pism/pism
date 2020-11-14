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

#ifndef BLATTERTESTHALFAR_H
#define BLATTERTESTHALFAR_H

#include "pism/stressbalance/blatter/Blatter.hh"

namespace pism {
namespace stressbalance {

/*!
 * Implements the analytical source term and Dirichlet BC for the X-Z Halfar dome setup.
 */
class BlatterTestHalfar : public Blatter {
public:
  BlatterTestHalfar(IceGrid::ConstPtr grid, int Mz, int n_levels, int coarsening_factor);

private:
  bool vertical_cliff_face(int face, const int *node_type);

  bool dirichlet_node(const DMDALocalInfo &info, const fem::Element3::GlobalIndex& I);

  Vector2 u_bc(double x, double y, double z);

  void residual_source_term(const fem::Q1Element3 &element,
                            const double *surface,
                            Vector2 *residual);

  double m_B;

  double m_H0;

  double m_R0;

  double m_rho;

  double m_g;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* BLATTERTESTHALFAR_H */

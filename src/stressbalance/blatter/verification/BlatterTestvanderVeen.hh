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

#ifndef BLATTERTESTVANDERVEEN_H
#define BLATTERTESTVANDERVEEN_H

#include "pism/stressbalance/blatter/Blatter.hh"

namespace pism {
namespace stressbalance {

/*!
 * Implements the Dirichlet BC at the inflow boundary and the stress BC at the calving
 * front boundary for the X-Z van der Veen shelf setup.
 *
 * The basal drag coefficient beta is x-dependent but its values are passed in using the
 * basal yield stress input field.
 */
class BlatterTestvanderVeen : public Blatter {
public:
  BlatterTestvanderVeen(IceGrid::ConstPtr grid, int Mz, int n_levels, int coarsening_factor);

  // Exact ice velocity
  Vector2 u_exact(double x) const;

  // Exact ice thickness
  double H_exact(double x) const;

  // Exact bed elevation
  double b_exact(double x) const;

  // Exact sliding coefficient beta
  double beta_exact(double x) const;
private:
  bool dirichlet_node(const DMDALocalInfo &info, const fem::Element3::GlobalIndex& I);

  Vector2 u_bc(double x, double y, double z) const;

  void residual_lateral(const fem::Q1Element3 &element,
                        const fem::Q1Element3Face &face,
                        const double *surface_nodal,
                        const double *z_nodal,
                        const double *sl_nodal,
                        Vector2 *residual);

  // constant ice hardness
  double m_B;

  // thickness at the inflow boundary
  double m_H0;

  // velocity at the inflow boundary
  double m_V0;

  // ratio of surface elevation to ice thickness
  double m_alpha;

  // just a constant used to simplify formulas (see the code for the definition)
  double m_C0;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* BLATTERTESTVANDERVEEN_H */

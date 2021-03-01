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

#ifndef BLATTERISMIPHOM_H
#define BLATTERISMIPHOM_H


#include "pism/stressbalance/blatter/Blatter.hh"

namespace pism {
namespace stressbalance {

enum ISMIPHOMTest {HOM_A, HOM_B, HOM_C, HOM_D};

/*!
 * This class implements periodic geometry experiments from the ISMIP-HOM
 * inter-comparison.
 */
class BlatterISMIPHOM : public Blatter {
public:
  BlatterISMIPHOM(IceGrid::ConstPtr grid, int Mz, int coarsening_factor,
                  ISMIPHOMTest test);

protected:
  void nodal_parameter_values(const fem::Q1Element3 &element,
                              Parameters **P,
                              int i,
                              int j,
                              int *node_type,
                              double *bottom,
                              double *thickness,
                              double *surface,
                              double *sea_level) const;
  ISMIPHOMTest m_test;

  typedef double (*geometry)(double x, double y, double L);

  geometry m_b, m_s;

  // domain length scale
  double m_L;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* BLATTERISMIPHOM_H */

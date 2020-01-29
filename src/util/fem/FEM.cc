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

#include "FEM.hh"

namespace pism {
namespace fem {

namespace q0 {
/*!
 * Piecewise-constant shape functions.
 */
Germ chi(unsigned int k, const QuadPoint &pt) {
  assert(k < q0::n_chi);

  Germ result;

  if ((k == 0 and pt.xi <= 0.0 and pt.eta <= 0.0) or
      (k == 1 and pt.xi > 0.0 and pt.eta <= 0.0) or
      (k == 2 and pt.xi > 0.0 and pt.eta > 0.0) or
      (k == 3 and pt.xi <= 0.0 and pt.eta > 0.0)) {
    result.val = 1.0;
  } else {
    result.val = 0.0;
  }

  result.dx  = 0.0;
  result.dy  = 0.0;

  return result;
}
} // end of namespace q0

namespace q1 {

//! Q1 basis functions on the reference element with nodes (-1,-1), (1,-1), (1,1), (-1,1).
Germ chi(unsigned int k, const QuadPoint &pt) {
  // coordinates of reference element nodes
  static const double xi[n_chi]  = {-1.0,  1.0, 1.0, -1.0};
  static const double eta[n_chi] = {-1.0, -1.0, 1.0,  1.0};

  assert(k < q1::n_chi);

  return {0.25 * (1.0 + xi[k] * pt.xi) * (1.0 + eta[k] * pt.eta),
          0.25 * xi[k] * (1.0 + eta[k] * pt.eta),
          0.25 * eta[k] * (1.0 + xi[k] * pt.xi)};
}

} // end of namespace q1

namespace p1 {

//! P1 basis functions on the reference element with nodes (0,0), (1,0), (0,1).
Germ chi(unsigned int k, const QuadPoint &pt) {
  assert(k < q1::n_chi);

  switch (k) {
  default:
  case 0:
    return {1.0 - pt.xi - pt.eta, -1.0, -1.0};
  case 1:
    return {pt.xi, 1.0, 0.0};
  case 2:
    return {pt.eta, 0.0, 1.0};
 }
}

} // end of namespace p1

} // end of namespace fem
} // end of namespace pism

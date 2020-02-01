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
#include "pism/util/node_types.hh"

namespace pism {
namespace fem {

namespace linear {

//! Linear basis functions on the interval [-1, -1]
Germ chi(unsigned int k, const QuadPoint &pt) {
  // coordinates of reference element nodes
  static const double xis[n_chi]  = {-1.0,  1.0};

  assert(k < linear::n_chi);

  return {0.5 * (1.0 + xis[k] * pt.xi),
          0.5 * xis[k],
          0.0};                         // unused
}

} // end of namespace linear

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

ElementType element_type(int node_type[q1::n_chi]) {

  // number of exterior nodes in this element
  const int n_exterior_nodes = ((node_type[0] == NODE_EXTERIOR) +
                                (node_type[1] == NODE_EXTERIOR) +
                                (node_type[2] == NODE_EXTERIOR) +
                                (node_type[3] == NODE_EXTERIOR));

  // an element is a "Q1 interior" if all its nodes are interior or boundary
  if (n_exterior_nodes == 0) {
    return ELEMENT_Q;
  }

  if (n_exterior_nodes == 1) {
    // an element is a "P1 interior" if it has exactly 1 exterior node

    for (unsigned int k = 0; k < q1::n_chi; ++k) {
      // Consider a Q1 element with one exterior node and three interior (or boundary)
      // nodes. This is not an interior element by itself, but it contains an embedded P1
      // element that *is* interior and should contribute. We need to find which of the four
      // types of embedded P1 elements to use, but P1 elements are numbered using the node
      // at the right angle of the reference element, not the "missing" node. Here we map
      // "missing" nodes to "opposite" nodes, i.e. nodes used to choose P1 element types.
      if (node_type[k] == NODE_EXTERIOR) {
        return ElementType((k + 2) % 4);
      }
    }
  }

  return ELEMENT_EXTERIOR;
}

} // end of namespace fem
} // end of namespace pism

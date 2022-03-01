/* Copyright (C) 2016, 2017, 2018, 2020, 2021, 2022 PISM Authors
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

#include "node_types.hh"

#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/error_handling.hh"

namespace pism {

/**
   Identify node types of all the nodes (grid points).

   Uses width=1 ghosts of `ice_thickness` (a "box" stencil).

   A node is can be *interior*, *boundary*, or *exterior*.

   An element is considered *icy* if at least three of its nodes have ice thickness above
   `thickness_threshold`.

   A node is considered *interior* if all of the elements it belongs to are icy.

   A node is considered *boundary* if it is not interior and at least one element it belongs to is
   icy.

   A node is considered *exterior* if it is neither interior nor boundary.

   Now an element "face" (side) is a part of a boundary if and only if both of its nodes are
   boundary nodes.

Cell layout:
~~~
(i-1,j+1) +-------N--------+ (i+1,j+1)
          |       |        |
          |  NW   |   NE   |
          |       |        |
          W-----(i,j)------E
          |       |        |
          |  SW   |   SE   |
          |       |        |
(i-1,j-1) +-------S--------+ (i+1,j-1)
~~~
 */
void compute_node_types(const array::Scalar &ice_thickness,
                        double thickness_threshold,
                        array::Scalar &result) {

  IceGrid::ConstPtr grid = ice_thickness.grid();

  const double &H_min = thickness_threshold;

  array::AccessScope list{&ice_thickness, &result};

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      auto H = ice_thickness.box(i, j);

      // flags indicating whether the current node and its neighbors are "icy"
      stencils::Box<int> icy;

      icy.ij = static_cast<int>(H.ij >= H_min);
      icy.nw = static_cast<int>(H.nw >= H_min);
      icy.n  = static_cast<int>(H.n  >= H_min);
      icy.ne = static_cast<int>(H.ne >= H_min);
      icy.e  = static_cast<int>(H.e  >= H_min);
      icy.se = static_cast<int>(H.se >= H_min);
      icy.s  = static_cast<int>(H.s  >= H_min);
      icy.sw = static_cast<int>(H.sw >= H_min);
      icy.w  = static_cast<int>(H.w  >= H_min);

      // flags indicating whether neighboring elements are "icy" (an element is icy if at
      // least three of its nodes are icy)
      const bool
        ne_element_is_icy = (icy.ij + icy.e + icy.ne + icy.n) >= 3,
        nw_element_is_icy = (icy.ij + icy.n + icy.nw + icy.w) >= 3,
        sw_element_is_icy = (icy.ij + icy.w + icy.sw + icy.s) >= 3,
        se_element_is_icy = (icy.ij + icy.s + icy.se + icy.e) >= 3;

      if (ne_element_is_icy and nw_element_is_icy and
          sw_element_is_icy and se_element_is_icy) {
        // all four elements are icy: we are at an interior node
        result(i, j) = NODE_INTERIOR;
      } else if (icy.ij != 0) {
        // the current node is icy: we are at a boundary
        result(i, j) = NODE_BOUNDARY;
      } else {
        // all elements are ice-free: we are at an exterior node
        result(i, j) = NODE_EXTERIOR;
      }

    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result.update_ghosts();
}

} // end of namespace pism

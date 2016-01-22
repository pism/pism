/* Copyright (C) 2016 PISM Authors
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

enum NodeType {
  NODE_INTERIOR = 0,
  NODE_BOUNDARY = 1,
  NODE_EXTERIOR = 2
};

/**
   Identify node types of all the nodes (grid points).

   Uses width=1 ghosts of `ice_thickness` (a "box" stencil).

   An element if either interior, a boundary, or exterior.

   FIXME: we need to use temporary storage and scatter ghosts to it instead of expecting a ghosted
   ice_thickness.

   Now an element side is a part of a boundary if and only if both of its nodes are on a boundary.


Cell layout:
~~~
+-------N--------+
|       |        |
|  NW   |   NE   |
|       |        |
W-----(i,j)------E
|       |        |
|  SW   |   SE   |
|       |        |
+-------S--------+
~~~
 */
void compute_node_types(const IceModelVec2S &ice_thickness,
                        double thickness_threshold,
                        IceModelVec2Int &result) {

  IceGrid::ConstPtr grid = ice_thickness.get_grid();

  const IceModelVec2S &H     = ice_thickness;
  const double        &H_min = thickness_threshold;

  IceModelVec::AccessList list(H);
  list.add(result);

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.();

      // indexing shortcuts (to reduce chances of making typos below)
      const int
        N = j + 1,
        E = i + 1,
        S = j - 1,
        W = i - 1;

      // booleans indicating whether the current node and its neighbors are "icy"
      const bool
        current_is_icy = H(i, j) >= H_min,
        nw_is_icy      = H(W, N) >= H_min,
        n_is_icy       = H(i, N) >= H_min,
        ne_is_icy      = H(E, N) >= H_min,
        e_is_icy       = H(E, j) >= H_min,
        se_is_icy      = H(E, S) >= H_min,
        s_is_icy       = H(i, S) >= H_min,
        sw_is_icy      = H(W, S) >= H_min,
        w_is_icy       = H(W, j) >= H_min;

      // booleans indicating whether neighboring elements are "icy" (and element is icy if all its
      // nodes are icy)
      const bool
        ne_elt_is_icy = (current_is_icy and e_is_icy and ne_is_icy and n_is_icy),
        nw_elt_is_icy = (current_is_icy and n_is_icy and nw_is_icy and w_is_icy),
        sw_elt_is_icy = (current_is_icy and w_is_icy and sw_is_icy and s_is_icy),
        se_elt_is_icy = (current_is_icy and s_is_icy and se_is_icy and e_is_icy);

      if (ne_elt_is_icy and nw_elt_is_icy and sw_elt_is_icy and se_elt_is_icy) {
        // all four elements are icy: we are at an interior node
        result(i, j) = NODE_INTERIOR;
      } else if (ne_elt_is_icy or nw_elt_is_icy or sw_elt_is_icy or se_elt_is_icy) {
        // at least one element is icy: we are at a boundary node
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
}

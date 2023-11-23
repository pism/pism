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
#ifndef PISM_ELEMENTITERATOR_H
#define PISM_ELEMENTITERATOR_H

namespace pism {

class Grid;

namespace fem {

//! Manages iterating over element indices.
/*! When computing residuals and Jacobians, there is a loop over all the elements in the
  Grid, and computations are done on each element. The Grid has an underlying PETSc
  DM, and a process typically does not own all of the nodes in the grid. Therefore we
  should perform a computation on a subset of the elements. In general, an element will
  have ghost (o) and real (*) vertices:

  \verbatim
  o---*---*---*---o
  |   |   |   |   |
  o---*---*---*---o
  |   |   |   |   |
  o---o---o---o---o
  \endverbatim

  The strategy is to perform computations on this process on every element that has
  a vertex that is owned by this processor.  But we only update entries in the
  global residual and Jacobian matrices if the corresponding row corresponds to a
  vertex that we own.  In the worst case, where each vertex of an element is owned by
  a different processor, the computations for that element will be repeated four times,
  once for each processor.

  This same strategy also correctly deals with periodic boundary conditions. The way PETSc
  deals with periodic boundaries can be thought of as using a kind of ghost. So the rule
  still works: compute on all elements containg a real vertex, but only update rows
  corresponding to that real vertex.

  The calculation of what elements to index over needs to account for ghosts and the
  presence or absense of periodic boundary conditions in the Grid. The ElementIterator
  performs that computation for you (see ElementIterator::xs and friends).
*/
class ElementIterator {
public:
  ElementIterator(const Grid &g);

  /*!\brief The total number of elements to be iterated over.  Useful for creating per-element storage.*/
  int element_count() {
    return xm*ym;
  }

  /*!\brief Convert an element index (`i`,`j`) into a flattened (1-d) array index, with the first
    element (`i`, `j`) to be iterated over corresponding to flattened index 0. */
  int flatten(int i, int j) {
    return (i-xs) + (j-ys)*xm;
  }

  //! x-coordinate of the first element to loop over.
  int xs;
  //! total number of elements to loop over in the x-direction.
  int xm;

  //! y-coordinate of the first element to loop over.
  int ys;
  //! total number of elements to loop over in the y-direction.
  int ym;

  //! x-index of the first local element.
  int lxs;
  //! total number local elements in x direction.
  int lxm;

  //! y-index of the first local element.
  int lys;
  //! total number local elements in y direction.
  int lym;
};

} // end of namespace fem
} // end of namespace pism

#endif /* PISM_ELEMENTITERATOR_H */

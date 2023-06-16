/* Copyright (C) 2020, 2023 PISM Authors
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

#include "ElementIterator.hh"

#include "pism/util/IceGrid.hh"

namespace pism {
namespace fem {

ElementIterator::ElementIterator(const IceGrid &grid) {
  // Start by assuming ghost elements exist in all directions.
  // Elements are indexed by their lower left vertex.  If there is a ghost
  // element on the right, its i-index will be the same as the maximum
  // i-index of a non-ghost vertex in the local grid.
  xs = grid.xs() - 1;                    // Start at ghost to the left.
  int xf = grid.xs() + grid.xm() - 1; // End at ghost to the right.
  ys = grid.ys() - 1;                    // Start at ghost at the bottom.
  int yf = grid.ys() + grid.ym() - 1; // End at ghost at the top.

  lxs = grid.xs();
  int lxf = lxs + grid.xm() - 1;
  lys = grid.ys();
  int lyf = lys + grid.ym() - 1;

  // Now correct if needed. The only way there will not be ghosts is if the
  // grid is not periodic and we are up against the grid boundary.

  if (!(grid.periodicity() & grid::X_PERIODIC)) {
    // Leftmost element has x-index 0.
    if (xs < 0) {
      xs = 0;
    }
    // Rightmost vertex has index grid.Mx-1, so the rightmost element has index grid.Mx-2
    if (xf > (int)grid.Mx() - 2) {
      xf  = grid.Mx() - 2;
      lxf = grid.Mx() - 2;
    }
  }

  if (!(grid.periodicity() & grid::Y_PERIODIC)) {
    // Bottom element has y-index 0.
    if (ys < 0) {
      ys = 0;
    }
    // Topmost vertex has index grid.My - 1, so the topmost element has index grid.My - 2
    if (yf > (int)grid.My() - 2) {
      yf  = grid.My() - 2;
      lyf = grid.My() - 2;
    }
  }

  // Tally up the number of elements in each direction
  xm  = xf - xs + 1;
  ym  = yf - ys + 1;
  lxm = lxf - lxs + 1;
  lym = lyf - lys + 1;
}

} // end of namespace fem
} // end of namespace pism

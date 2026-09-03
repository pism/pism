// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef PISM_DEBRIS_TERMINUS_KERNELS_HH
#define PISM_DEBRIS_TERMINUS_KERNELS_HH

#include <functional>
#include <vector>

namespace pism {
namespace debris {

/*!
 * Grid-walking helpers for the off-glacier debris removal scheme (Verhaegen and
 * Huybrechts, 2026, equation 23).
 *
 * Fields are passed as callables `f(i, j)` so that these functions are free of PISM, PETSc
 * and MPI and can be unit-tested with a plain executable.
 */
namespace terminus {

//! Offset of a cell relative to a starting cell.
struct Offset {
  int di, dj;
};

/*!
 * Walk up to `n` steps from `(i, j)`, each time moving to the ice-covered neighbor (of
 * eight) in the direction of the steepest ascent, i.e. with the largest increase of the
 * surface per unit distance. Stops early at a local maximum.
 *
 * This selects the cells "upstream (in the direction of the steepest slope uphill)" over
 * which the debris thickness and surface slope are averaged when the marginal length scale
 * exceeds the grid spacing.
 *
 * @return offsets (relative to `(i, j)`) of the cells visited, in order
 */
std::vector<Offset> upstream_chain(const std::function<double(int, int)> &surface,
                                   const std::function<bool(int, int)> &icy,
                                   int i, int j, int n, double dx, double dy);

/*!
 * Magnitude of the surface slope from the ice-covered cell `(i, j)` down into the ice-free
 * foreland: one-sided differences toward ice-free neighbors that are lower than `(i, j)`,
 * combined as `sqrt(s_x^2 + s_y^2)` (largest drop in each direction).
 *
 * @return 0 if no ice-free neighbor is lower than `(i, j)`
 */
double foreland_slope(const std::function<double(int, int)> &surface,
                      const std::function<bool(int, int)> &ice_free,
                      int i, int j, double dx, double dy);

} // end of namespace terminus
} // end of namespace debris
} // end of namespace pism

#endif /* PISM_DEBRIS_TERMINUS_KERNELS_HH */

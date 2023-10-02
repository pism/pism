// Copyright (C) 2007-2020, 2023 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cstddef>              // size_t
#include <cstring>
#include <cstdlib>
#include <algorithm>            // std::min
#include <gsl/gsl_interp.h>

#include "pism/util/io/File.hh"
#include "pism/util/io/LocalInterpCtx.hh"
#include "pism/util/Grid.hh"

#include "pism/util/interpolation.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Context.hh"

namespace pism {

//! Compute start and count for getting a subset of x.
/*! Given a grid `x` we find `x_start` and `x_count` so that `[subset_x_min, subset_x_max]` is in
 *  `[x[x_start], x[x_start + x_count]]`.
 *
 *  `x_start` and `x_count` define the smallest subset of `x` with this property.
 *
 *  Note that `x_start + x_count <= x.size()` as long as `x` is strictly increasing.
 *
 * @param[in] x input grid (defined interpolation domain)
 * @param[in] subset_x_min minimum `x` of a subset (interpolation range)
 * @param[in] subset_x_max maxumum `x` of a subset (interpolation range)
 * @param[out] x_start starting index
 * @param[out] x_count number of elements required
 */
static void subset_start_and_count(const std::vector<double> &x,
                                   double subset_x_min, double subset_x_max,
                                   int &x_start, int &x_count) {
  auto x_size = (int)x.size();

  x_start = (int)gsl_interp_bsearch(x.data(), subset_x_min, 0, x_size - 1);

  auto x_end = (int)gsl_interp_bsearch(x.data(), subset_x_max, 0, x_size - 1) + 1;

  x_end = std::min(x_size - 1, x_end);

  x_count = x_end - x_start + 1;
}

//! Construct a local interpolation context.
/*!
  The essential quantities to compute are where each processor should start within the NetCDF file grid
  (`start[]`) and how many grid points, from the starting place, the processor has.  The resulting
  portion of the grid is stored in array `a` (a field of the `LocalInterpCtx`).

  We make conservative choices about `start[]` and `count[]`.  In particular, the portions owned by
  processors \e must overlap at one point in the NetCDF file grid, but they \e may overlap more than that
  (as computed here).

  Note this constructor doesn't extract new information from the NetCDF file or
  do communication. The information from the NetCDF file must already be
  extracted, validly stored in a grid_info structure `input`.

  The `Grid` is used to determine what ranges of the target arrays (i.e. \c
  Vecs into which NetCDF information will be interpolated) are owned by each
  processor.
*/
LocalInterpCtx::LocalInterpCtx(const grid::InputGridInfo &input, const Grid &grid,
                               const std::vector<double> &z_output,
                               InterpolationType type) {

  grid.ctx()->log()->message(4, "\nRegridding file grid info:\n");
  input.report(*grid.ctx()->log(), 4, grid.ctx()->unit_system());

  // limits of the processor's part of the target computational domain
  const double
    x_min_proc = grid.x(grid.xs()),
    x_max_proc = grid.x(grid.xs() + grid.xm() - 1),
    y_min_proc = grid.y(grid.ys()),
    y_max_proc = grid.y(grid.ys() + grid.ym() - 1);

  // T
  start[T] = input.t_len - 1;       // use the latest time.
  count[T] = 1;                     // read only one record

  // X
  subset_start_and_count(input.x, x_min_proc, x_max_proc, start[X], count[X]);

  // Y
  subset_start_and_count(input.y, y_min_proc, y_max_proc, start[Y], count[Y]);

  // Z
  start[Z] = 0;                    // always start at the base
  count[Z] = std::max((int)input.z.size(), 1); // read at least one level

  if (type == LINEAR or type == NEAREST) {
    x.reset(new Interpolation(type, &input.x[start[X]], count[X], &grid.x()[grid.xs()], grid.xm()));

    y.reset(new Interpolation(type, &input.y[start[Y]], count[Y], &grid.y()[grid.ys()], grid.ym()));

    z.reset(new Interpolation(type, input.z, z_output));
  } else {
    throw RuntimeError(PISM_ERROR_LOCATION, "invalid interpolation type in LocalInterpCtx");
  }
}

int LocalInterpCtx::buffer_size() const {
  const int T = 0, X = 1, Y = 2, Z = 3; // indices, just for clarity
  return count[X] * count[Y] * std::max(count[Z], 1);
}


} // end of namespace pism

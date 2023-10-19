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
#include <memory>
#include <vector>

#include "IO_Flags.hh"
#include "pism/util/io/LocalInterpCtx.hh"
#include "pism/util/Grid.hh"

#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/interpolation.hh"

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
static void subset_start_and_count(const std::vector<double> &x, double subset_x_min,
                                   double subset_x_max, int &x_start, int &x_count) {
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
LocalInterpCtx::LocalInterpCtx(const grid::InputGridInfo &input_grid, const Grid &internal_grid,
                               const std::vector<double> &z_internal, InterpolationType type)
  : LocalInterpCtx(input_grid, internal_grid, type) {

  // Z
  start[Z_AXIS] = 0;                                     // always start at the base
  count[Z_AXIS] = std::max((int)input_grid.z.size(), 1); // read at least one level

  if (type == LINEAR or type == NEAREST) {
    z.reset(new Interpolation(type, input_grid.z, z_internal));
  } else {
    throw RuntimeError(PISM_ERROR_LOCATION, "invalid interpolation type in LocalInterpCtx");
  }
}

/*!
 * The two-dimensional version of the interpolation context.
 */
LocalInterpCtx::LocalInterpCtx(const grid::InputGridInfo &input_grid, const Grid &internal_grid,
                               InterpolationType type) {

  // limits of the processor's part of the target computational domain
  const double x_min_proc = internal_grid.x(internal_grid.xs()),
               x_max_proc = internal_grid.x(internal_grid.xs() + internal_grid.xm() - 1),
               y_min_proc = internal_grid.y(internal_grid.ys()),
               y_max_proc = internal_grid.y(internal_grid.ys() + internal_grid.ym() - 1);

  // T
  start[T_AXIS] = (int)input_grid.t_len - 1; // use the latest time.
  count[T_AXIS] = 1;                         // read only one record

  // X
  subset_start_and_count(input_grid.x, x_min_proc, x_max_proc, start[X_AXIS], count[X_AXIS]);

  // Y
  subset_start_and_count(input_grid.y, y_min_proc, y_max_proc, start[Y_AXIS], count[Y_AXIS]);

  // Z
  start[Z_AXIS] = 0;
  count[Z_AXIS] = 1;

  if (type == LINEAR or type == NEAREST) {
    x = std::make_shared<Interpolation>(type, &input_grid.x[start[X_AXIS]], count[X_AXIS],
                                        &internal_grid.x()[internal_grid.xs()], internal_grid.xm());

    y = std::make_shared<Interpolation>(type, &input_grid.y[start[Y_AXIS]], count[Y_AXIS],
                                        &internal_grid.y()[internal_grid.ys()], internal_grid.ym());

    std::vector<double> zz = {0.0};
    z = std::make_shared<Interpolation>(type, zz, zz);
  } else {
    throw RuntimeError(PISM_ERROR_LOCATION, "invalid interpolation type in LocalInterpCtx");
  }
}


int LocalInterpCtx::buffer_size() const {
  return count[X_AXIS] * count[Y_AXIS] * std::max(count[Z_AXIS], 1);
}


} // end of namespace pism

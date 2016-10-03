// Copyright (C) 2007-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cstdlib>
#include <algorithm>            // std::min
#include <gsl/gsl_interp.h>

#include "PIO.hh"
#include "base/util/pism_const.hh"
#include "LocalInterpCtx.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/IceGrid.hh"

#include "base/util/interpolation.hh"
#include "base/util/error_handling.hh"
#include "base/util/Logger.hh"

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
                                   unsigned int &x_start, unsigned int &x_count) {
  unsigned int x_size = x.size();

  x_start = gsl_interp_bsearch(&x[0], subset_x_min, 0, x_size - 1);

  unsigned int x_end = gsl_interp_bsearch(&x[0], subset_x_max, 0, x_size - 1) + 1;

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

  The `IceGrid` is used to determine what ranges of the target arrays (i.e. \c
  Vecs into which NetCDF information will be interpolated) are owned by each
  processor.
*/
LocalInterpCtx::LocalInterpCtx(const grid_info &input, const IceGrid &grid,
                               double z_min, double z_max) {
  const int T = 0, X = 1, Y = 2, Z = 3; // indices, just for clarity

  com = grid.com;
  rank = grid.rank();
  report_range = true;

  grid.ctx()->log()->message(3, "\nRegridding file grid info:\n");
  input.report(*grid.ctx()->log(), 3, grid.ctx()->unit_system());

  // Grid spacing (assume that the grid is equally-spaced) and the
  // extent of the domain. To compute the extent of the domain, assume
  // that the grid is cell-centered, so edge of the domain is half of
  // the grid spacing away from grid points at the edge.
  const double
    x_min = grid.x0() - grid.Lx(),
    x_max = grid.x0() + grid.Lx(),
    y_min = grid.y0() - grid.Ly(),
    y_max = grid.y0() + grid.Ly();
  const double
    input_x_min = input.x0 - input.Lx,
    input_x_max = input.x0 + input.Lx,
    input_y_min = input.y0 - input.Ly,
    input_y_max = input.y0 + input.Ly;

  // tolerance (one micron)
  double eps = 1e-6;
  if (not (x_min >= input_x_min - eps and x_max <= input_x_max + eps and
           y_min >= input_y_min - eps and y_max <= input_y_max + eps and
           z_min >= input.z_min - eps and z_max <= input.z_max + eps)) {

    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "target computational domain not a subset of source (in NetCDF file)\n"
                                  "computational domain:\n"
                                  "target domain: [%3.3f, %3.3f] x [%3.3f, %3.3f] x [%3.3f, %3.3f] meters\n"
                                  "source domain: [%3.3f, %3.3f] x [%3.3f, %3.3f] x [%3.3f, %3.3f] meters",
                                  x_min, x_max, y_min, y_max, z_min, z_max,
                                  input_x_min, input_x_max, input_y_min, input_y_max, input.z_min, input.z_max);
  }

  // limits of the processor's part of the target computational domain
  double
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
  count[Z] = std::max(input.z_len, 1u); // read at least one level
  zlevels = input.z;

  // We need a buffer for the local data, but node 0 needs to have as much
  // storage as the node with the largest block (which may be anywhere), hence
  // we perform a reduce so that node 0 has the maximum value.
  unsigned int buffer_size = count[X] * count[Y] * std::max(count[Z], 1u);
  unsigned int proc0_buffer_size = buffer_size;
  MPI_Reduce(&buffer_size, &proc0_buffer_size, 1, MPI_UNSIGNED, MPI_MAX, 0, com);

  ParallelSection allocation(grid.com);
  try {
    if (rank == 0) {
      buffer.resize(proc0_buffer_size);
    } else {
      buffer.resize(buffer_size);
    }
  } catch (...) {
    allocation.failed();
  }
  allocation.check();

  // Compute indices of neighbors and map-plane interpolation coefficients.
  LinearInterpolation x_interp(&input.x[start[X]], count[X],
                               &grid.x()[grid.xs()], grid.xm());
  x_left = x_interp.left();
  x_right = x_interp.right();
  x_alpha = x_interp.alpha();

  LinearInterpolation y_interp(&input.y[start[Y]], count[Y],
                               &grid.y()[grid.ys()], grid.ym());
  y_left = y_interp.left();
  y_right = y_interp.right();
  y_alpha = y_interp.alpha();
}

LocalInterpCtx::~LocalInterpCtx() {
  // empty
}

} // end of namespace pism

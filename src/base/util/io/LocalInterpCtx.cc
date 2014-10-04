// Copyright (C) 2007-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "PIO.hh"
#include "pism_const.hh"
#include "LocalInterpCtx.hh"

namespace pism {

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
  rank = grid.rank;
  report_range = true;

  print_grid_info(input, grid.get_unit_system(), 3);

  // Grid spacing (assume that the grid is equally-spaced) and the
  // extent of the domain. To compute the extent of the domain, assume
  // that the grid is cell-centered, so edge of the domain is half of
  // the grid spacing away from grid points at the edge.
  const double
    x_min = grid.x0 - grid.Lx,
    x_max = grid.x0 + grid.Lx,
    y_min = grid.y0 - grid.Ly,
    y_max = grid.y0 + grid.Ly;
  const double
    input_x_min = input.x0 - input.Lx,
    input_x_max = input.x0 + input.Lx,
    input_y_min = input.y0 - input.Ly,
    input_y_max = input.y0 + input.Ly;

  // tolerance (one micron)
  double eps = 1e-6;
  if (not (x_min >= input_x_min - eps && x_max <= input_x_max + eps &&
           y_min >= input_y_min - eps && y_max <= input_y_max + eps &&
           z_min >= input.z_min - eps && z_max <= input.z_max + eps)) {

    PetscPrintf(com,
                "target computational domain not a subset of source (in NetCDF file)\n"
                "  computational domain:\n");
    PetscPrintf(grid.com, "target domain: [%3.3f, %3.3f] x [%3.3f, %3.3f] x [%3.3f, %3.3f] meters\n",
                x_min, x_max,
                y_min, y_max,
                z_min, z_max);
    PetscPrintf(grid.com, "source domain: [%3.3f, %3.3f] x [%3.3f, %3.3f] x [%3.3f, %3.3f] meters\n",
                input_x_min, input_x_max,
                input_y_min, input_y_max,
                input.z_min, input.z_max);
    PISMEnd();
  }

  // limits of the processor's part of the target computational domain
  double
    x_min_proc = grid.x[grid.xs],
    x_max_proc = grid.x[grid.xs + grid.xm - 1],
    y_min_proc = grid.y[grid.ys],
    y_max_proc = grid.y[grid.ys + grid.ym - 1];

  // T
  start[T] = input.t_len - 1;       // use the latest time.
  count[T] = 1;                     // read only one record

  // X
  start[X] = 0;
  while (start[X] + 1 < input.x_len && input.x[start[X] + 1] < x_min_proc)
    start[X]++;

  count[X] = 1;
  while (start[X] + count[X] < input.x_len && input.x[start[X] + count[X] - 1] <= x_max_proc)
    count[X]++;

  // Y
  start[Y] = 0;
  while (start[Y] + 1 < input.y_len && input.y[start[Y] + 1] < y_min_proc)
    start[Y]++;

  count[Y] = 1;
  while (start[Y] + count[Y] < input.y_len && input.y[start[Y] + count[Y] - 1] <= y_max_proc)
    count[Y]++;

  // Z
  start[Z] = 0;                    // always start at the base
  count[Z] = PetscMax(input.z_len, 1); // read at least one level
  zlevels = input.z;

  // We need a buffer for the local data, but node 0 needs to have as much
  // storage as the node with the largest block (which may be anywhere), hence
  // we perform a reduce so that node 0 has the maximum value.
  a_len = count[X] * count[Y] * PetscMax(count[Z], 1);
  int my_a_len = a_len;
  MPI_Reduce(&my_a_len, &(a_len), 1, MPI_INT, MPI_MAX, 0, com);
  PetscMalloc(a_len * sizeof(double), &(a));
  // FIXME: we need error checking here

  // Compute indices of neighbors and map-plane interpolation coefficients.
  x_left.resize(grid.xm);
  x_right.resize(grid.xm);
  x_alpha.resize(grid.xm);

  // x-direction; loop over all grid points in this processor's sub-domain
  for (int i = 0; i < grid.xm; ++i) {
    // i is the index in this processor's sub-domain
    // x is the x-coordinate in the total domain
    double x = grid.x[grid.xs + i];

    // This is here to make it crash and burn if something goes wrong, instead
    // of quietly doing the wrong thing.
    x_left[i]  = -1;
    x_right[i] = -1;
    x_alpha[i] = -1;

    // go through all the grid intervals in the x-direction in this
    // processor's sub-domain:
    for (unsigned int k = 0; k < count[X] - 1; ++k) {
      unsigned int kk = k + start[X]; // global grid point index
      if (input.x[kk] <= x && input.x[kk + 1] >= x) {
        // if the i-th point is in the current interval
        x_left[i]  = std::min(k,     count[X] - 1);
        x_right[i] = std::min(k + 1, count[X] - 1);
        x_alpha[i] = (x - input.x[kk]) / (input.x[kk + 1] - input.x[kk]);
        break;
      }
    }
  }

  // Check all values of x_left to make sure that the edge of the
  // computational domain is handled correctly. In PISM each grid
  // point is at the center of a grid cell, so it may be impossible to
  // use *interpolation* to transfer a field from a coarse grid to a
  // finer grid: some points of the fine grid may be *inside* the
  // domain covered by the coarse grid, but *outside* the convex hull
  // of all the coarse grid points. So, we modify x_left, x_right, and
  // x_alpha to use constant extrapolation for these fine grid points.

  for (int i = 0; i < grid.xm; ++i) {
    // for all points in the x-direction in this sub-domain

    // get the coordinate
    double x = grid.x[grid.xs + i];
    if (x_left[i] == -1 && x < input.x[0]) {
      // if x_left was not assigned and the point is to the left of
      // the left-most column in the input grid
      x_left[i]  = 0;
      x_right[i] = 0;
      x_alpha[i] = 0;
    }

    if (x_left[i] == -1 && x > input.x.back()) {
      x_left[i]  = count[X] - 1;
      x_right[i] = count[X] - 1;
      x_alpha[i] = 0;
    }
  }

  // y-direction
  y_left.resize(grid.ym);
  y_right.resize(grid.ym);
  y_alpha.resize(grid.ym);

  for (int j = 0; j < grid.ym; ++j) {
    double y = grid.y[grid.ys + j];

    // This is here to make it crash and burn if something goes wrong, instead
    // of quietly doing the wrong thing.
    y_left[j]  = -1;
    y_right[j] = -1;
    y_alpha[j] = -1;

    for (unsigned int k = 0; k < count[Y] - 1; ++k) {
      unsigned int kk = k + start[Y];
      if (input.y[kk] <= y && input.y[kk + 1] >= y) {
        y_left[j]  = std::min(k,     count[Y] - 1);
        y_right[j] = std::min(k + 1, count[Y] - 1);
        y_alpha[j] = (y - input.y[kk]) / (input.y[kk + 1] - input.y[kk]);
        break;
      }
    }
  }

  // Take care of the edge if the domain in the y-direction (see above
  // for details).
  for (int i = 0; i < grid.ym; ++i) {
    // for all points in the y-direction in this sub-domain

    // get the coordinate
    double y = grid.y[grid.ys + i];
    if (y_left[i] == -1 && y < input.y[0]) {
      // if y_left was not assigned and the point is below the
      // bottom row in the input grid
      y_left[i]  = 0;
      y_right[i] = 0;
      y_alpha[i] = 0;
    }

    if (y_left[i] == -1 && y > input.y.back()) {
      y_left[i]  = count[Y] - 1;
      y_right[i] = count[Y] - 1;
      y_alpha[i] = 0;
    }
  }
}


//! Deallocate memory.
LocalInterpCtx::~LocalInterpCtx() {
  PetscFreeVoid(a);
}

//! Print out the actual array information stored in the local interpolation context.
/*!
Every processor in the communicator `com` must call this for it to work, I think.
 */
PetscErrorCode LocalInterpCtx::printArray() {
  PetscErrorCode ierr;

  ierr = PetscSynchronizedPrintf(com,"\nLocalInterpCtx::printArray():  rank = %d, a_len = %d\n",
             rank, a_len); CHKERRQ(ierr);
  for (unsigned int k = 0; k < a_len; k++) {
    ierr = PetscSynchronizedPrintf(com," %5.4f,",a[k]); CHKERRQ(ierr);
  }
#if PETSC_VERSION_LT(3,5,0)
  ierr = PetscSynchronizedFlush(com); CHKERRQ(ierr);
#else
  ierr = PetscSynchronizedFlush(com, NULL); CHKERRQ(ierr);
#endif
  return 0;
}

void LocalInterpCtx::print_grid_info(const grid_info &g, const UnitSystem &s, int threshold) {

  verbPrintf(threshold, com,
             "\nRegridding file grid info:\n");

  verbPrintf(threshold, com,
             "  x:  %5d points, [%10.3f, %10.3f] km, x0 = %10.3f km, Lx = %10.3f km\n",
             g.x_len,
             (g.x0 - g.Lx)/1000.0,
             (g.x0 + g.Lx)/1000.0,
             g.x0/1000.0,
             g.Lx/1000.0);

  verbPrintf(threshold, com,
             "  y:  %5d points, [%10.3f, %10.3f] km, y0 = %10.3f km, Ly = %10.3f km\n",
             g.y_len,
             (g.y0 - g.Ly)/1000.0,
             (g.y0 + g.Ly)/1000.0,
             g.y0/1000.0,
             g.Ly/1000.0);

  verbPrintf(threshold, com,
             "  z:  %5d points, [%10.3f, %10.3f] m\n",
             g.z_len, g.z_min, g.z_max);

  verbPrintf(threshold, com,
             "  t:  %5d points, last time = %.3f years\n\n",
             g.t_len, s.convert(g.time, "seconds", "years"));
}

} // end of namespace pism

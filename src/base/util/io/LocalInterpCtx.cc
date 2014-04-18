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
#include "PIO.hh"
#include "pism_const.hh"
#include "LocalInterpCtx.hh"

namespace pism {

//! Construct a local interpolation context from arrays of parameters.
/*!
  This method constructs a class from existing information already read from a NetCDF file and stored
  in arrays.

  The essential quantities to compute are where each processor should start within the NetCDF file grid
  (`start[]`) and how many grid points, from the starting place, the processor has.  The resulting
  portion of the grid is stored in array `a` (a field of the `LocalInterCtx`).

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
LocalInterpCtx::LocalInterpCtx(grid_info input, const IceGrid &grid,
                               double z_min, double z_max) {
  const int T = 0, X = 1, Y = 2, Z = 3; // indices, just for clarity

  com = grid.com;
  rank = grid.rank;
  report_range = true;

  print_grid_info(input, grid.get_unit_system(), 3);

  double eps = 1e-6;         // tolerance (one micron)
  if (!(grid.x.front() >= input.x_min - eps && grid.x.back() <= input.x_max + eps &&
        grid.y.front() >= input.y_min - eps && grid.y.back() <= input.y_max + eps &&
        z_min >= input.z_min - eps && z_max <= input.z_max + eps)) {

    PetscPrintf(com,
                "target computational domain not a subset of source (in NetCDF file)\n"
                "  computational domain:\n");
    PetscPrintf(grid.com, "target domain: [%3.3f, %3.3f] x [%3.3f, %3.3f] x [%3.3f, %3.3f] meters\n",
                grid.x.front(), grid.x.back(),
                grid.y.front(), grid.y.back(),
                z_min, z_max);
    PetscPrintf(grid.com, "source domain: [%3.3f, %3.3f] x [%3.3f, %3.3f] x [%3.3f, %3.3f] meters\n",
                input.x_min, input.x_max,
                input.y_min, input.y_max,
                input.z_min, input.z_max);
    PISMEnd();
  }

  // limits of the processor's part of the target computational domain
  double x_min = grid.x[grid.xs],
    x_max = grid.x[grid.xs + grid.xm - 1],
    y_min = grid.y[grid.ys],
    y_max = grid.y[grid.ys + grid.ym - 1];

  // T
  start[T] = input.t_len - 1;       // use the latest time.
  count[T] = 1;                     // read only one record

  // X
  start[X] = 0;
  while (start[X] + 1 < input.x_len && input.x[start[X] + 1] < x_min)
    start[X]++;

  count[X] = 1;
  while (start[X] + count[X] < input.x_len && input.x[start[X] + count[X] - 1] <= x_max)
    count[X]++;

  // Y
  start[Y] = 0;
  while (start[Y] + 1 < input.y_len && input.y[start[Y] + 1] < y_min)
    start[Y]++;

  count[Y] = 1;
  while (start[Y] + count[Y] < input.y_len && input.y[start[Y] + count[Y] - 1] <= y_max)
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

  // Compute indices of neighbors and map-plane interpolation coefficients.
  x_left.resize(grid.xm);
  x_right.resize(grid.xm);
  x_alpha.resize(grid.xm);

  // x-direction
  for (int i = 0; i < grid.xm; ++i) {
    double x = grid.x[grid.xs + i];

    // This is here to make it crash and burn if something goes wrong, instead
    // of quietly doing the wrong thing.
    x_left[i]  = -1;
    x_right[i] = -1;
    x_alpha[i] = -1;

    for (unsigned int k = 0; k < count[X] - 1; ++k) {
      unsigned kk = k + start[X];
      if (input.x[kk] <= x && input.x[kk + 1] >= x) {
        x_left[i]  = PetscMin(k,     count[X] - 1);
        x_right[i] = PetscMin(k + 1, count[X] - 1);
        x_alpha[i] = (x - input.x[kk]) / (input.x[kk + 1] - input.x[kk]);
        break;
      }
    }
  }

  // Check the smallest and the biggest x (all other indices are guaranteed to
  // be initialized; these may remain "-1" because of rounding errors).

  // (Note that would not be necessary if we didn't have "eps" at the top of
  // this constructor, but that would make our code less robust. Also note that
  // the code below does not cover the case of dx < eps with eps defined above.
  // Something is seriously wrong if dx < "one micron", though.)
  if (x_left[0] == -1)
    x_alpha[0] = x_left[0] = x_right[0] = 0;

  if (x_left[grid.xm - 1] == -1) {
    int k = grid.xm - 1;
    x_left[k]  = count[X] - 1;
    x_right[k] = count[X] - 1;
    x_alpha[k] = 0;
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
        y_left[j]  = PetscMin(k,     count[Y] - 1);
        y_right[j] = PetscMin(k + 1, count[Y] - 1);
        y_alpha[j] = (y - input.y[kk]) / (input.y[kk + 1] - input.y[kk]);
        break;
      }
    }
  }

  // Check the smallest and the biggest y (all other indices are guaranteed to
  // be initialized; these may remain "-1" because of rounding errors).
  if (y_left[0] == -1)
    y_alpha[0] = y_left[0]  = y_right[0] = 0;

  if (y_left[grid.ym - 1] == -1) {
    int k = grid.ym - 1;
    y_left[k]  = count[Y] - 1;
    y_right[k] = count[Y] - 1;
    y_alpha[k] = 0;
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
  for (int k = 0; k < a_len; k++) {
    ierr = PetscSynchronizedPrintf(com," %5.4f,",a[k]); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush(com); CHKERRQ(ierr);
  return 0;
}

void LocalInterpCtx::print_grid_info(grid_info g, PISMUnitSystem s, int threshold) {

  verbPrintf(threshold, com,
             "\nRegridding file grid info:\n");

  verbPrintf(threshold, com,
             "  x:  %5d points, [%10.3f, %10.3f] km, x0 = %10.3f km, Lx = %10.3f km\n",
             g.x_len, g.x_min/1000.0, g.x_max/1000.0, (g.x_min + g.x_max)/2.0/1000.0,
             (g.x_max - g.x_min)/1000.0);

  verbPrintf(threshold, com,
             "  y:  %5d points, [%10.3f, %10.3f] km, y0 = %10.3f km, Ly = %10.3f km\n",
             g.y_len, g.y_min/1000.0, g.y_max/1000.0,
             (g.y_min + g.y_max)/2.0/1000.0,
             (g.y_max - g.y_min)/1000.0);

  verbPrintf(threshold, com,
             "  z:  %5d points, [%10.3f, %10.3f] m\n",
             g.z_len, g.z_min, g.z_max);

  verbPrintf(threshold, com,
             "  t:  %5d points, last time = %.3f years\n\n",
             g.t_len, s.convert(g.time, "seconds", "years"));
}

} // end of namespace pism

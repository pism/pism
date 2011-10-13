// Copyright (C) 2007-2011 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
#include "NCTool.hh"
#include "pism_const.hh"
#include "LocalInterpCtx.hh"


//! Construct a local interpolation context from arrays of parameters.
/*!
  This method constructs a class from existing information already read from a NetCDF file and stored
  in arrays.  
  
  The essential quantities to compute are where each processor should start within the NetCDF file grid
  (<tt>start[]</tt>) and how many grid points, from the starting place, the processor has.  The resulting
  portion of the grid is stored in array \c a (a field of the \c LocalInterCtx).

  We make conservative choices about \c start[] and \c count[].  In particular, the portions owned by 
  processors \e must overlap at one point in the NetCDF file grid, but they \e may overlap more than that
  (as computed here).

  Note this constructor doesn't extract new information from the NetCDF file or
  do communication. The information from the NetCDF file must already be
  extracted, validly stored in a grid_info structure \c input.

  The \c IceGrid is used to determine what ranges of the target arrays (i.e. \c
  Vecs into which NetCDF information will be interpolated) are owned by each
  processor.
*/
LocalInterpCtx::LocalInterpCtx(grid_info input, IceGrid &grid,
                               PetscReal z_min, PetscReal z_max) {
  const int T = 0, X = 1, Y = 2, Z = 3; // indices, just for clarity

  com = grid.com;
  rank = grid.rank;
  report_range = true;

  input.print(com);

  if (!(grid.x.front() >= input.x_min && grid.x.back() <= input.x_max &&
        grid.y.front() >= input.y_min && grid.y.back() <= input.y_max &&
        z_min >= input.z_min && z_max <= input.z_max)) {

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
  for (PetscInt i = 0; i < grid.xm; ++i) {
    double x = grid.x[grid.xs + i];

    for (unsigned int k = start[X]; k < start[X] + count[X]; ++k) {
      if (input.x[k] <= x && input.x[k + 1] >= x) {
        x_left[i]  = PetscMin(k     - start[X], count[X] - 1);
        x_right[i] = PetscMin(k + 1 - start[X], count[X] - 1);
        x_alpha[i] = (x - input.x[k]) / (input.x[k + 1] - input.x[k]);
        break;
      }
    }
  }

  // y-direction
  y_left.resize(grid.ym);
  y_right.resize(grid.ym);
  y_alpha.resize(grid.ym);

  for (PetscInt j = 0; j < grid.ym; ++j) {
    double y = grid.y[grid.ys + j];

    for (unsigned int k = start[Y]; k < start[Y] + count[Y]; ++k) {
      if (input.y[k] <= y && input.y[k + 1] >= y) {
        y_left[j]  = PetscMin(k     - start[Y], count[Y] - 1);
        y_right[j] = PetscMin(k + 1 - start[Y], count[Y] - 1);
        y_alpha[j] = (y - input.y[k]) / (input.y[k + 1] - input.y[k]);
        break;
      }
    }
  }
}


//! Deallocate memory.
LocalInterpCtx::~LocalInterpCtx() {
  PetscFreeVoid(a);
}

//! Print out the actual array information stored in the local interpolation context.
/*!
Every processor in the communicator \c com must call this for it to work, I think.
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

grid_info::grid_info() {
  t_len = 0;
  time  = 0;

  x_len = 0;
  x_max = 0;
  x_min = 0;

  y_len = 0;
  y_max = 0;
  y_min = 0;

  z_len = 0;
  z_min = 0;
  z_max = 0;
}

PetscErrorCode grid_info::print(MPI_Comm com, int threshold) {
  PetscErrorCode ierr;
  double zero = 0;
  ierr = verbPrintf(threshold, com, "\nRegridding file grid info:\n"); CHKERRQ(ierr);

  ierr = verbPrintf(threshold, com, "  x:  %5d points, [%10.3f, %10.3f] km, x0 = %10.3f km, Lx = %10.3f km\n",
		    x_len, x_min/1000.0, x_max/1000.0, (x_min + x_max)/2.0/1000.0,
                    (x_max - x_min)/1000.0); CHKERRQ(ierr);
  ierr = verbPrintf(threshold, com, "  y:  %5d points, [%10.3f, %10.3f] km, y0 = %10.3f km, Ly = %10.3f km\n",
		    y_len, y_min/1000.0, y_max/1000.0,
                    (y_min + y_max)/2.0/1000.0,
                    (y_max - y_min)/1000.0); CHKERRQ(ierr);
  ierr = verbPrintf(threshold, com, "  z:  %5d points, [%10.3f, %10.3f] m\n",
		    z_len, zero, z_max); CHKERRQ(ierr);
  ierr = verbPrintf(threshold, com, "  t:  %5d points, last time = %.3f years\n\n",
		    t_len, convert(time, "seconds", "years")); CHKERRQ(ierr);
  return 0;
}


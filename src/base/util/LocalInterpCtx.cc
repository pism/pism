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
  extracted, validly stored in a grid_info structure \c g and preallocated arrays
  \c zlevsIN[], \c zblevsIN[], and each processor must have the same values for
  all of them.

  The \c IceGrid is used to determine what ranges of the target arrays (i.e. \c Vecs into which NetCDF
  information will be interpolated) are owned by each processor.
*/
LocalInterpCtx::LocalInterpCtx(grid_info g, IceGrid &grid) {
  const int T = 0, X = 1, Y = 2, Z = 3; // indices, just for clarity
  const double Lx = grid.Lx,
               Ly = grid.Ly,
               slop = 1.000001; // allowed slop; grids must be subsets within this factor
  com = grid.com;
  rank = grid.rank;
  report_range = true;

  g.print(com);

  PetscScalar dx0 = grid.x0 - g.x0,
              dy0 = grid.y0 - g.y0;
  if (sqrt(dx0*dx0 + dy0*dy0) > 1e-6) {
    PetscPrintf(com, "PISM ERROR: Grid centers do not match!\n");
    PISMEnd();
  }

  if (g.Lx*slop < Lx) {
    PetscPrintf(com,
		"target computational domain not a subset of source (in NetCDF file)\n"
		"  computational domain:\n");
    PetscPrintf(com,
		"    need  Lx < g.Lx,  but Lx = %5.4f km while g.Lx = %5.4f km\n"
		"    ENDING ...\n",
		Lx, g.Lx);
    PISMEnd();
  }
  if (g.Ly*slop < Ly) {
    PetscPrintf(com,
		"target computational domain not a subset of source (in NetCDF file)\n"
		"  computational domain:\n");
    PetscPrintf(com,
		"    need  Ly < g.Ly,  but Ly = %5.4f km while g.Ly = %5.4f km\n"
		"    ENDING ...\n",
		Ly, g.Ly);
   PISMEnd();
  }

  // FIXME: put back vertical extent check

  // limits of the processor's part of the target computational domain
  double xbdy_tgt[2] = {grid.x[grid.xs] - grid.x0,
                        grid.x[grid.xs + grid.xm - 1] - grid.x0};
  double ybdy_tgt[2] = {grid.y[grid.ys] - grid.y0,
                        grid.y[grid.ys + grid.ym - 1] - grid.y0};

/*
To make this work with unequal spacing <i>in the horizontal dimension</i>, IF EVER NEEDED,
we have some choices.
Suppose a processor owns indices \f$\{i_m, \dots, i_{m'}\}\f$ and we know nothing about the
spacing.  Then to find the interpolated value at \f$x(j)\f$ where \f$i_m \le j \le i_{m'}\f$
we need the index \f$J\f$ such that \f$X(J) \le x(j) \le X(J+1)\f$.

Note that we could just loop through an array of \f$X(\cdot)\f$ to find \f$J\f$, and it would
not be a performance bottleneck.  It would also be more general.  Of course, for
a special case like Chebyshev-Gauss-Lobatto we can compute the indices.  A good
approach would be to have a structure representing the layout in each dimension.
Then we can have a function which takes a floating point value and returns the
largest index which is not greater than that value.

Note that \c lic.start and \c lic.count is all that is necessary to pull the correct
data from the netCDF file, so if we implement this general scheme, the \c fstart
and \c delta entries in the struct will not be meaningful.
 */
 
  // Distances between entries (i.e. dx and dy and dz) in the netCDF file (floating point).
  delta[T] = GSL_NAN; // Delta probably will never make sense in the time dimension.
  delta[X] = (g.x_max - g.x_min) / (g.x_len - 1);
  delta[Y] = (g.y_max - g.y_min) / (g.y_len - 1);

  // start[i] = index of the first needed entry in the source netCDF file; start[i] is of type int
  // We use the latest time.
  start[T] = g.t_len - 1;

  // The following intentionally under-counts the number of times delta[X] fits in the
  // interval (xbdy_tgt[0], -g.Lx), which is the number of points to skip in
  // the X-direction, and therefore the index of the first needed point
  // (because indices start at zero).
  start[X] = (int)floor((xbdy_tgt[0] - (-g.Lx)) / delta[X] - 0.5);

  // Same in the Y-direction.
  start[Y] = (int)floor((ybdy_tgt[0] - (-g.Ly)) / delta[Y] - 0.5);

  // be sure the start[X] and start[Y] are not too small:
  if (start[X] < 0) start[X] = 0;
  if (start[Y] < 0) start[Y] = 0;

  start[Z] = 0;			// start at base of ice

  fstart[T] = g.time;
  fstart[X] = -g.Lx + start[X] * delta[X];
  fstart[Y] = -g.Ly + start[Y] * delta[Y];

  count[T] = 1;			// Only take one time.
  count[X] = 1 + (int)ceil((xbdy_tgt[1] - fstart[X]) / delta[X]);
  count[Y] = 1 + (int)ceil((ybdy_tgt[1] - fstart[Y]) / delta[Y]);
  // make count[] smaller if start[]+count[] would go past end of range:
  while (start[X]+count[X] > g.x_len)
    count[X]--;
  while (start[Y]+count[Y] > g.y_len)
    count[Y]--;

  count[Z] = PetscMax(g.z_len, 1); // read at least one level

  zlevels = g.zlevels;

  // We need a buffer for the local data, but node 0 needs to have as much
  // storage as the node with the largest block (which may be anywhere), hence
  // we perform a reduce so that node 0 has the maximum value.
  a_len = count[X] * count[Y] * PetscMax(count[Z], 1);
  int my_a_len = a_len;
  MPI_Reduce(&my_a_len, &(a_len), 1, MPI_INT, MPI_MAX, 0, com);
  PetscMalloc(a_len * sizeof(double), &(a));
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
  t_len  = 0;
  x_len  = 0;
  y_len  = 0;
  z_len  = 0;
  zb_len = 0;
  x0     = 0;
  y0     = 0;
  Lx     = 0;
  Ly     = 0;
  x_min  = 0;
  x_max  = 0;
  y_min  = 0;
  y_max  = 0;
  z_max  = 0;
  zb_min = 0;
}

PetscErrorCode grid_info::print(MPI_Comm com, int threshold) {
  PetscErrorCode ierr;
  double zero = 0;
  ierr = verbPrintf(threshold, com, "\nRegridding file grid info:\n"); CHKERRQ(ierr);

  ierr = verbPrintf(threshold, com, "  x:  %5d points, [%10.3f, %10.3f] km, x0 = %10.3f km, Lx = %10.3f km\n",
		    x_len, x_min/1000.0, x_max/1000.0, x0/1000.0, Lx/1000.0); CHKERRQ(ierr);
  ierr = verbPrintf(threshold, com, "  y:  %5d points, [%10.3f, %10.3f] km, y0 = %10.3f km, Ly = %10.3f km\n",
		    y_len, y_min/1000.0, y_max/1000.0, y0/1000.0, Ly/1000.0); CHKERRQ(ierr);
  ierr = verbPrintf(threshold, com, "  z:  %5d points, [%10.3f, %10.3f] m\n",
		    z_len, zero, z_max); CHKERRQ(ierr);
  ierr = verbPrintf(threshold, com, "  zb: %5d points, [%10.3f, %10.3f] m\n",
		    zb_len, zb_min, zero); CHKERRQ(ierr);
  ierr = verbPrintf(threshold, com, "  t:  %5d points, last time = %.3f years\n\n",
		    t_len, convert(time, "seconds", "years")); CHKERRQ(ierr);
  return 0;
}


// Copyright (C) 2007-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <netcdf.h>
#include "nc_util.hh"
#include "pism_const.hh"
#include "iceModel.hh"
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
LocalInterpCtx::LocalInterpCtx(grid_info g,
			       const double zlevsIN[], const double zblevsIN[], IceGrid &grid) {
  const int T = 0, X = 1, Y = 2, Z = 3, ZB = 4; // indices, just for clarity
  const double Lx = grid.Lx,
               Ly = grid.Ly,
               Lz = grid.Lz,
               Lbz = grid.Lbz,
               dx = grid.dx,
               dy = grid.dy,
               slop = 1.000001; // allowed slop; grids must subsets within this factor
  com = grid.com;
  rank = grid.rank;
  regrid_2d_only = false;
  no_regrid_bedrock = false;
  report_range = true;

  g.print(com);

  if ((zlevsIN == NULL) || (zblevsIN == NULL))
    regrid_2d_only = true;

  PetscScalar dx0 = grid.x0 - g.x0,
              dy0 = grid.y0 - g.y0;
  if (sqrt(dx0*dx0 + dy0*dy0) > 1e-6) {
    PetscPrintf(com, "PISM ERROR: Grid centers do not match!\n");
    PetscEnd();
  }

  if (g.Lx*slop < Lx) {
    PetscPrintf(com,
		"target computational domain not a subset of source (in NetCDF file)\n"
		"  computational domain:\n");
    PetscPrintf(com,
		"    need  Lx < g.Lx,  but Lx = %5.4f km while g.Lx = %5.4f km\n"
		"    ENDING ...\n",
		Lx, g.Lx);
    PetscEnd();
  }
  if (g.Ly*slop < Ly) {
    PetscPrintf(com,
		"target computational domain not a subset of source (in NetCDF file)\n"
		"  computational domain:\n");
    PetscPrintf(com,
		"    need  Ly < g.Ly,  but Ly = %5.4f km while g.Ly = %5.4f km\n"
		"    ENDING ...\n",
		Ly, g.Ly);
   PetscEnd();
  }
  
  if (regrid_2d_only == false) {
    if (g.z_max*slop < Lz) {
      verbPrintf(3,com,
		 "  WARNING: vertical dimension of target computational domain\n"
		 "    not a subset of source (in NetCDF file) computational domain;\n"
		 "    g.z_max = %5.4f < Lz = %5.4f; ALLOWING ONLY 2D REGRIDDING ...\n",
		 g.z_max, Lz);
      regrid_2d_only = true;
    }
    if (-g.zb_min*slop < Lbz) {
      verbPrintf(3,com,
		 "  LIC WARNING: vertical dimension of target BEDROCK computational domain\n"
		 "    not a subset of source (in NetCDF file) BEDROCK computational domain;\n"
		 "    -g.zb_min = %5.4f < Lbz = %5.4f; NOT ALLOWING BEDROCK REGRIDDING ...\n",
		 -g.zb_min, Lbz);
      no_regrid_bedrock = true;
    } else {
      no_regrid_bedrock = false;
    }

    // This disables regridding bedrock temperature if an input file has only
    // one bedrock layer.
    if (g.zb_len < 2) no_regrid_bedrock = true;
  }

  verbPrintf(5, com, "LIC INFO: regrid_2d_only = %d, no_regrid_bedrock = %d\n",
	     regrid_2d_only, no_regrid_bedrock);

  // limits of the processor's part of the target computational domain
  double xbdy_tgt[2] = {-Lx + dx * grid.xs, -Lx + dx * (grid.xs + grid.xm - 1)};
  double ybdy_tgt[2] = {-Ly + dy * grid.ys, -Ly + dy * (grid.ys + grid.ym - 1)};

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
  start[ZB] = 0;		// start at lowest bedrock level

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

  // This allows creating a local interpolation context for 2D regridding
  // without specifying dummy z or zb-related information in the g argument.
  if (regrid_2d_only) {
    count[Z] = 1;
    count[ZB] = 1;
  } else {
    count[Z] = g.z_len;
    count[ZB] = g.zb_len;
  }

  zlevs = new double[count[Z]];
  zblevs = new double[count[ZB]];

  if (regrid_2d_only) {
    zlevs[0] = 0.0;
    zblevs[0] = 0.0;
  } else {
    for (int k = 0; k < count[Z]; k++) {
      zlevs[k] = zlevsIN[k];
    }
    for (int k = 0; k < count[ZB]; k++) {
      zblevs[k] = zblevsIN[k];
    }
  }

  nz = count[Z]; nzb = count[ZB];

  // We need a buffer for the local data, but node 0 needs to have as much
  // storage as the node with the largest block (which may be anywhere), hence
  // we perform a reduce so that node 0 has the maximum value.
  a_len = count[X] * count[Y] * PetscMax(count[Z], count[ZB]);
  int my_a_len = a_len;
  MPI_Reduce(&my_a_len, &(a_len), 1, MPI_INT, MPI_MAX, 0, com);
  PetscMalloc(a_len * sizeof(double), &(a));
}


//! Deallocate memory.
LocalInterpCtx::~LocalInterpCtx() {
  delete[] zlevs;
  delete[] zblevs;
  PetscFreeVoid(a);
}


//! Print out the grid information stored in the local interpolation context.
/*!
Every processor in the communicator \c com must call this for it to work, I think.
 */
PetscErrorCode LocalInterpCtx::printGrid() {
  PetscErrorCode ierr;
  const int T = 0, X = 1, Y = 2, Z = 3, ZB = 4; // indices, just for clarity

  ierr = PetscSynchronizedPrintf(com,"\nLocalInterpCtx::printGrid():  rank = %d\n",
                    rank); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(com,"  delta[1,2] = %5.4f, %5.4f\n",
                    delta[X],delta[Y]); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(com,"  start[0,..,4] = %d, %d, %d, %d, %d\n",
                    start[T],start[X],start[Y],start[Z],start[ZB]); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(com,"  fstart[0,1,2] = %5.4f, %5.4f, %5.4f\n",
                    fstart[T],fstart[X],fstart[Y]); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(com,"  count[0,..,4] = %d, %d, %d, %d, %d\n",
                    count[T],count[X],count[Y],count[Z],count[ZB]); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(com,"  zlevs[] = "); CHKERRQ(ierr);
  for (int k = 0; k < nz; k++) {
    ierr = PetscSynchronizedPrintf(com," %5.4f,",zlevs[k]); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedPrintf(com,"\n  zblevs[] = "); CHKERRQ(ierr);
  for (int k = 0; k < nzb; k++) {
    ierr = PetscSynchronizedPrintf(com," %5.4f,",zblevs[k]); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedPrintf(com,"\n"); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(com); CHKERRQ(ierr);
  return 0;
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


//! Get the index for the source grid below a given vertical elevation within the ice.
/*!
This code duplicates that in IceGrid::kBelowHeight().
 */
int LocalInterpCtx::kBelowHeight(const double height) {
  if (height < 0.0 - 1.0e-6) {
    PetscPrintf(com, 
       "LocalInterpCtx kBelowHeight(): height = %5.4f is below base of ice (height must be nonnegative)\n",
       height);
    PetscEnd();
  }
  const double Lz = zlevs[nz-1];
  if (height > Lz + 1.0e-6) {
    PetscPrintf(com, 
       "LocalInterpCtx kBelowHeight(): height = %5.4f is above top of computational grid Lz = %5.4f\n",
       height,Lz);
    PetscEnd();
  }
  PetscInt mcurr = 0;
  while (zlevs[mcurr+1] < height) {
    mcurr++;
  }
  return mcurr;
}


//! Get the index for the source grid below a given vertical elevation within the bedrock
int LocalInterpCtx::kbBelowHeight(const double elevation) {
  if (elevation < zblevs[0] - 1.0e-6) {
    PetscPrintf(com, 
       "LocalInterpCtx kbBelowHeight(): elevation = %5.4f is below base of bedrock -Lbz = %5.4f\n",
       elevation, zblevs[0]);
    PetscEnd();
  }
  if (elevation > 0.0 + 1.0e-6) {
    PetscPrintf(com, 
       "LocalInterpCtx kbBelowHeight(): elevation = %5.4f is above top of bedrock at z=0\n",
       elevation);
    PetscEnd();
  }
  double myelevation = PetscMin(0.0,elevation);
  PetscInt mcurr = 0;
  while (zblevs[mcurr+1] < myelevation) {
    mcurr++;
  }
  return mcurr;
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
		    t_len, time/secpera); CHKERRQ(ierr);
  return 0;
}


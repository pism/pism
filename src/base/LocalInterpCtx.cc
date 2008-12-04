// Copyright (C) 2007-2008 Jed Brown, Ed Bueler and Constantine Khroulev
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

Note this constructor doesn't extract new information from the NetCDF file or do communication.
The information from the NetCDF file must already be extracted, validly stored in preallocated arrays 
\c dim[], \c bdy[], \c zlevsIN[], \c zblevsIN[], and each processor must have the same values for 
all of these arrays.

The \c IceGrid is used to determine what ranges of the target arrays (i.e. \c Vecs into which NetCDF
information will be interpolated) are owned by each processor.
 */
LocalInterpCtx::LocalInterpCtx(const size_t dim[], const double bdy[],
                               const double zlevsIN[], const double zblevsIN[], IceGrid &grid) {
  PetscErrorCode ierr;
  const double Lx = grid.Lx,
               Ly = grid.Ly,
               Lz = grid.Lz,
               Lbz = grid.Lbz,
               dx = grid.dx,
               dy = grid.dy;

  if (bdy[1] > -Lx || bdy[2] < Lx) {
    PetscPrintf(grid.com,
        "target computational domain not a subset of source (in NetCDF file)\n"
        "  computational domain:\n");
    PetscPrintf(grid.com,
        "    need  [-Lx,Lx] contained in [bdy[1],bdy[2]],  but Lx = %5.4f km while\n"
        "    [bdy[1],bdy[2]]=[%5.4f,%5.4f] km;  ENDING ...\n",
        Lx,bdy[1],bdy[2]);
    PetscEnd();
  }
  if (bdy[3] > -Ly || bdy[4] < Ly) {
    PetscPrintf(grid.com,
        "target computational domain not a subset of source (in NetCDF file)\n"
        "  computational domain:\n");
    PetscPrintf(grid.com,
        "    need  [-Ly,Ly] contained in [bdy[3],bdy[4]],  but Ly = %5.4f km while\n"
        "    [bdy[3],bdy[4]]=[%5.4f,%5.4f] km;  ENDING ...\n",
        Ly,bdy[3],bdy[4]);
    PetscEnd();
  }
  if (bdy[6] < Lz) {
    verbPrintf(3,grid.com,
        "  WARNING: vertical dimension of target computational domain\n"
        "    not a subset of source (in NetCDF file) computational domain;\n"
        "    bdy[6] = %5.4f < Lz = %5.4f; ALLOWING ONLY 2D REGRIDDING ...\n",
        bdy[6], Lz);
    regrid_2d_only = true;
  } else {
    regrid_2d_only = false;
  }
  if (-bdy[5] < Lbz) {
    verbPrintf(3,grid.com,
        "  LIC WARNING: vertical dimension of target BEDROCK computational domain\n"
        "    not a subset of source (in NetCDF file) BEDROCK computational domain;\n"
        "    -bdy[5] = %5.4f < Lbz = %5.4f; NOT ALLOWING BEDROCK REGRIDDING ...\n",
        -bdy[5], Lbz);
    no_regrid_bedrock = true;
  } else {
    no_regrid_bedrock = false;
  }

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
largest index which is not greater than that valie.

Note that \c lic.start and \c lic.count is all that is necessary to pull the correct
data from the netCDF file, so if we implement this general scheme, the \c fstart
and \c delta entries in the struct will not be meaningful.
 */
 
  // Distances between entries (i.e. dx and dy and dz) in the netCDF file (floating point).
  delta[0] = NAN; // Delta probably will never make sense in the time dimension.
  delta[1] = (bdy[2] - bdy[1]) / (dim[1] - 1);
  delta[2] = (bdy[4] - bdy[3]) / (dim[2] - 1);

  // start[i] = index of the first needed entry in the source netCDF file; start[i] is of type int
  start[0] = dim[0] - 1; // We use the latest time
  start[1] = (int)floor((xbdy_tgt[0] - bdy[1]) / delta[1] - 0.5);
  start[2] = (int)floor((ybdy_tgt[0] - bdy[3]) / delta[2] - 0.5);
  // be sure the start[] are not too small:
  for (int m = 1; m < 3; m++) {
    if (start[m] < 0)
      start[m] = 0;
  }
  start[3] = 0; // start at base of ice
  start[4] = 0;  // start at lowest bedrock level

  fstart[0] = bdy[0];
  fstart[1] = bdy[1] + start[1] * delta[1];
  fstart[2] = bdy[3] + start[2] * delta[2];

  count[0] = 1; // Only take one time.
  count[1] = 1 + (int)ceil((xbdy_tgt[1] - fstart[1]) / delta[1]);
  count[2] = 1 + (int)ceil((ybdy_tgt[1] - fstart[2]) / delta[2]);
  // make count[] smaller if start[]+count[] would go past end of range:
  for (int m = 1; m < 3; m++) {  
    while (start[m]+count[m] > (int)(dim[m]))
      count[m]--;
  }
  count[3] = dim[3];
  count[4] = dim[4];

  nz = dim[3];
  ierr = PetscMalloc(dim[3] * sizeof(double), &(zlevs)); //CHKERRQ(ierr);
  for (size_t k = 0; k < dim[3]; k++) {
    zlevs[k] = zlevsIN[k];
  }
  nzb = dim[4];
  ierr = PetscMalloc(dim[4] * sizeof(double), &(zblevs)); //CHKERRQ(ierr);
  for (size_t k = 0; k < dim[4]; k++) {
    zblevs[k] = zblevsIN[k];
  }

  // We need a buffer for the local data, but node 0 needs to have as much
  // storage as the node with the largest block (which may be anywhere), hence
  // we perform a reduce so that node 0 has the maximum value.
  const int myzcount = (count[3] >= count[4]) ? count[3] : count[4];
  a_len = count[1] * count[2] * myzcount;
  int my_a_len = a_len;
  MPI_Reduce(&my_a_len, &(a_len), 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(a_len * sizeof(double), &(a)); //CHKERRQ(ierr);

  //return 0;  // can't return; its a constructor
}


//! Deallocate memory.
LocalInterpCtx::~LocalInterpCtx() {
  PetscErrorCode ierr;
  ierr = PetscFree(zlevs); 
  ierr = PetscFree(zblevs);
  ierr = PetscFree(a);
}


//! Print out the grid information stored in the local interpolation context.
/*!
Every processor in the communicator \c com must call this for it to work, I think.
 */
PetscErrorCode LocalInterpCtx::printGrid(MPI_Comm com) {
  PetscErrorCode ierr;

  PetscMPIInt rank;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(com,"\nLocalInterpCtx::printGrid():  rank = %d\n",
                    rank); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(com,"  delta[1,2] = %5.4f, %5.4f\n",
                    delta[1],delta[2]); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(com,"  start[0,..,4] = %d, %d, %d, %d, %d\n",
                    start[0],start[1],start[2],start[3],start[4]); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(com,"  fstart[0,1,2] = %5.4f, %5.4f, %5.4f\n",
                    fstart[0],fstart[1],fstart[2]); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(com,"  count[0,..,4] = %d, %d, %d, %d, %d\n",
                    count[0],count[1],count[2],count[3],count[4]); CHKERRQ(ierr);
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
PetscErrorCode LocalInterpCtx::printArray(MPI_Comm com) {
  PetscErrorCode ierr;

  PetscMPIInt rank;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
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
int LocalInterpCtx::kBelowHeight(const double height, MPI_Comm com) {
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
int LocalInterpCtx::kbBelowHeight(const double elevation, MPI_Comm com) {
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


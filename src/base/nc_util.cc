// Copyright (C) 2007-2008 Jed Brown and Ed Bueler
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


int nc_check(int stat) {
  if (stat)
    SETERRQ1(1, "NC_ERR: %s\n", nc_strerror(stat));
  return 0;
}


int check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    (void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    SETERRQ1(1, "NC_ERR: %s\n", nc_strerror(stat));
    //exit(1);
  }
  return 0;
}


//FIXME:  put a flag in LocalInterpCtx to indicate that only 2d stuff can be regridded; then change warning below

//! Construct a local interpolation context from arrays of parameters.
/*!
This method constructs a class from existing information already read from a NetCDF file and stored
in arrays.  It doesn't extract new information from the NetCDF file or do communication.
 */
LocalInterpCtx::LocalInterpCtx(int ncidIN, const size_t dim[], const double bdy[],
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
        "target computational domain not a subset of source (in NetCDF file) computational domain:\n");
    PetscPrintf(grid.com,
        "    need  [-Lx,Lx] contained in [bdy[1],bdy[2]],  but Lx = %5.4f km while\n"
        "    [bdy[1],bdy[2]]=[%5.4f,%5.4f] km;  ENDING ...\n",
        Lx,bdy[1],bdy[2]);
    PetscEnd();
  }
  if (bdy[3] > -Ly || bdy[4] < Ly) {
    PetscPrintf(grid.com,
        "target computational domain not a subset of source (in NetCDF file) computational domain:\n");
    PetscPrintf(grid.com,
        "    need  [-Ly,Ly] contained in [bdy[3],bdy[4]],  but Ly = %5.4f km while\n"
        "    [bdy[3],bdy[4]]=[%5.4f,%5.4f] km;  ENDING ...\n",
        Ly,bdy[3],bdy[4]);
    PetscEnd();
  }
  if (-bdy[5] < Lbz || bdy[6] < Lz) {
    verbPrintf(2,grid.com,
        "  WARNING: vertical  dimension of target computational domain not a subset of\n"
        "    source (in NetCDF file) computational domain;  either  -bdy[5] = %5.4f < Lbz = %5.4f\n"
        "    or  bdy[6] = %5.4f < Lz = %5.4f; ALLOWING ONLY 2D REGRIDDING ...\n");
    regrid_2d_only = true;
  } else {
    regrid_2d_only = false;
  }

  // limits of the processor's part of the target computational domain
  double xbdy_tgt[2] = {-Lx + dx * grid.xs, -Lx + dx * (grid.xs + grid.xm - 1)};
  double ybdy_tgt[2] = {-Ly + dy * grid.ys, -Ly + dy * (grid.ys + grid.ym - 1)};
//  double zbdy_tgt[2] = {-Lbz, Lz};

/*
To make this work with unequal spacing <i>in the horizontal dimension</i>, we have some choices.
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
 
  ncid = ncidIN;

  // Distances between entries (i.e. dx and dy and dz) in the netCDF file (floating point).
  delta[0] = NAN; // Delta probably will never make sense in the time dimension.
  delta[1] = (bdy[2] - bdy[1]) / (dim[1] - 1);
  delta[2] = (bdy[4] - bdy[3]) / (dim[2] - 1);

  // Index of the first needed entry in the source netCDF file; int type.
  start[0] = dim[0] - 1; // We use the latest time
  start[1] = (int)floor((xbdy_tgt[0] - bdy[1]) / delta[1]);
  start[2] = (int)floor((ybdy_tgt[0] - bdy[3]) / delta[2]);
  start[3] = 0; // start at base of ice
  start[4] = 0;  // start at lowest bedrock level

  fstart[0] = bdy[0];
  fstart[1] = bdy[1] + start[1] * delta[1];
  fstart[2] = bdy[3] + start[2] * delta[2];

  count[0] = 1; // Only take one time.
  count[1] = (int)ceil((xbdy_tgt[1] - fstart[1]) / delta[1] + 1);
  count[2] = (int)ceil((ybdy_tgt[1] - fstart[2]) / delta[2] + 1);
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
  ierr = PetscMalloc(a_len * sizeof(float), &(a)); //CHKERRQ(ierr);

}


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


NCTool::NCTool() {
  //FIXME: does there need to be a default MaskInterp?
  myMaskInterp = PETSC_NULL;
}


PetscErrorCode NCTool::set_MaskInterp(MaskInterp *mi_in) {
  myMaskInterp = mi_in;
  return 0;
}


//! Put a \c DA -managed local \c Vec \c v into a variable in a NetCDF file; a global \c Vec \c g is used for storage space.
PetscErrorCode NCTool::put_local_var(const IceGrid *grid, int ncid, const int var_id, nc_type type,
                                     DA da, Vec v, Vec g, const int *s, const int *c,
                                     int dims, void *a_mpi, int a_size) {

  PetscErrorCode ierr;
  ierr = DALocalToGlobal(da, v, INSERT_VALUES, g); CHKERRQ(ierr);
  ierr = put_global_var(grid, ncid, var_id, type, da, g, s, c, dims, a_mpi, a_size); CHKERRQ(ierr);
  return 0;
}


//! Put a \c DA -managed global \c Vec \c g into a variable in a NetCDF file.  \e In \e parallel.
PetscErrorCode NCTool::put_global_var(const IceGrid *grid, int ncid, const int var_id, nc_type type,
                                      DA da, Vec g, const int *s, const int *c,
                                      int dims, void *a_mpi, int a_size) {
  const int lim_tag = 1; // MPI tag for limits block
  const int var_tag = 2; // MPI tag for data block
  const int sc_size = 8;
  PetscErrorCode ierr;
  MPI_Status mpi_stat;
  int stat;
  int sc[sc_size]; // buffer to hold both `s' and `c'
  size_t sc_nc[sc_size];
  float *a_float = NULL;
  unsigned char *a_uchar = NULL;

  if (type == NC_FLOAT) {
    a_float = (float *)a_mpi;
  } else if (type == NC_BYTE) {
    a_uchar = (unsigned char *)a_mpi;
  } else {
    SETERRQ(1, "Unsupported type.");
  }

  for (int i = 0; i < 2 * dims; i++) {
    sc[i] = (i < dims) ? s[i] : c[i - dims];
  }

  int b_len = 1;
  for (int i = 0; i < dims; i++) b_len *= c[i];

  // convert IceModel Vec containing PetscScalar to array of float or char for NetCDF
  PetscScalar *a_petsc;
  ierr = VecGetArray(g, &a_petsc); CHKERRQ(ierr);
  for (int i = 0; i < b_len; i++) {
    if (type == NC_FLOAT) {
      a_float[i] = (float)a_petsc[i];
    } else if (type == NC_BYTE) {
      a_uchar[i] = (unsigned char)a_petsc[i];
    } else {
      SETERRQ(1, "Unsupported type.");
    }
  }
  ierr = VecRestoreArray(g, &a_petsc); CHKERRQ(ierr);

  if (grid->rank == 0) { // on rank 0 processor, receive messages from every other
                         //    processor, then write it out to the NC file
    for (int proc = 0; proc < grid->size; proc++) { // root will write itself last
      if (proc != 0) {
        MPI_Recv(sc, sc_size, MPI_INT, proc, lim_tag, grid->com, &mpi_stat);
        if (type == NC_FLOAT) {
          MPI_Recv(a_float, a_size, MPI_FLOAT, proc, var_tag, grid->com, &mpi_stat);
        } else if (type == NC_BYTE) {
          MPI_Recv(a_uchar, a_size, MPI_UNSIGNED_CHAR, proc, var_tag, grid->com, &mpi_stat);
        } else {
          SETERRQ(1, "Unsupported type.");
        }
      }

      /* {
        printf("[%1d] writing %4d [", proc, var_id);
        for (int i = 0; i < 2 * dims; i++) {
          if (i == dims) printf("] [");
          printf("%5d", sc[i]);
        }
        printf("]\n");
      } */

      for (int i = 0; i < 2 * dims; i++) sc_nc[i] = (size_t)sc[i]; // we need size_t
      if (type == NC_FLOAT) {
        stat = nc_put_vara_float(ncid, var_id, &sc_nc[0], &sc_nc[dims], a_float);
      } else if (type == NC_BYTE) {
        stat = nc_put_vara_uchar(ncid, var_id, &sc_nc[0], &sc_nc[dims], a_uchar);
      } else {
        SETERRQ(1, "Unsupported type.");
      }
      CHKERRQ(check_err(stat,__LINE__,__FILE__));
    }
  } else {  // all other processors send to rank 0 processor
    MPI_Send(sc, 2 * dims, MPI_INT, 0, lim_tag, grid->com);
    if (type == NC_FLOAT) {
      MPI_Send(a_float, a_size, MPI_FLOAT, 0, var_tag, grid->com);
    } else if (type == NC_BYTE) {
      MPI_Send(a_uchar, a_size, MPI_UNSIGNED_CHAR, 0, var_tag, grid->com);
    } else {
      SETERRQ(1, "Unsupported type.");
    }
  }
  return 0;
}


//! Put the variable for a dimension in a NetCDF file.  Uses starting value and a spacing for regularly-spaced values.
/*!
Note the variable corresponding to a dimension is always \c double in a PISM NetCDF file.
 */
PetscErrorCode NCTool::put_dimension_regular(int ncid, int v_id, int len, double start, double delta) {
  PetscErrorCode ierr;
  int stat;
  double *v;

  ierr = PetscMalloc(len * sizeof(double), &v); CHKERRQ(ierr);
  for (int i = 0; i < len; i++) {
    v[i] = start + i * delta;
  }
  stat = nc_put_var_double(ncid, v_id, v); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  ierr = PetscFree(v); CHKERRQ(ierr);

  return 0;
}


//! Put the variable for a dimension in a NetCDF file.  Makes no assumption about spacing.
/*!
Note the variable corresponding to a dimension is always \c double in a PISM NetCDF file.
 */
PetscErrorCode NCTool::put_dimension(int ncid, int v_id, int len, PetscScalar *vals) {
  PetscErrorCode ierr;
  int stat;
  double *v;

  ierr = PetscMalloc(len * sizeof(double), &v); CHKERRQ(ierr);
  for (int i = 0; i < len; i++) {
    v[i] = (double)vals[i];
  }
  stat = nc_put_var_double(ncid, v_id, v); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  ierr = PetscFree(v); CHKERRQ(ierr);
  return 0;
}


//! Read the first and last values, and the lengths, of the x,y,z,zb dimensions from a NetCDF file.  Read the last t.
/*!
Correspondence between the parameters read into arrays \c dim[] and \c bdy[] by this procedure,
and the values in \c IceGrid:
  - <tt>dim[0]</tt> is the number of t values in the NetCDF file; it is ignored for many purposes
  - <tt>dim[1]</tt> is the length of the x dimension; it becomes <tt>grid.Mx</tt>
  - <tt>dim[2]</tt> is the length of the y dimension; it becomes <tt>grid.My</tt>
  - <tt>dim[3]</tt> is the length of the z (ice) dimension; it becomes <tt>grid<tt>dim[0..2]</tt> 
             and <tt>bdy[0..4]</tt>.Mz</tt>
  - <tt>dim[4]</tt> is the length of the zb (bedrock) dimension; it becomes <tt>grid.Mbz</tt>
  - <tt>bdy[0]</tt> is current time (and becomes <tt>grid.year</tt>)
  - <tt>bdy[2]</tt>=-<tt>bdy[1]</tt> is \f$x\f$ half-length of computational domain; becomes <tt>grid.Lx</tt>
  - <tt>bdy[4]</tt>=-<tt>bdy[3]</tt> is \f$y\f$ half-length of computational domain; becomes <tt>grid.Ly</tt>
  - <tt>bdy[5]</tt> is the negative of the thickness (positive) of bedrock layer (for thermal model); 
             <tt>grid.Lbz</tt>=-<tt>bdy[5]</tt>
  - <tt>bdy[6]</tt> is thickness (positive) of ice layer; becomes <tt>grid.Lz</tt>
 */
PetscErrorCode NCTool::get_dims_limits_lengths(int ncid, size_t dim[], double bdy[], MPI_Comm com) {
  PetscErrorCode ierr;
  
  // first fill dim[0..2] and bdy[0..4]
  ierr = get_dims_limits_lengths_2d(ncid, dim, bdy, com); CHKERRQ(ierr);

  PetscMPIInt rank;
  int stat;
  int z_dim, zb_dim;
  int z_id, zb_id;
  size_t z_len, zb_len;

  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);

  if (rank == 0) {
    stat = nc_inq_dimid(ncid, "z", &z_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "zb", &zb_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_inq_dimlen(ncid, z_dim, &z_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimlen(ncid, zb_dim, &zb_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    dim[3] = z_len; 
    dim[4] = zb_len;

    stat = nc_inq_varid(ncid, "z", &z_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_varid(ncid, "zb", &zb_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    size_t z_bdy[] = {0, z_len - 1}; // Start at 0 in `zb', end at z_len - 1 of `z'

    stat = nc_get_var1_double(ncid, zb_id, &z_bdy[0], &bdy[5]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_var1_double(ncid, z_id, &z_bdy[1], &bdy[6]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  MPI_Bcast(dim, 5, MPI_LONG, 0, com);
  MPI_Bcast(bdy, 7, MPI_DOUBLE, 0, com);

  return 0;
}


//! Read the first and last values, and the lengths, of the x,y dimensions from a NetCDF file.  Read the last t.
/*!
See get_dims_limits_lengths() for the meaning of the returned entries of \c dim[] and \c bdy[].

This version only sets sets <tt>dim[0..2]</tt> and <tt>bdy[0..4]</tt>.  It ignors the 3d information,
if that information is present at all (e.g. when calling bootstrapFromFile_netCDF() it is not usually present).  
This method does not modify <tt>dim[3,4]</tt> and <tt>bdy[5,6]</tt>, which may even be invalid.
 */
PetscErrorCode NCTool::get_dims_limits_lengths_2d(int ncid, size_t dim[], double bdy[], MPI_Comm com) {
  PetscErrorCode ierr;

  ierr = get_last_time(ncid, &(bdy[0]), com); CHKERRQ(ierr);

  PetscMPIInt rank;
  int stat;
  int t_dim, x_dim, y_dim;
  int x_id, y_id;
  size_t t_len, x_len, y_len;

  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);

  if (rank == 0) {
    stat = nc_inq_dimid(ncid, "t", &t_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "x", &x_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "y", &y_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_inq_dimlen(ncid, t_dim, &t_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimlen(ncid, x_dim, &x_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimlen(ncid, y_dim, &y_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    dim[0] = t_len; 
    dim[1] = x_len; 
    dim[2] = y_len;

    stat = nc_inq_varid(ncid, "x", &x_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_varid(ncid, "y", &y_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    // get extent of grid by looking for last and first of variables x,y
    size_t x_bdy[] = {0, x_len - 1};
    size_t y_bdy[] = {0, y_len - 1};
    stat = nc_get_var1_double(ncid, x_id, &x_bdy[0], &bdy[1]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_var1_double(ncid, x_id, &x_bdy[1], &bdy[2]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_var1_double(ncid, y_id, &y_bdy[0], &bdy[3]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_var1_double(ncid, y_id, &y_bdy[1], &bdy[4]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  
  MPI_Bcast(dim, 5, MPI_LONG, 0, com);
  MPI_Bcast(bdy, 7, MPI_DOUBLE, 0, com);

  return 0;
}


//! Read the last value of the time variable t from a NetCDF file.
PetscErrorCode NCTool::get_last_time(int ncid, double *time, MPI_Comm com) {
  PetscErrorCode ierr;
  PetscMPIInt rank;
  int stat;
  int t_dim;
  int t_id;
  size_t t_len;
  double mytime;

  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);

  if (rank == 0) {
    stat = nc_inq_dimid(ncid, "t", &t_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimlen(ncid, t_dim, &t_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_varid(ncid, "t", &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    size_t t_end = t_len - 1;
    stat = nc_get_var1_double(ncid, t_id, &t_end, &mytime);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  MPI_Bcast(&mytime, 1, MPI_DOUBLE, 0, com);

  *time = mytime;
  return 0;
}


// FIXME: shouldn't this be expecting type double, not type float, for dimensions?  or checking?
PetscErrorCode NCTool::get_ends_1d_var(int ncid, int vid, PetscScalar *gfirst, PetscScalar *glast,
                                       MPI_Comm com) {
  PetscErrorCode ierr;
  int stat;
  PetscScalar first, last;
  float       *f = NULL;
  int         *g = NULL;

  PetscMPIInt rank;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);

  if (rank == 0) {
    int dimids[NC_MAX_VAR_DIMS];
    int ndims, natts;
    nc_type xtype;
    char name[NC_MAX_NAME+1];
    stat = nc_inq_var(ncid, vid, name, &xtype, &ndims, dimids, &natts); CHKERRQ(nc_check(stat));
    if (ndims != 1) {
      SETERRQ2(1, "getFirstLast: number of dimensions = %d for %s\n",
               ndims, name);
    }

    // In the netCDF file,
    // we index 0:M in the x direction and 0:N in the y direction.  A
    // location $(i,j) \in [0,M] \times [0,N]$ is addressed as [i*N + j]
    size_t M;
    stat = nc_inq_dimlen(ncid, dimids[0], &M); CHKERRQ(nc_check(stat));

    switch (xtype) {
      case NC_INT:
        g = new int[M];
        stat = nc_get_var_int(ncid, vid, g); CHKERRQ(nc_check(stat));
        break;
      case NC_FLOAT:
        f = new float[M];
        stat = nc_get_var_float(ncid, vid, f); CHKERRQ(nc_check(stat));
        break;
      default:
        SETERRQ1(1, "NC_VAR `%s' not of type NC_INT or NC_FLOAT.\n", name);
    }

    if (g != NULL) {
      first = g[0];
      last = g[M-1];
    } else if (f != NULL) {
      first = f[0];
      last = f[M-1];
    } else {
      SETERRQ(1, "This should not happen.\n");
    }
  } else {
    first = 1.0e30;
    last = -1.0e30;
  }

  ierr = PetscGlobalMin(&first,gfirst,com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&last,glast,com); CHKERRQ(ierr);
  return 0;
}


//! Read in the variables \c z and \c zb from the NetCDF file; <i>do not</i> assume they are equally-spaced.
/*!
Arrays z_read[] and zb_read[] must be pre-allocated arrays of length at least z_len, zb_len, respectively.

Get values of z_len, zb_len by using get_dims_limits_lengths() before this routine; z_len = dim[3] and zb_len = dim[4].
 */
PetscErrorCode NCTool::get_vertical_dims(int ncid, const int z_len, const int zb_len,
                                         double z_read[], double zb_read[], MPI_Comm com) {
  PetscErrorCode ierr;
  PetscMPIInt rank;
  int stat;
  int z_id, zb_id;
  size_t zeroST  = 0,
         zlenST  = (size_t) z_len,
         zblenST = (size_t) zb_len;

  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);

  if (rank == 0) {
    stat = nc_inq_varid(ncid, "z", &z_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_varid(ncid, "zb", &zb_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_get_vara_double(ncid, z_id, &zeroST, &zlenST, z_read);
             CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_vara_double(ncid, zb_id, &zeroST, &zblenST, zb_read);
             CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  MPI_Bcast(z_read, z_len, MPI_DOUBLE, 0, com);
  MPI_Bcast(zb_read, zb_len, MPI_DOUBLE, 0, com);
  return 0;
}


//! Read a NetCDF variable into a \c DA -managed local \c Vec \c v; a global \c Vec \c g is used for storage.
PetscErrorCode NCTool::get_local_var(const IceGrid *grid, int ncid, const char *name, nc_type type,
                                     DA da, Vec v, Vec g, const int *s, const int *c,
                                     int dims, void *a_mpi, int a_size) {

  PetscErrorCode ierr;
  ierr = get_global_var(grid, ncid, name, type, da, g, s, c, dims, a_mpi, a_size); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, g, INSERT_VALUES, v); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, g, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


//! Read a variable in a NetCDF file into a \c DA -managed global \c Vec \c g.  \e In \e parallel.
PetscErrorCode NCTool::get_global_var(const IceGrid *grid, int ncid, const char *name, nc_type type,
                                      DA da, Vec g, const int *s, const int *c,
                                      int dims, void *a_mpi, int a_size) {
  const int req_tag = 1; // MPI tag for request block
  const int var_tag = 2; // MPI tag for data block
  const int sc_size = 8;
  PetscErrorCode ierr;
  MPI_Status mpi_stat;
  int stat;
  int sc[sc_size]; // buffer to hold both `s' and `c'
  size_t sc_nc[sc_size];
  float *a_float = NULL;
  unsigned char *a_uchar = NULL;

  if (type == NC_FLOAT) {
    a_float = (float *)a_mpi;
  } else if (type == NC_BYTE) {
    a_uchar = (unsigned char *)a_mpi;
  } else {
    SETERRQ(1, "Unsupported type.");
  }

  for (int i = 0; i < 2 * dims; i++) {
    sc[i] = (i < dims) ? s[i] : c[i - dims];
  }

  if (grid->rank == 0) {
    int sc0[sc_size];
    for (int i = 0; i < sc_size; i++) sc0[i] = sc[i]; // root needs to save its range
    for (int proc = grid->size - 1; proc >= 0; proc--) { // root will read itself last
      if (proc == 0) {
        for (int i = 0; i < sc_size; i++) sc[i] = sc0[i];
      } else {
        MPI_Recv(sc, sc_size, MPI_INT, proc, req_tag, grid->com, &mpi_stat);
      }
      for (int i = 0; i < 2 * dims; i++) sc_nc[i] = (size_t)sc[i]; // we need size_t

      /* {
        printf("[%1d] reading %10s [", proc, name);
        for (int i = 0; i < 2 * dims; i++) {
          if (i == dims) printf("] [");
          printf("%5d", (int)sc_nc[i]);
        }
        printf("]\n");
      } */

      int var_id;
      stat = nc_inq_varid(ncid, name, &var_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
      if (type == NC_FLOAT) {
        stat = nc_get_vara_float(ncid, var_id, &sc_nc[0], &sc_nc[dims], a_float);
      } else if (type == NC_BYTE) {
        stat = nc_get_vara_uchar(ncid, var_id, &sc_nc[0], &sc_nc[dims], a_uchar);
      } else {
        SETERRQ(1, "Unsupported type.");
      }
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      if (proc != 0) {
        int b_len = 1;
        for (int i = dims; i < 2 * dims; i++) b_len *= sc[i];
        if (type == NC_FLOAT) {
          MPI_Send(a_float, b_len, MPI_FLOAT, proc, var_tag, grid->com);
        } else if (type == NC_BYTE) {
          MPI_Send(a_uchar, b_len, MPI_UNSIGNED_CHAR, proc, var_tag, grid->com);
        } else {
          SETERRQ(1, "Unsupported type.");
        }
      }
    }
  } else {
    MPI_Send(sc, 2 * dims, MPI_INT, 0, req_tag, grid->com);
    if (type == NC_FLOAT) {
      MPI_Recv(a_float, a_size, MPI_FLOAT, 0, var_tag, grid->com, &mpi_stat);
    } else if (type == NC_BYTE) {
      MPI_Recv(a_uchar, a_size, MPI_UNSIGNED_CHAR, 0, var_tag, grid->com, &mpi_stat);
    } else {
      SETERRQ(1, "Unsupported type.");
    }
  }

  int b_len = 1;
  for (int i = dims; i < 2 * dims; i++) b_len *= sc[i];
  PetscScalar *a_petsc;
  ierr = VecGetArray(g, &a_petsc); CHKERRQ(ierr);
  for (int i = 0; i < b_len; i++) {
    if (type == NC_FLOAT) {
      a_petsc[i] = (PetscScalar)a_float[i];
    } else if (type == NC_BYTE) {
      a_petsc[i] = (PetscScalar)a_uchar[i];
    } else {
      SETERRQ(1, "Unsupported type.");
    }
  }

  ierr = VecRestoreArray(g, &a_petsc); CHKERRQ(ierr);
  return 0;
}


//! Find a 2D or 3D variable in a NetCDF file and regrid it using a global \c Vec for storage.
/*!
Simply calls regrid_global_var().  Then transfers the global \c Vec \c g to the local \c Vec \c vec.
 */
PetscErrorCode NCTool::regrid_local_var(const char *name, int dim_flag, LocalInterpCtx &lic,
                                        IceGrid &grid, DA da, Vec vec, Vec g, bool useMaskInterp) {
  PetscErrorCode ierr;
  ierr = regrid_global_var(name, dim_flag, lic, grid, da, g, useMaskInterp); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, g, INSERT_VALUES, vec); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, g, INSERT_VALUES, vec); CHKERRQ(ierr);
  return 0;
}


//! Find a 2D or 3D variable in a NetCDF file and regrid it onto the current grid.  \e In \e parallel.
/*!
We need to move a 2D or 3D variable from within a NetCDF file, with its "source" grid, to the
current grid, the "target" grid.  The variable on the source grid is only read from processor 
zero.  The target grid is spread across all processors.

The source grid may be coarser or finer than the target
grid.  In fact the source grid may be coarser in one dimension and finer in another.

We <i>require</i>,
however, that the source grid have greater extent than the target.  That is, the values
of the parameters \c Lx, \c Ly, \c Lz, and \c Lbz for the source grid must exceed those of the
target grid.

This method checks whether the single character flag \c c is in the string \c vars.

Note that \c dim_flag is 2 for 2-D quantities, 3 for 3-D ice quantities, and 4 for 3-D
bedrock quantities.

We interpolate as needed for functions of two or three variables.  The interpolation is on 
rectangles and rectangular solids, respectively, and the interpolant families have four and eight
free parameters, respectively.

In two variables we do interpolation using 
functions like \f$F(x,y) = A + Bx + Cy + Dxy\f$.  In particular, we take four values of the source 
variable \f$a_{ij}, a_{i+1,j}, a_{i,j+1}, a_{i+1,j+1}\f$ 
at the corners of a rectangle in the source grid.  
The corners themselves would be \f$(x_i,y_j)\f$, \f$(x_i+\Delta x,y_j)\f$, \f$(x_i,y_j+\Delta y)\f$, 
\f$(x_i+\Delta x,y_j+\Delta y)\f$, where \f$\Delta x\f$ is stored in \c lic.delta[1] and
\f$\Delta y\f$ is stored in \c lic.delta[2].  Denote the four values from the source grid as 
\f$a_{mm}, a_{pm}, a_{mp}, a_{pp}\f$.  Then
   \f{align*}
      F(x,y) &= \left[a_{mm} * (1 - jj) + a_{mp} * jj\right] * (1 - ii) \\
             &\qquad + \left[a_{pm} * (1 - jj) + a_{pp} * jj\right] * ii
   \f}
where \f$ii = (x - x_i)/\Delta x\f$ and \f$jj = (y - y_j)/\Delta y\f$.  Note that the compact
way to write this, and the way it is implemented, is to define
   \f{align*}
      a_m &:= a_{mm} * (1 - jj) + a_{mp} * jj, \\
      a_p &:= a_{pm} * (1 - jj) + a_{pp} * jj \\
   \f}
and then compute
   \f{align*}
      F(x,y) &= a_m * (1 - ii) + a_p * ii
   \f}

In three variables we do interpolation in a completely analogous way using 
functions like \f$F(x,y,z) = A + Bx + Cy + Dxy + Pz + Qxz + Ryz + Sxyz\f$.  With analogous notation
for the values of the source variable at the \e eight corners of the rectangular solid source grid 
cell, the first step is to compute values at the four corners of the rectangle which comes 
(conceptually) from projecting the solid source grid cell into the plane:
   \f{align*}
      a_{mm} &:= a_{mmm} * (1 - kk) + a_{mmp} * kk, \\
      a_{mp} &:= a_{mpm} * (1 - kk) + a_{mpp} * kk, \\
      a_{pm} &:= a_{pmm} * (1 - kk) + a_{pmp} * kk, \\
      a_{pp} &:= a_{ppm} * (1 - kk) + a_{ppp} * kk. \\
   \f}
Here \f$kk\f$ is the fraction of the vertical space above level \c k, between the two grid
levels on the source grid.  That is,

Then we just do the two variable interpolation as before, finding \f$a_{m}\f$ and \f$a_p\f$ before
computing \f$F(x,y,z)\f$.
 */
PetscErrorCode NCTool::regrid_global_var(const char *name, int dim_flag, LocalInterpCtx &lic,
                                         IceGrid &grid, DA da, Vec g, bool useMaskInterp) {
  PetscErrorCode ierr;

  const int req_tag = 1; // MPI tag for request block
  const int var_tag = 2; // MPI tag for data block
  const int sc_len = 8;
  MPI_Status mpi_stat;
  int stat, dims;
  int sc[sc_len];

  switch (dim_flag) {
    case 2:
      dims = 3; // time, x, y
      break;
    case 3:
    case 4:
      dims = 4; // time, x, y, {z|zb}
      if (lic.regrid_2d_only) {
        SETERRQ(2, "regrid_2d_only is set, so dim_flag == 3 or 4 is not allowed");
      }
      break;
    default:
      SETERRQ(1, "Invalid value for `dim_flag'.");
  }

  // pack start and count into a single buffer
  for (int i = 0; i < 4; i++) sc[i] = lic.start[i];
  for (int i = 0; i < 4; i++) sc[4 + i] = lic.count[i];

  // At this point, sc[] is set up correctly for ice 3-D quantities.  In order
  // to avoid duplicating lots of code below, we play some tricks here.
  if (dim_flag == 2) { // 2-D quantity
    // Treat it as a 3-D quantity with extent 1 in the z-direction.  The netCDF
    // read actually will not use this information (it is passed in, but it
    // won't index these entries in the array, because it knows how many
    // dimensions it there are), but sc[7] will be used when node 0 computes how
    // much data it has to send back.  sc[3] should never be used in this case.
    sc[3] = 0; sc[7] = 1;
  } else if (dim_flag == 4) { // Bedrock quantity
    // Just fill in the appropriate values
    sc[3] = lic.start[4]; sc[7] = lic.count[4];
  }

  if (grid.rank == 0) {
    int var_id;
    stat = nc_inq_varid(lic.ncid, name, &var_id);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    // Node 0 will service all the other nodes before itself.  We need to save
    // sc[] so that it knows how to get its block at the end.
    int sc0[sc_len];
    for (int i = 0; i < sc_len; i++) sc0[i] = sc[i];

    for (int proc = grid.size - 1; proc >= 0; proc--) {
      if (proc == 0) {// Get the bounds.
        for (int i = 0; i < sc_len; i++) sc[i] = sc0[i];
      } else {
        MPI_Recv(sc, sc_len, MPI_INT, proc, req_tag, grid.com, &mpi_stat);
      }

      // It is not safe to cast memory.  In particular on amd64 int and size_t are
      // different sizes.  Since netCDF uses size_t for the offsets, we need to here too.
      size_t sc_nc[sc_len];
      for (int i = 0; i < sc_len; i++) sc_nc[i] = (size_t)sc[i]; // we need size_t

      // Actually read the block into the buffer.
      stat = nc_get_vara_float(lic.ncid, var_id, &sc_nc[0], &sc_nc[4], lic.a);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      /*{
        printf("lic.ncid, var_id = %d %d\n", lic.ncid, var_id);
        printf("a[] ni {%f %f %f %f}\n", lic.a[50], lic.a[150], lic.a[250], lic.a[350]);

        const int blen = 9;
        float *buf = (float *)malloc(blen * sizeof(float));
        sc_nc[5] = sc_nc[6] = 3; sc_nc[7] = 1;
        stat = nc_get_vara_float(lic.ncid, var_id, &sc_nc[0], &sc_nc[4], buf);
        CHKERRQ(check_err(stat,__LINE__,__FILE__));
        printf("buf = ");
        for (int i = 0; i < blen; i++) printf(" %7.2e ", buf[i]);
        printf("\n");
        free(buf);
      }*/

      // Find out how big the buffer actually is.  Remember that node 0 has a
      // buffer that will only be filled by if the process it is serving has a
      // maximal sized local domain.
      int a_len = 1;
      for (int i = 0; i < dims; i++) a_len *= sc[4 + i];

      // send the filled buffer
      if (proc != 0) {
        MPI_Send(lic.a, a_len, MPI_FLOAT, proc, var_tag, grid.com);
      }
    }
  } else { // not process 0:
    MPI_Send(sc, sc_len, MPI_INT, 0, req_tag, grid.com);  // send out my bounds
    MPI_Recv(lic.a, lic.a_len, MPI_FLOAT, 0, var_tag, grid.com, &mpi_stat); // get back filled buffer
  }

  // At this point, the buffer lic.a[] should contain lic.a_len floats.  This is the 
  // local processor's part of the source variable.  
  // That is, it should be enough of the source variable so that \e interpolation
  // (not extrapolation) onto the local processor's part of the target grid is possible.
  //ierr = lic.printArray(grid.com); CHKERRQ(ierr);
  
  // indexing parameters
  const int ycount = lic.count[2];
  int myMz = 1, zcount = 1; // indexing trick so that we don't have to duplicate code for the 2-D case.
  if (dim_flag == 3) {
    myMz = grid.Mz;
    zcount = lic.count[3];
  } else if (dim_flag == 4) {
    myMz = grid.Mbz;
    zcount = lic.count[4];
  }

  // We'll work with the raw storage here so that the array we are filling is
  // indexed the same way as the buffer we are pulling from (lic.a)
  PetscScalar *vec_a;
  ierr = VecGetArray(g, &vec_a); CHKERRQ(ierr);

  for (int i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (int j = grid.ys; j < grid.ys + grid.ym; j++) {

      for (int k = 0; k < myMz; k++) {
        // location (x,y,z) is in target computational domain
        const double x = -grid.Lx + i * grid.dx,
                     y = -grid.Ly + j * grid.dy,
                     z = (dim_flag == 3) ? grid.zlevels[k] : grid.zblevels[k];

        // We need to know how the point (x,y,z) sits within the local block we
        // pulled from the netCDF file.  This part is specialized to a regular
        // grid.  In particular floor(ic) is the index of the 'left neighbor'
        // and ceil(ic) is the index of the 'right neighbor'.  For the irregular
        // case, we would also need the physical locations of these points, so
        // that we could compute ii,jj,kk (described below).  We should just
        // compute the 'left index' and 'right index' explicitly (but in the
        // same way as we find the domain bounds; see note below)
        const double ic = (x - lic.fstart[1]) / lic.delta[1];
        const double jc = (y - lic.fstart[2]) / lic.delta[2];

        double a_mm, a_mp, a_pm, a_pp;  // filled differently in 2d and 3d cases

        if (dim_flag == 3 || dim_flag == 4) {
          // get the index into the source grid, for just below the level z
          const int kc = (dim_flag == 3) ? lic.kBelowHeight(z,grid.com)
                                           : lic.kbBelowHeight(z,grid.com);

          // We pretend that there are always 8 neighbors.  And compute the
          // indices into the buffer for those neighbors.  
          // It is important to
          // note that floor(ic) + 1 = ceil(ic) does not hold when ic is an
          // integer.  Computation of the domain (in constructor of LocalInterpCtx; note
          // that lic.count uses ceil) must be done in a compatible way,
          // otherwise we can index improperly here.  When it is in the
          // interior, it should not matter (since the coefficient of the
          // erroneous point will be 0), but if it occurs at the end of the
          // buffer, it will index off the array, possibly causing a
          // segmentation fault.
          //
          // In light of this degenerate case, we observe that not all of these
          // neighbors are necessarily unique, but doing it this way enables us
          // to not handle all the cases explicitly.
          // 
          // (These comments do not apply to the z case.)
          const int mmm = ((int)floor(ic) * ycount + (int)floor(jc)) * zcount + kc;
          const int mmp = ((int)floor(ic) * ycount + (int)floor(jc)) * zcount + kc + 1;
          const int mpm = ((int)floor(ic) * ycount + (int)ceil(jc)) * zcount + kc;
          const int mpp = ((int)floor(ic) * ycount + (int)ceil(jc)) * zcount + kc + 1;
          const int pmm = ((int)ceil(ic) * ycount + (int)floor(jc)) * zcount + kc;
          const int pmp = ((int)ceil(ic) * ycount + (int)floor(jc)) * zcount + kc + 1;
          const int ppm = ((int)ceil(ic) * ycount + (int)ceil(jc)) * zcount + kc;
          const int ppp = ((int)ceil(ic) * ycount + (int)ceil(jc)) * zcount + kc + 1;

          // We know how to index the neighbors, but we don't yet know where the
          // point lies within this box.  This is represented by kk in [0,1].
          // For the irregular case, with left index km and right index kp in the source grid, we
          // would have
          //   kk = (km == kp) ? 0.0 : (z - Z(km)) / (Z(kp) - Z(km))
          // where Z(.) are the physical coordinates on the source grid.  Note
          // that any value in [0,1] would be okay when km == kp.
          const double zkc = (dim_flag == 3) ? lic.zlevs[kc] : lic.zblevs[kc];
          double dz;
          if (dim_flag == 3) {
            if (kc == zcount - 1) {
              dz = lic.zlevs[kc] - lic.zlevs[kc-1];
            } else {
              dz = lic.zlevs[kc+1] - lic.zlevs[kc];
            }
          } else {
            if (kc == zcount - 1) {
              dz = lic.zblevs[kc] - lic.zblevs[kc-1];
            } else {
              dz = lic.zblevs[kc+1] - lic.zblevs[kc];
            }
          }
          const double kk = (z - zkc) / dz;

          // linear interpolation in the z-direction
          a_mm = (double)(lic.a[mmm]) * (1.0 - kk) + (double)(lic.a[mmp]) * kk;
          a_mp = (double)(lic.a[mpm]) * (1.0 - kk) + (double)(lic.a[mpp]) * kk;
          a_pm = (double)(lic.a[pmm]) * (1.0 - kk) + (double)(lic.a[pmp]) * kk;
          a_pp = (double)(lic.a[ppm]) * (1.0 - kk) + (double)(lic.a[ppp]) * kk;
        } else {
          // we don't need to interpolate vertically for the 2-D case
          a_mm = (double)(lic.a[(int)floor(ic) * ycount + (int)floor(jc)]);
          a_mp = (double)(lic.a[(int)floor(ic) * ycount + (int)ceil(jc)]);
          a_pm = (double)(lic.a[(int)ceil(ic) * ycount + (int)floor(jc)]);
          a_pp = (double)(lic.a[(int)ceil(ic) * ycount + (int)ceil(jc)]);
        }

        const double jj = jc - floor(jc);

        // interpolate in y direction
        const double a_m = a_mm * (1.0 - jj) + a_mp * jj;
        const double a_p = a_pm * (1.0 - jj) + a_pp * jj;

        const double ii = ic - floor(ic);
        int index = ((i - grid.xs) * grid.ym + (j - grid.ys)) * myMz + k;

        // index into the new array and interpolate in x direction
        vec_a[index] = a_m * (1.0 - ii) + a_p * ii;
        if (useMaskInterp) {
          if (myMaskInterp == PETSC_NULL) {
            SETERRQ(3,"NCTool::myMaskInterp needed, but not initialized correctly");
          }
          double val = vec_a[index];
          if (   (myMaskInterp->number_allowed == 1)
              || (val <= (double)(myMaskInterp->allowed_levels[0])) ) {
            val = (double)(myMaskInterp->allowed_levels[0]);
          } else {
            int k=1;
            while (k < myMaskInterp->number_allowed) {
             const double mid = (  (double)(myMaskInterp->allowed_levels[k-1])
                                 + (double)(myMaskInterp->allowed_levels[k])   ) / 2.0;
              if (val < mid) {
                val = (double)(myMaskInterp->allowed_levels[k-1]);
                break;
              }
              k++;
            }
            if (k >= myMaskInterp->number_allowed) {
              val = (double)(myMaskInterp->allowed_levels[k-1]);
            }
          }
        }
        // done with the point at (x,y,z)
      }
    }
  }

  ierr = VecRestoreArray(g, &vec_a); CHKERRQ(ierr);

  return 0;
}

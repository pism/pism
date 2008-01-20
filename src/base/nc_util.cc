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


//! Nothin' here because NCTool is a class created so that doxygen can document it.
NCTool::NCTool() {
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


//! Put a \c DA -managed global \c Vec \c g into a variable in a NetCDF file.
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


//! Read the first and last values, and the lengths, of the x,y,z,zb dimensions from a NetCDF file.  Also read final time.
/*!
Correspondence between the parameters to this procedure and the values in \c IceGrid:
  - <tt>bdy[0]</tt> is current time and becomes <tt>grid.year</tt>
  - <tt>-bdy[1]</tt>=<tt>bdy[2]</tt> is \f$x\f$ half-length of computational domain and becomes <tt>grid.Lx</tt>
  - <tt>-bdy[3]</tt>=<tt>bdy[4]</tt> is \f$y\f$ half-length of computational domain and becomes <tt>grid.Ly</tt>
  - <tt>-bdy[5]</tt> is thickness (positive) of bedrock layer (for thermal model); becomes <tt>grid.Lbz</tt>
  - <tt>bdy[6]</tt> is thickness (positive) of ice layer and becomes <tt>grid.Lz</tt>
 */
PetscErrorCode NCTool::get_dims_limits_lengths(int ncid, size_t dim[], double bdy[], MPI_Comm com) {
  PetscErrorCode ierr;
  PetscMPIInt rank;
  int stat;
  int t_dim, x_dim, y_dim, z_dim, zb_dim;
  int t_id, x_id, y_id, z_id, zb_id;
  size_t t_len, x_len, y_len, z_len, zb_len;

  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);

  if (rank == 0) {
    stat = nc_inq_dimid(ncid, "t", &t_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "x", &x_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "y", &y_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "z", &z_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "zb", &zb_dim); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_inq_dimlen(ncid, t_dim, &t_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimlen(ncid, x_dim, &x_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimlen(ncid, y_dim, &y_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimlen(ncid, z_dim, &z_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimlen(ncid, zb_dim, &zb_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    dim[0] = t_len; dim[1] = x_len; dim[2] = y_len; dim[3] = z_len; dim[4] = zb_len;

    stat = nc_inq_varid(ncid, "t", &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_varid(ncid, "x", &x_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_varid(ncid, "y", &y_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_varid(ncid, "z", &z_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_varid(ncid, "zb", &zb_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    
    size_t t_end = t_len - 1;
    // get extent of grid by looking for last and first of variables x,y;
    // for z get first of zb and last of z
    size_t x_bdy[] = {0, x_len - 1};
    size_t y_bdy[] = {0, y_len - 1};
    size_t z_bdy[] = {0, z_len - 1}; // Start at 0 in `zb', end at z_len - 1 of `z'

    stat = nc_get_var1_double(ncid, t_id, &t_end, &bdy[0]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_var1_double(ncid, x_id, &x_bdy[0], &bdy[1]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_var1_double(ncid, x_id, &x_bdy[1], &bdy[2]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_var1_double(ncid, y_id, &y_bdy[0], &bdy[3]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_var1_double(ncid, y_id, &y_bdy[1], &bdy[4]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_var1_double(ncid, zb_id, &z_bdy[0], &bdy[5]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_var1_double(ncid, z_id, &z_bdy[1], &bdy[6]);  CHKERRQ(check_err(stat,__LINE__,__FILE__));    
  }
  MPI_Bcast(dim, 5, MPI_LONG, 0, com);
  MPI_Bcast(bdy, 7, MPI_DOUBLE, 0, com);

  return 0;
}


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
    // we index 0:M in the x direction and 0:N in the y direction.  Such a
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


//! Read in the variables \c z and \c zb from the NetCDF file; don't assume they are equally-spaced.
PetscErrorCode NCTool::get_vertical_dims(int ncid, int z_len, int zb_len, 
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


//! Read from a variable in a NetCDF file into a \c DA -managed local \c Vec \c v; a global \c Vec \c g is used for storage.
PetscErrorCode NCTool::get_local_var(const IceGrid *grid, int ncid, const char *name, nc_type type,
                                     DA da, Vec v, Vec g, const int *s, const int *c,
                                     int dims, void *a_mpi, int a_size) {

  PetscErrorCode ierr;
  ierr = get_global_var(grid, ncid, name, type, da, g, s, c, dims, a_mpi, a_size); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, g, INSERT_VALUES, v); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, g, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


//! Read from a variable in a NetCDF file into a \c DA -managed global \c Vec \c g.
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


PetscErrorCode NCTool::var_to_da_vec(IceGrid &grid, int ncid, int vid, DA da, Vec vecl,
                                     Vec vecg, Vec vindzero) {
  PetscErrorCode  ierr;
  MaskInterp defaultmasktool;
  
  defaultmasktool.number_allowed = 4;
  //FIXME: these need to come from IceModel, not be hard coded
  defaultmasktool.allowed_levels[0] = 1;  // MASK_SHEET
  defaultmasktool.allowed_levels[1] = 2;  // MASK_DRAGGING;
  defaultmasktool.allowed_levels[2] = 3;  // MASK_FLOATING;
  defaultmasktool.allowed_levels[3] = 7;  // MASK_FLOATING_OCEAN0;
  ierr = var_to_da_vec(grid,ncid,vid,da,vecl,vecg,vindzero,defaultmasktool); CHKERRQ(ierr);
  return 0;
}


//! Read a 2D variable from a NetCDF file into a \c DA -managed local \c Vec.  Only used in bootstrapping.
/*!
This procedure should probably be DEPRECATED.  It duplicates much of the function of get_local_var() and
regrid_local_var().  It uses the VecScatter created by IceModel::getIndZero(), while the other NCTool routines
directly use MPI calls for parallel communication.

This procedure is only called by IceModel::bootstrapFromFile_netCDF() and IceModel::readShelfStreamBCFromFile_netCDF().

A special ability appears here, namely the ability to take floating point values on one grid, for what is actually an
integer mask (in the NetCDF file), and then produce an integer mask, with certain allowed values, in the target 
\c DA -managed \c Vec.  This ability is not available in get_local_var().
 */
PetscErrorCode NCTool::var_to_da_vec(IceGrid &grid, int ncid, int vid, DA da, Vec vecl,
                                     Vec vecg, Vec vindzero, MaskInterp masktool) {
  PetscErrorCode  ierr;
  int stat;
  PetscScalar **ind;
  float       *f = NULL;
  int         *g = NULL;

  if (masktool.number_allowed < 1) {
    SETERRQ(99, "var_to_da_vec: number of allowed levels in masktool must be at least one");
  }
  if (grid.rank == 0) {
    int ndims, natts, dimids[NC_MAX_VAR_DIMS];
    nc_type xtype;
    char name[NC_MAX_NAME+1];
    stat = nc_inq_var(ncid, vid, name, &xtype, &ndims, dimids, &natts); CHKERRQ(nc_check(stat));
    if (ndims != 2) {
      SETERRQ2(1, "var_to_da_vec: number of dimensions = %d for %s\n",ndims, name);
    }
    // in netCDF file, we index 0:M in the x direction and 0:N in the y direction
    // location (i,j) is addressed as [i*N + j]
    size_t M, N;
    stat = nc_inq_dimlen(ncid, dimids[0], &M); CHKERRQ(nc_check(stat));
    stat = nc_inq_dimlen(ncid, dimids[1], &N); CHKERRQ(nc_check(stat));
    
    switch (xtype) {
      case NC_FLOAT:
        f = new float[M*N];
        stat = nc_get_var_float(ncid, vid, f); CHKERRQ(nc_check(stat));
        break;
      case NC_INT:
        g = new int[M*N];
        stat = nc_get_var_int(ncid, vid, g); CHKERRQ(nc_check(stat));
        break;
      default:
        SETERRQ1(1, "var_to_da_vec: NC_VAR `%s' not of type NC_INT or NC_FLOAT.\n", name);
    }

    ierr = VecGetArray2d(vindzero, grid.Mx, grid.My, 0, 0, &ind); CHKERRQ(ierr);
    
    // netCDF concepts of $\Delta x$ and $\Delta y$
    // We have rescaled the grid early in bootstrapFromFile_netCDF() to match the 
    // physical extent of the netCDF file.
    const float ncdx = 2 * grid.Lx / (M - 1);
    const float ncdy = 2 * grid.Ly / (N - 1);

    for (PetscInt i=0; i < grid.Mx; i++) {
      for (PetscInt j=0; j < grid.My; j++) {
        const float x = grid.dx * (i - grid.Mx/2);
        const float y = grid.dy * (j - grid.My/2);
        if (PetscAbs(x) > grid.Lx) {
          SETERRQ1(2, "var_to_da_vec: x=%f not in bounds.  Grid corrupted.\n", x);
        }
        if (PetscAbs(y) > grid.Ly) {
          SETERRQ1(3, "var_to_da_vec: y=%f not in bounds.  Grid corrupted.\n", y);
        }
            
        const float ii = M / 2 + x / ncdx;
        const float jj = N / 2 + y / ncdy;
        // These live in [0,1)
        const float xx = ii - floor(ii);
        const float yy = jj - floor(jj);
        // Define weights for bilinear interpolation
        const float w11 = (1 - xx) * (1 - yy);
        const float w12 = xx * (1 - yy);
        const float w21 = (1 - xx) * yy;
        const float w22 = xx * yy;
        // Locations to sample from.  These should be on the grid.
        const int i1 = int(floor(ii)) % M;
        const int i2 = int(ceil(ii)) % M;
        const int j1 = int(floor(jj)) % N;
        const int j2 = int(ceil(jj)) % N;

        PetscScalar val;
        if (g != NULL) { // an integer array
          val = w11 * g[i1*N + j1] + w21 * g[i1*N + j2]
                 + w12 * g[i2*N + j1] + w22 * g[i2*N + j2];
          if ((masktool.number_allowed == 1) || (val <= (float)masktool.allowed_levels[0])) {
            val = (float)masktool.allowed_levels[0];
          } else {
            int k=1;        
            while (k < masktool.number_allowed) {
              float mid = ( (float)masktool.allowed_levels[k-1] 
                            + (float)masktool.allowed_levels[k] ) / 2.0;
              if (val < mid) {
                val = (float)masktool.allowed_levels[k-1];
                break;
              }
              k++;
            }
            if (k >= masktool.number_allowed) {
              val = (float)masktool.allowed_levels[k-1];
            }
          }
        } else if (f != NULL) { // a float array
          val = w11 * f[i1*N + j1] + w21 * f[i1*N + j2]
                 + w12 * f[i2*N + j1] + w22 * f[i2*N + j2];
// this way in rev 151 and earlier:
//          val = w11 * f[i1*N + j1] + w12 * f[i1*N + j2]
//                 + w21 * f[i2*N + j1] + w22 * f[i2*N + j2];
        } else {
          SETERRQ(4, "var_to_da_vec: this should not happen");
        }
        
        // The backward indexing is merely to make the plots look upright with
        // the default plotting methods.  When I can make the axes work in a
        // sane manner, this can be improved.
        ierr = VecSetValue(vecg, (PetscInt) ind[grid.Mx - 1 - i][j],
                           val, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
   
    if (f != NULL) delete [] f;
    if (g != NULL) delete [] g;
    ierr = VecRestoreArray2d(vindzero, grid.Mx, grid.My, 0, 0, &ind); CHKERRQ(ierr);
  }

  ierr = VecAssemblyBegin(vecg); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vecg); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, vecg, INSERT_VALUES, vecl); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, vecg, INSERT_VALUES, vecl); CHKERRQ(ierr);
  return 0;
}


//! The "local interpolation context" holds various parameters describing the source NetCDF file for regridding.
/*! 
This procedure merely puts various information into a struct; it doesn't extract new information or do communication.
 */
PetscErrorCode NCTool::form_LocalInterpCtx(int ncid, const size_t dim[], const double bdy[],
                                           const double zlevsIN[], const double zblevsIN[],
                                           LocalInterpCtx &lic, IceGrid &grid) {
  PetscErrorCode ierr;
  const double Lx = grid.Lx,
               Ly = grid.Ly, 
               Lz = grid.Lz, 
               Lbz = grid.Lbz,
               dx = grid.dx,
               dy = grid.dy;
  
  double xbdy[2] = {-Lx + dx * grid.xs, -Lx + dx * (grid.xs + grid.xm - 1)};
  double ybdy[2] = {-Ly + dy * grid.ys, -Ly + dy * (grid.ys + grid.ym - 1)};
  double zbdy[2] = {-Lbz, Lz};

  if (bdy[1] > -Lx || bdy[2] < Lx || bdy[3] > -Ly || bdy[4] < Ly
      || -bdy[5] < Lbz || bdy[6] < Lz)
    SETERRQ(1, "Grid not a subset");

  lic.ncid = ncid;
  
  lic.delta[0] = NAN; // Delta probably will never make sense in the time dimension.
  lic.delta[1] = (bdy[2] - bdy[1]) / (dim[1] - 1);
  lic.delta[2] = (bdy[4] - bdy[3]) / (dim[2] - 1);
  lic.delta[3] = bdy[6] / (dim[3] - 1);  // corresponds to grid.dzEQ and grid.dzbEQ; IGNOR it in general!
  
  lic.start[0] = dim[0] - 1; // We use the latest time
  lic.start[1] = (int)floor((xbdy[0] - bdy[1]) / lic.delta[1]);
  lic.start[2] = (int)floor((ybdy[0] - bdy[3]) / lic.delta[2]);
  lic.start[3] = 0; // We start at the bed.
// replace:  lic.start[4] = regridFile.kbBelowHeight(zbdy[0]), so to speak
  lic.start[4] = (int)floor((zbdy[0] - bdy[5]) / lic.delta[3]);
  
  lic.fstart[0] = bdy[0];
  lic.fstart[1] = bdy[1] + lic.start[1] * lic.delta[1];
  lic.fstart[2] = bdy[3] + lic.start[2] * lic.delta[2];

  lic.count[0] = 1; // Only take one time.
  lic.count[1] = (int)ceil((xbdy[1] - lic.fstart[1]) / lic.delta[1] + 1);
  lic.count[2] = (int)ceil((ybdy[1] - lic.fstart[2]) / lic.delta[2] + 1);
// replace: lic.count[3] = regridFile.kBelowHeight(zbdy[1]) - lic.start[3] + 1;
  lic.count[3] = (int)ceil(Lz / lic.delta[3] + 1);
  lic.count[4] = dim[4] - lic.start[4];

  ierr = PetscMalloc(dim[3] * sizeof(double), &(lic.zlevs)); CHKERRQ(ierr);
  for (size_t k = 0; k < dim[3]; k++) {
    lic.zlevs[k] = zlevsIN[k];
  }

  ierr = PetscMalloc(dim[4] * sizeof(double), &(lic.zblevs)); CHKERRQ(ierr);
  for (size_t k = 0; k < dim[4]; k++) {
    lic.zblevs[k] = zblevsIN[k];
  }

  int a_len = lic.a_len = lic.count[1] * lic.count[2] * lic.count[3];
  MPI_Reduce(&a_len, &(lic.a_len), 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(lic.a_len * sizeof(float), &(lic.a)); CHKERRQ(ierr);

  return 0;
}


//! Find a 2D or 3D variable in a NetCDF file and regrid it onto the current grid; a global \c Vec is used for storage.
/*!
Simply calls regrid_global_var() after transferring the local \c Vec to a global \c Vec.
 */
PetscErrorCode NCTool::regrid_local_var(const char *vars, char c, const char *name,
                                        int dim_flag, LocalInterpCtx &lic,
                                        IceGrid &grid, DA da, Vec vec, Vec g) {
  PetscErrorCode ierr;
  ierr = regrid_global_var(vars, c, name, dim_flag, lic, grid, da, g); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, g, INSERT_VALUES, vec); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, g, INSERT_VALUES, vec); CHKERRQ(ierr);
  return 0;
}


//! Find a 2D or 3D variable in a NetCDF file and regrid it onto the current grid.
/*! 
We need to move a 2D or 3D variable from within a NetCDF file, with its "source" grid, to the
current grid, the "target" grid.  The target grid is spread across all processors.  The variable 
on the source grid is only read from processor zero.

The source grid may be coarser or finer than the target 
grid, and it may even be coarser in one dimension and finer in another.  We do require,
however, that the source grid have greater extent than the target.  That is, the values
of the \c IceParam parameters \c Lx, \c Ly, \c Lz, and \c Lbz must exceed those of the 
target grid.

We interpolate as needed for functions of two or three variables.  For two variables we take values
\f$f_{ij}, f_{i+1,j}, f_{i,j+1}, f_{i+1,j+1}\f$


Note that the procedure checks whether the single character flag \c c is in the string
\c vars.  Note that \c dim_flag is 2 for 2-D quantities, 3 for 3-D ice quantities, and 4 for 3-D 
bedrock quantities.
 */
PetscErrorCode NCTool::regrid_global_var(const char *vars, char c, const char *name,
                                         int dim_flag, LocalInterpCtx &lic,
                                         IceGrid &grid, DA da, Vec g) {
  PetscErrorCode ierr;
                                
  if (!grid.isEqualVertSpacing()) {
    SETERRQ(604,"only implemented for grids with equal dz spacing in vertical\n");
  }

  if (! strchr(vars, c)) {
    return 0;
  }
  ierr = verbPrintf(2, grid.com, "\n   %c: regridding `%s' ... ", c, name); CHKERRQ(ierr);
  
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
      break;
    default:
      SETERRQ(1, "Invalid value for `dim_flag'.");
  }
        
  for (int i = 0; i < 4; i++) sc[i] = lic.start[i];
  for (int i = 0; i < 4; i++) sc[4 + i] = lic.count[i];

  // At this point, sc[] is set up correctly for ice 3-D quantities.
  if (dim_flag == 2) { // 2-D quantity
    sc[3] = 0; sc[7] = 1;
  } else if (dim_flag == 4) { // Bedrock quantity
    sc[3] = lic.start[4]; sc[7] = lic.count[4];
  }
  
  if (grid.rank == 0) {
    int sc0[sc_len];
    for (int i = 0; i < sc_len; i++) sc0[i] = sc[i];
    for (int proc = grid.size - 1; proc >= 0; proc--) {
      if (proc == 0) {
        for (int i = 0; i < sc_len; i++) sc[i] = sc0[i];
      } else {
        MPI_Recv(sc, sc_len, MPI_INT, proc, req_tag, grid.com, &mpi_stat);
      }
      
      size_t sc_nc[sc_len];
      for (int i = 0; i < sc_len; i++) sc_nc[i] = (size_t)sc[i]; // we need size_t
      int var_id;
      stat = nc_inq_varid(lic.ncid, name, &var_id);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      stat = nc_get_vara_float(lic.ncid, var_id, &sc_nc[0], &sc_nc[4], lic.a);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      /* {
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
      } */

      int a_len = 1;
      for (int i = 0; i < dims; i++) a_len *= sc[4 + i];

      if (proc != 0) {
        MPI_Send(lic.a, a_len, MPI_FLOAT, proc, var_tag, grid.com);
      }
    }
  } else {
    MPI_Send(sc, sc_len, MPI_INT, 0, req_tag, grid.com);
    MPI_Recv(lic.a, lic.a_len, MPI_FLOAT, 0, var_tag, grid.com, &mpi_stat);
  }
  
  PetscScalar *vec_a;
  ierr = VecGetArray(g, &vec_a); CHKERRQ(ierr);

  const int ycount = lic.count[2];
  for (int i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (int j = grid.ys; j < grid.ys + grid.ym; j++) {
      int myMz, zcount;
      float bottom = 0.0;
      float zfstart = 0.0;
      if (dim_flag == 2) {
        myMz = 1;
        zcount = 1;
      } else if (dim_flag == 3) {
        myMz = grid.Mz;
        zcount = lic.count[3];
      } else if (dim_flag == 4) {
        myMz = grid.Mbz;
// replace: remove bottom
        bottom = -grid.Lbz;
        zcount = lic.count[4];
// replace: zfstart = regridFile.zblevel[(regridFile.Mbz - 1) - (zcount - 1)], so to speak
        zfstart = -(zcount - 1) * lic.delta[3];
      }

      for (int k = 0; k < myMz; k++) {
        float a_mm, a_mp, a_pm, a_pp;
        
        const float x = -grid.Lx + i * grid.dx;
        const float y = -grid.Ly + j * grid.dy;
// replace:  const float z = (dim_flag == 4) ? grid.zlevels[k] : grid.zblevels[k];
        const float z = k * grid.dzMIN + bottom;  // assumes grid.dzMIN = grid.dzMAX, for now

        const float ic = (x - lic.fstart[1]) / lic.delta[1];
        const float jc = (y - lic.fstart[2]) / lic.delta[2];
        if (dim_flag == 3 || dim_flag == 4) {
// replace: ???
          const float kc = (z - zfstart) / lic.delta[3];

          int mmm = ((int)floor(ic) * ycount + (int)floor(jc)) * zcount + (int)floor(kc);
          int mmp = ((int)floor(ic) * ycount + (int)floor(jc)) * zcount + (int)ceil(kc);
          int mpm = ((int)floor(ic) * ycount + (int)ceil(jc)) * zcount + (int)floor(kc);
          int mpp = ((int)floor(ic) * ycount + (int)ceil(jc)) * zcount + (int)ceil(kc);
          int pmm = ((int)ceil(ic) * ycount + (int)floor(jc)) * zcount + (int)floor(kc);
          int pmp = ((int)ceil(ic) * ycount + (int)floor(jc)) * zcount + (int)ceil(kc);
          int ppm = ((int)ceil(ic) * ycount + (int)ceil(jc)) * zcount + (int)floor(kc);
          int ppp = ((int)ceil(ic) * ycount + (int)ceil(jc)) * zcount + (int)ceil(kc);

          const float kk = kc - floor(kc);
          a_mm = lic.a[mmm] * (1.0 - kk) + lic.a[mmp] * kk;
          a_mp = lic.a[mpm] * (1.0 - kk) + lic.a[mpp] * kk;
          a_pm = lic.a[pmm] * (1.0 - kk) + lic.a[pmp] * kk;
          a_pp = lic.a[ppm] * (1.0 - kk) + lic.a[ppp] * kk;

          // printf("indices = %d %d %d %d\n", mmp, mpp, pmp, ppp);
          // printf("values  = %f %f %f %f\n", a_mm, a_mp, a_pm, a_pp);
        } else {
          a_mm = lic.a[(int)floor(ic) * ycount + (int)floor(jc)];
          a_mp = lic.a[(int)floor(ic) * ycount + (int)ceil(jc)];
          a_pm = lic.a[(int)ceil(ic) * ycount + (int)floor(jc)];
          a_pp = lic.a[(int)ceil(ic) * ycount + (int)ceil(jc)];
        }

        const float jj = jc - floor(jc);
        const float a_m = a_mm * (1.0 - jj) + a_mp * jj;
        const float a_p = a_pm * (1.0 - jj) + a_pp * jj;
        
        const float ii = ic - floor(ic);
        int index = ((i - grid.xs) * grid.ym + (j - grid.ys)) * myMz + k;
        vec_a[index] = a_m * (1.0 - ii) + a_p * ii;
      }
    }
  }
  
  ierr = VecRestoreArray(g, &vec_a); CHKERRQ(ierr);

  return 0;
}


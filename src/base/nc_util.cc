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

NCTool::NCTool(IceGrid *my_grid) {
  //FIXME: does there need to be a default MaskInterp?
  myMaskInterp = PETSC_NULL;
  grid = my_grid;
  ncid = -1;

  // Initialize UDUNITS if needed
  if (utIsInit() == 0) {
    if (utInit(NULL) != 0) {
      PetscPrintf(grid->com, "PISM ERROR: UDUNITS initialization failed.\n");
      PetscEnd();
    }
  }
}

PetscErrorCode  NCTool::find_dimension(const char short_name[], int *dimid, bool &exists) {
  PetscErrorCode ierr;
  int stat, found = 0, my_dimid;
  if (grid->rank == 0) {
      stat = nc_inq_dimid(ncid, short_name, &my_dimid);
      if (stat == NC_NOERR)
	found = 1;
  }
  ierr = MPI_Bcast(&found, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&my_dimid, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);


  if (found) {
    exists = true;
    if (dimid != NULL)
      *dimid = my_dimid;
  } else {
    exists = false;
    // dimid is not modified
  }

  return 0;
}

//! Finds the variable by its standard_name attribute.
/*!
  Here's how it works:
  <ol>
  <li> Check if the current IceModelVec has a standard_name. If it does, go to
  step 2, otherwise go to step 4.

  <li> Find the variable with this standard_name. Bail out if two
  variables have the same standard_name, otherwise go to step 3.

  <li> If no variable was found, go to step 4, otherwise go to step 5.

  <li> Find the variable with the right variable name. Go to step 5.

  <li> Broadcast the existence flag and the variable ID.
  </ol>
 */
PetscErrorCode  NCTool::find_variable(string short_name, string standard_name,
				      int *varidp, bool &exists) {
  int ierr;
  int stat, found = 0, my_varid = -1, nvars;

  if (standard_name != "") {
    if (grid->rank == 0) {
      ierr = nc_inq_nvars(ncid, &nvars); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }
    ierr = MPI_Bcast(&nvars, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);

    string attribute;
    for (int j = 0; j < nvars; j++) {
      ierr = get_att_text(j, "standard_name", attribute); CHKERRQ(ierr);
      if (attribute == "")
	continue;

      if (attribute == standard_name) {
	if (!found) {
	  found = true;
	  my_varid = j;
	} else {
	  ierr = PetscPrintf(grid->com,
			     "PISM ERROR: Inconsistency in the input file: "
			     "Variables #%d and #%d have the same standard_name ('%s').\n",
			     my_varid, j, attribute.c_str());
	  CHKERRQ(ierr);
	  PetscEnd();
	}
      }
    } // end of for (int j = 0; j < nvars; j++)
  } // end of if (standard_name != NULL)

  // Check the short name:
  if (!found) {
    if (grid->rank == 0) {
      stat = nc_inq_varid(ncid, short_name.c_str(), &my_varid);
      if (stat == NC_NOERR)
	found = true;
    }
    ierr = MPI_Bcast(&found,    1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
    ierr = MPI_Bcast(&my_varid, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  }

  if (found) {
    exists = true;
    if (varidp != NULL)
      *varidp = my_varid;
  } else {
    exists = false;
    // *varidp is not modified
  }

  return 0;
}

PetscErrorCode NCTool::find_variable(string short_name, int *varid, bool &exists) {
  return find_variable(short_name, "", varid, exists);
}

//! Read the first and last values, and the lengths, of the x,y,z,zb dimensions from a NetCDF file.  Read the last t.
PetscErrorCode NCTool::get_grid_info(grid_info &g) {
  PetscErrorCode ierr;

  ierr = get_grid_info_2d(g); CHKERRQ(ierr);

  ierr = get_dim_length("z",  &g.z_len);  CHKERRQ(ierr);
  ierr = get_dim_length("zb", &g.zb_len); CHKERRQ(ierr);
  ierr = get_dim_limits("zb", &g.zb_min, NULL);     CHKERRQ(ierr);
  ierr = get_dim_limits("z",  NULL,      &g.z_max); CHKERRQ(ierr);
  
  return 0;
}

//! Read the first and last values, and the lengths, of the x,y dimensions from a NetCDF file.  Read the last t.
PetscErrorCode NCTool::get_grid_info_2d(grid_info &g) {
  PetscErrorCode ierr;

  ierr = get_dim_length("t",  &g.t_len); CHKERRQ(ierr);
  ierr = get_dim_length("y",  &g.x_len); CHKERRQ(ierr); // transpose
  ierr = get_dim_length("x",  &g.y_len); CHKERRQ(ierr); // transpose

  ierr = get_dim_limits("t", NULL, &g.time); CHKERRQ(ierr);
  ierr = get_dim_limits("y", &g.x_min, &g.x_max); CHKERRQ(ierr); // transpose
  ierr = get_dim_limits("x", &g.y_min, &g.y_max); CHKERRQ(ierr); // transpose

  g.x0 = (g.x_max + g.x_min) / 2.0;
  g.y0 = (g.y_max + g.y_min) / 2.0;

  g.Lx = (g.x_max - g.x_min) / 2.0;
  g.Ly = (g.y_max - g.y_min) / 2.0;

  return 0;
}

//! Read the last value of the time variable t from a NetCDF file.
PetscErrorCode NCTool::get_last_time(double *time) {
  PetscErrorCode ierr;
  ierr = get_dim_limits("t", NULL, time); CHKERRQ(ierr);
  return 0;
}


//! Read in the variables \c z and \c zb from the NetCDF file; <i>do not</i> assume they are equally-spaced.
/*!
  This function allocates arrays z_levels and zb_levels, and they have to be
  freed by the caller (using delete[]).
 */
PetscErrorCode NCTool::get_vertical_dims(double* &z_levels, double* &zb_levels) {
  int stat;
  int z_id, zb_id, z_len, zb_len;
  size_t zero  = 0, nc_z_len, nc_zb_len;

  stat = get_dim_length("z", &z_len); CHKERRQ(stat);
  stat = get_dim_length("zb", &zb_len); CHKERRQ(stat);

  z_levels = new double[z_len];
  zb_levels = new double[zb_len];

  nc_z_len  = (size_t) z_len;
  nc_zb_len = (size_t) zb_len;

  if (grid->rank == 0) {
    stat = nc_inq_varid(ncid, "z", &z_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_varid(ncid, "zb", &zb_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_get_vara_double(ncid, z_id, &zero, &nc_z_len, z_levels);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_vara_double(ncid, zb_id, &zero, &nc_zb_len, zb_levels);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  MPI_Bcast(z_levels, z_len, MPI_DOUBLE, 0, grid->com);
  MPI_Bcast(zb_levels, zb_len, MPI_DOUBLE, 0, grid->com);
  return 0;
}


//! Put the variable for a dimension in a NetCDF file.
/*! Uses starting and ending values and a grid length for regularly-spaced
values.
 */
PetscErrorCode NCTool::put_dimension_regular(int varid, int len, double start, double end) {
  PetscErrorCode ierr;
  int stat;
  double *v, delta;

  if (end <= start)
    SETERRQ2(1, "Can't write dimension: start = %f, end = %f",
	     start, end);

  delta = (end - start) / (len - 1);

  ierr = PetscMalloc(len * sizeof(double), &v); CHKERRQ(ierr);
  for (int i = 0; i < len; i++)
    v[i] = start + i * delta;

  // Sometimes v[len - 1] turns out to be greater than end (because of rounding
  // errors). If that happens, we need to fix v[len - 1] by setting it equal to
  // end.
  if (v[len - 1] > end) v[len - 1] = end;
  
  stat = nc_put_var_double(ncid, varid, v); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  ierr = PetscFree(v); CHKERRQ(ierr);

  return 0;
}


//! Put the variable for a dimension in a NetCDF file.  Makes no assumption about spacing.
PetscErrorCode NCTool::put_dimension(int varid, int len, PetscScalar *vals) {
  PetscErrorCode ierr;
  int stat;
  double *v;

  ierr = PetscMalloc(len * sizeof(double), &v); CHKERRQ(ierr);
  for (int i = 0; i < len; i++) {
    v[i] = (double)vals[i];
  }
  stat = nc_put_var_double(ncid, varid, v); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  ierr = PetscFree(v); CHKERRQ(ierr);
  return 0;
}


//! Read a NetCDF variable into a \c DA -managed local \c Vec \c v; a global \c Vec \c g is used for storage.
/*!
Just calls get_global_var().  Then transfers the global \c Vec \c g to the local \c Vec \c vec.
 */
PetscErrorCode NCTool::get_local_var(const int varid, DA da, Vec v, GridType dims, int t) {

  PetscErrorCode ierr;
  Vec g;
  ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);

  ierr = get_global_var(varid, g, dims, t); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, g, INSERT_VALUES, v); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, g, INSERT_VALUES, v); CHKERRQ(ierr);

  ierr = VecDestroy(g); CHKERRQ(ierr);
  return 0;
}


//! Read a variable in a NetCDF file into a \c DA -managed global \c Vec \c g.  \e In \e parallel.
PetscErrorCode NCTool::get_global_var(const int varid, Vec g, GridType dims, int t) {
  const int N = 5;
  const int start_tag = 1;
  const int count_tag = 2;
  const int data_tag =  3;
  PetscErrorCode ierr;
  MPI_Status mpi_stat;
  int stat;
  double *a_double = NULL, *a_raw = NULL;
  
  int start[N] = {t, grid->xs, grid->ys, 0,        0};
  int count[N] = {1, grid->xm, grid->ym, grid->Mz, grid->Mbz};

  // Find the local size:
  int block_size, buffer_size;
  block_size = compute_block_size(dims, count);
  // And the maximum size of the data block:
  buffer_size = block_size;
  ierr = MPI_Reduce(&block_size, &buffer_size, 1, MPI_INT, MPI_MAX, 0, grid->com); CHKERRQ(ierr);

  // Memory allocation:
  ierr = PetscMalloc(buffer_size * sizeof(double), &a_double); CHKERRQ(ierr);

  if (grid->rank == 0) {
    // Allocate an array to store input in (before the transpose)
    ierr = PetscMalloc(buffer_size * sizeof(double), &a_raw); CHKERRQ(ierr);

    int start0[N], count0[N];
    for (int j = 0; j < N; j++) {
      // root needs to save its range
      start0[j] = start[j];
      count0[j] = count[j];
    }
    for (int proc = grid->size - 1; proc >= 0; proc--) { // root will read itself last
      if (proc == 0) {
        for (int j = 0; j < N; j++) {
	  start[j] = start0[j];
	  count[j] = count0[j];
	}
      } else {
        MPI_Recv(start, N, MPI_INT, proc, start_tag, grid->com, &mpi_stat);
        MPI_Recv(count, N, MPI_INT, proc, count_tag, grid->com, &mpi_stat);
      }

      size_t *nc_start, *nc_count;
      ierr = compute_start_and_count(varid, start, count, nc_start, nc_count); CHKERRQ(ierr);

      stat = nc_get_vara_double(ncid, varid, nc_start, nc_count, a_raw);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      delete[] nc_start;
      delete[] nc_count;

      ierr = transpose(varid, dims, count, a_raw, a_double); CHKERRQ(ierr);

      if (proc != 0) {
        int block_size;
	block_size = compute_block_size(dims, count);
	MPI_Send(a_double, block_size, MPI_DOUBLE, proc, data_tag, grid->com);
      }
    }
    // Clean up:
    ierr = PetscFree(a_raw); CHKERRQ(ierr);
  } else {
    MPI_Send(start, N, MPI_INT, 0, start_tag, grid->com);
    MPI_Send(count, N, MPI_INT, 0, count_tag, grid->com);
    MPI_Recv(a_double, buffer_size, MPI_DOUBLE, 0, data_tag, grid->com, &mpi_stat);
  }

  block_size = compute_block_size(dims, count);
  PetscScalar *a_petsc;
  ierr = VecGetArray(g, &a_petsc); CHKERRQ(ierr);
  for (int i = 0; i < block_size; i++) {
    a_petsc[i] = (PetscScalar)a_double[i];
  }
  ierr = VecRestoreArray(g, &a_petsc); CHKERRQ(ierr);

  // Cleanup:
  ierr = PetscFree(a_double); CHKERRQ(ierr);
  return 0;
}


//! Put a \c DA -managed local \c Vec \c v into a variable in a NetCDF file
/*!
Just calls put_global_var(), after transfering the local \c Vec called \c v into the global \c Vec \c g.
 */
PetscErrorCode NCTool::put_local_var(const int varid, DA da, Vec v, GridType dims) {

  PetscErrorCode ierr;
  Vec g;
  ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);
  ierr = DALocalToGlobal(da, v, INSERT_VALUES, g); CHKERRQ(ierr);

  ierr = put_global_var(varid, g, dims); CHKERRQ(ierr);

  ierr = VecDestroy(g); CHKERRQ(ierr);
  return 0;
}


//! Put a \c DA -managed global \c Vec \c g into a variable in a NetCDF file.  \e In \e parallel.
PetscErrorCode NCTool::put_global_var(const int varid, Vec g, GridType dims) {
  const int start_tag = 1;
  const int count_tag = 2;
  const int data_tag =  3;
  const int N = 5;
  PetscErrorCode ierr;
  MPI_Status mpi_stat;
  int stat;
  double *a_double = NULL;

  // Fill start and count arrays:
  int t;
  ierr = get_dim_length("t", &t); CHKERRQ(ierr);
  int start[N] = {t - 1, grid->xs, grid->ys, 0,        0};
  int count[N] = {1,     grid->xm, grid->ym, grid->Mz, grid->Mbz};

  // Find the local size:
  int block_size, buffer_size;
  block_size = compute_block_size(dims, count);
  // And the maximum size of the data block:
  buffer_size = block_size;
  ierr = MPI_Reduce(&block_size, &buffer_size, 1, MPI_INT, MPI_MAX, 0, grid->com); CHKERRQ(ierr);

  // Memory allocation:
  ierr = PetscMalloc(buffer_size * sizeof(double), &a_double); CHKERRQ(ierr);

  // convert IceModel Vec containing PetscScalar to array of double for NetCDF
  PetscScalar *a_petsc;
  ierr = VecGetArray(g, &a_petsc); CHKERRQ(ierr);
  for (int i = 0; i < block_size; i++) {
    a_double[i] = (double)a_petsc[i];
  }
  ierr = VecRestoreArray(g, &a_petsc); CHKERRQ(ierr);

  if (grid->rank == 0) { // on rank 0 processor, receive messages from every other
                         //    processor, then write it out to the NC file
    for (int proc = 0; proc < grid->size; proc++) { // root will write itself last
      if (proc != 0) {
        MPI_Recv(start, N, MPI_INT, proc, start_tag, grid->com, &mpi_stat);
        MPI_Recv(count, N, MPI_INT, proc, count_tag, grid->com, &mpi_stat);
	MPI_Recv(a_double, buffer_size, MPI_DOUBLE, proc, data_tag, grid->com, &mpi_stat);
      }
      
      size_t *nc_start, *nc_count;
      ierr = compute_start_and_count(varid, start, count, nc_start, nc_count); CHKERRQ(ierr);

      stat = nc_put_vara_double(ncid, varid, nc_start, nc_count, a_double);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      delete[] nc_start;
      delete[] nc_count;
    }
  } else {  // all other processors send to rank 0 processor
    MPI_Send(start, N, MPI_INT, 0, start_tag, grid->com);
    MPI_Send(count, N, MPI_INT, 0, count_tag, grid->com);
    MPI_Send(a_double, buffer_size, MPI_DOUBLE, 0, data_tag, grid->com);
  }

  // Cleanup:
  ierr = PetscFree(a_double); CHKERRQ(ierr);
  return 0;
}


//! Call this before calling regrid_local_var() or regrid_global_var() with argument useMaskInterp = true.
PetscErrorCode NCTool::set_MaskInterp(MaskInterp *mi_in) {
  myMaskInterp = mi_in;
  return 0;
}


//! Regrid a 2D or 3D NetCDF variable using a global \c Vec for storage.
/*!
Simply calls regrid_global_var().  Then transfers the global \c Vec \c g to the local \c Vec \c vec.
 */
PetscErrorCode NCTool::regrid_local_var(const int varid, GridType dim_flag, LocalInterpCtx &lic,
                                        DA da, Vec vec, 
                                        bool useMaskInterp) {
  PetscErrorCode ierr;
  Vec g;
  ierr = DACreateGlobalVector(da, &g); CHKERRQ(ierr);

  ierr = regrid_global_var(varid, dim_flag, lic, g, useMaskInterp); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, g, INSERT_VALUES, vec); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, g, INSERT_VALUES, vec); CHKERRQ(ierr);

  ierr = VecDestroy(g); CHKERRQ(ierr);
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
PetscErrorCode NCTool::regrid_global_var(const int varid, GridType dim_flag, LocalInterpCtx &lic,
                                         Vec g, bool useMaskInterp) {
  PetscErrorCode ierr;
  const int N = 5, X = 1, Y = 2, Z = 3, ZB = 4; // indices, just for clarity
  const int start_tag = 1; // MPI tag for the start array
  const int count_tag = 2; // MPI tag for the count array
  const int data_tag  = 3; // MPI tag for the data block
  MPI_Status mpi_stat;
  int stat, start[N], count[N];	// enough space for t, x, y, z, zb

  // make local copies of lic.start and lic.count
  for (int i = 0; i < N; i++) {
    start[i] = lic.start[i];
    count[i] = lic.count[i];
  }

  // This ensures that we don't read too many values (i.e. if varid refers to a
  // 3D variable and 2D-only regridding is requested) and that the size of the
  // buffer (a_len) is calculated correctly.
  switch (dim_flag) {
  case GRID_2D:
    count[Z] = 1;
    count[ZB] = 1;
    break;
  case GRID_3D:
    if (lic.regrid_2d_only) {
      SETERRQ(2, "regrid_2d_only is set, so dim_flag == GRID_3D is not allowed");
    }
    count[ZB] = 1;
    break;
  case GRID_3D_BEDROCK:
    if (lic.no_regrid_bedrock) {
      SETERRQ(2, "no_regrid_bedrock is set, so dim_flag == GRID_3D_BEDROCK is not allowed");
    }
    count[Z] = 1;
    break;
  }

  if (grid->rank == 0) {

    // Node 0 will service all the other nodes before itself.  We need to save
    // start[] and count[] so that it knows how to get its block at the end.
    int start0[N];
    int count0[N];

    for (int i = 0; i < N; i++) {
      start0[i] = start[i];
      count0[i] = count[i];
    }

    for (int proc = grid->size - 1; proc >= 0; proc--) {
      if (proc == 0) {// Get the bounds.
	for (int i = 0; i < N; i++) {
	  start[i] = start0[i];
	  count[i] = count0[i];
	}
      } else {
        MPI_Recv(start, N, MPI_INT, proc, start_tag, grid->com, &mpi_stat);
	MPI_Recv(count, N, MPI_INT, proc, count_tag, grid->com, &mpi_stat);
      }

      // Assemble nc_start and nc_count that are used below in the call to
      // nc_get_vara_double(...)
      size_t *nc_start, *nc_count;
      ierr = compute_start_and_count(varid, start, count, nc_start, nc_count); CHKERRQ(ierr);

      stat = nc_get_vara_double(ncid, varid, nc_start, nc_count, lic.a_raw);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      delete[] nc_start;
      delete[] nc_count;

      ierr = transpose(varid, dim_flag, count, lic.a_raw, lic.a); CHKERRQ(ierr);

      // Find out how big the buffer actually is.  Remember that node 0 has a
      // buffer that will only be filled by if the process it is serving has a
      // maximal sized local domain.
      // Also note that at least one of count[Z] and count[ZB] is 1.
      int a_len = count[X] * count[Y] * count[Z] * count[ZB];

      // send the filled buffer
      if (proc != 0) {
        MPI_Send(lic.a, a_len, MPI_DOUBLE, proc, data_tag, grid->com);
      }
    }
  } else { // not process 0:
    MPI_Send(start, N, MPI_INT, 0, start_tag, grid->com);  // send out my start
    MPI_Send(count, N, MPI_INT, 0, count_tag, grid->com);  // send out my count
    MPI_Recv(lic.a, lic.a_len, MPI_DOUBLE, 0, data_tag, grid->com, &mpi_stat); // get back filled buffer
  }

  // At this point, the buffer lic.a[] should contain lic.a_len doubles.  This is the 
  // local processor's part of the source variable.  
  // That is, it should be enough of the source variable so that \e interpolation
  // (not extrapolation) onto the local processor's part of the target grid is possible.
  //ierr = lic.printArray(grid.com); CHKERRQ(ierr);
  
  // indexing parameters
  int myMz = 1, zcount = 1; // indexing trick so that we don't have to duplicate code for the 2-D case.
  switch (dim_flag) {
  case GRID_3D:
    myMz = grid->Mz;
    zcount = count[Z];
    break;
  case GRID_3D_BEDROCK:
    myMz = grid->Mbz;
    zcount = count[ZB];
    break;
  default:
    break;
  }

  // We'll work with the raw storage here so that the array we are filling is
  // indexed the same way as the buffer we are pulling from (lic.a)
  PetscScalar *vec_a;
  ierr = VecGetArray(g, &vec_a); CHKERRQ(ierr);

  for (int i = grid->xs; i < grid->xs + grid->xm; i++) {
    for (int j = grid->ys; j < grid->ys + grid->ym; j++) {

      for (int k = 0; k < myMz; k++) {
        // location (x,y,z) is in target computational domain
        const double x = -grid->Lx + i * grid->dx,
                     y = -grid->Ly + j * grid->dy,
                     z = (dim_flag == GRID_3D) ? grid->zlevels[k] : grid->zblevels[k];

        // We need to know how the point (x,y,z) sits within the local block we
        // pulled from the netCDF file.  This part is special to a regular
        // grid.  In particular floor(ic) is the index of the 'left neighbor'
        // and ceil(ic) is the index of the 'right neighbor'.
        double ic = (x - lic.fstart[X]) / lic.delta[X];
        double jc = (y - lic.fstart[Y]) / lic.delta[Y];

        // bounds checking on ic and jc; cases where this is essential *have* been observed
        if ((int)floor(ic) < 0) {
          ic = 0.0;
          //SETERRQ2(101,"(int)floor(ic) < 0      [%d < 0; ic = %16.15f]",(int)floor(ic),ic);
        }
        if ((int)floor(jc) < 0) {
          jc = 0.0;
          //SETERRQ2(102,"(int)floor(jc) < 0      [%d < 0; jc = %16.15f]",(int)floor(jc),jc);
        }
        if ((int)ceil(ic) > lic.count[X]-1) {
          ic = (double)(lic.count[X]-1);
          //SETERRQ3(103,"(int)ceil(ic) > lic.count[1]-1      [%d > %d; ic = %16.15f]",
          //     (int)ceil(ic), lic.count[1]-1, ic);
        }
        if ((int)ceil(jc) > lic.count[Y]-1) {
          jc = (double)(lic.count[Y]-1);
          //SETERRQ3(104,"(int)ceil(jc) > lic.count[2]-1      [%d > %d; jc = %16.15f]",
          //     (int)ceil(jc), lic.count[2]-1, jc);
        }

        double a_mm, a_mp, a_pm, a_pp;  // filled differently in 2d and 3d cases

        if (dim_flag == GRID_3D || dim_flag == GRID_3D_BEDROCK) {
          // get the index into the source grid, for just below the level z
          const int kc = (dim_flag == GRID_3D) ? lic.kBelowHeight(z)
                                               : lic.kbBelowHeight(z);

          // We pretend that there are always 8 neighbors.  And compute the
          // indices into the buffer for those neighbors.  
          // Note that floor(ic) + 1 = ceil(ic) does not hold when ic is an
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
          const int mmm = ((int)floor(ic) * count[Y] + (int)floor(jc)) * zcount + kc;
          const int mmp = ((int)floor(ic) * count[Y] + (int)floor(jc)) * zcount + kc + 1;
          const int mpm = ((int)floor(ic) * count[Y] + (int)ceil(jc)) * zcount + kc;
          const int mpp = ((int)floor(ic) * count[Y] + (int)ceil(jc)) * zcount + kc + 1;
          const int pmm = ((int)ceil(ic) * count[Y] + (int)floor(jc)) * zcount + kc;
          const int pmp = ((int)ceil(ic) * count[Y] + (int)floor(jc)) * zcount + kc + 1;
          const int ppm = ((int)ceil(ic) * count[Y] + (int)ceil(jc)) * zcount + kc;
          const int ppp = ((int)ceil(ic) * count[Y] + (int)ceil(jc)) * zcount + kc + 1;

          // We know how to index the neighbors, but we don't yet know where the
          // point lies within this box.  This is represented by kk in [0,1].
          // For the irregular case, with left index km and right index kp in the source grid, we
          // would have
          //   kk = (km == kp) ? 0.0 : (z - Z(km)) / (Z(kp) - Z(km))
          // where Z(.) are the physical coordinates on the source grid.  Note
          // that any value in [0,1] would be okay when km == kp.
          const double zkc = (dim_flag == GRID_3D) ? lic.zlevs[kc] : lic.zblevs[kc];
          double dz;
          if (dim_flag == GRID_3D) {
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
          a_mm = (double)(lic.a[(int)floor(ic) * count[Y] + (int)floor(jc)]);
          a_mp = (double)(lic.a[(int)floor(ic) * count[Y] + (int)ceil(jc)]);
          a_pm = (double)(lic.a[(int)ceil(ic) * count[Y] + (int)floor(jc)]);
          a_pp = (double)(lic.a[(int)ceil(ic) * count[Y] + (int)ceil(jc)]);
        }

        const double jj = jc - floor(jc);

        // interpolate in y direction
        const double a_m = a_mm * (1.0 - jj) + a_mp * jj;
        const double a_p = a_pm * (1.0 - jj) + a_pp * jj;

        const double ii = ic - floor(ic);
        int index = ((i - grid->xs) * grid->ym + (j - grid->ys)) * myMz + k;

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
            int kk;
            for (kk=1; kk < myMaskInterp->number_allowed; kk++) {
             const double mid = (  (double)(myMaskInterp->allowed_levels[kk-1])
                                 + (double)(myMaskInterp->allowed_levels[kk])   ) / 2.0;
              if (val < mid) {
                val = (double)(myMaskInterp->allowed_levels[k-1]);
                break;
              }
              kk++;
            }
            if (kk >= myMaskInterp->number_allowed) {
              val = (double)(myMaskInterp->allowed_levels[kk-1]);
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

//! Writes \c history to the history global attribute of a NetCDF dataset.
/*!
  Appends if overwrite == false (default).
 */
PetscErrorCode NCTool::write_history(const char history[], bool overwrite) {
  int stat;
  string old_history, new_history;

  // Produce the new history string:
  stat = get_att_text(NC_GLOBAL, "history", old_history); CHKERRQ(stat);
  
  if (overwrite) {
    new_history = history;
  } else {
    new_history = history + old_history;
  }

  // Write it:
  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    
    stat = nc_put_att_text(ncid, NC_GLOBAL, "history",
			   new_history.size(), new_history.c_str());
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
      
    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  return 0;
}

PetscErrorCode NCTool::read_polar_stereographic(PolarStereoParams &ps, bool report) {
  PetscErrorCode ierr;
  vector<double> tmp;
  int varid;
  bool ps_exists, lon_exists = false, lat_exists = false, par_exists = false;

  ierr = find_variable("polar_stereographic", &varid, ps_exists); CHKERRQ(ierr);
  if (!ps_exists && report) {
    ierr = verbPrintf(2,grid->com,
		      "  polar stereographic variable not found, using defaults: svlfp=%6.2f, lopo=%6.2f, sp=%6.2f\n",
		      ps.svlfp, ps.lopo, ps.sp); CHKERRQ(ierr);
    return 0;
  }

  ierr = get_att_double(varid, "straight_vertical_longitude_from_pole", tmp); CHKERRQ(ierr);
  if (tmp.size() == 1) {
    ps.svlfp = tmp[0];
    lon_exists = true;
  }

  ierr = get_att_double(varid, "latitude_of_projection_origin", tmp); CHKERRQ(ierr);
  if (tmp.size() == 1) {
    ps.lopo = tmp[0];
    lat_exists = true;
  }

  ierr = get_att_double(varid, "standard_parallel", tmp); CHKERRQ(ierr);
  if (tmp.size() == 1) {
    ps.sp = tmp[0];
    par_exists = true;
  }

  if (report) {
    ierr = verbPrintf(2, grid->com,
		      "  polar stereographic variable found; attributes present: svlfp=%d, lopo=%d, sp=%d\n"
		      "     values: svlfp = %6.2f, lopo = %6.2f, sp = %6.2f\n",
		      lon_exists, lat_exists, par_exists,
		      ps.svlfp, ps.lopo, ps.sp); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode NCTool::write_polar_stereographic(PolarStereoParams &ps) {
  int stat, varid;

  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    // Check if polar_stereographic exists and define it if it does not.
    stat = nc_inq_varid(ncid, "polar_stereographic", &varid);
    if (stat == NC_ENOTVAR) {
      stat = nc_def_var(ncid, "polar_stereographic", NC_INT, 0, 0, &varid);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));
    }

    // Write attributes.
    stat = nc_put_att_text(ncid, varid, "grid_mapping_name", 19, "polar_stereographic");
    check_err(stat,__LINE__,__FILE__);

    stat = nc_put_att_double(ncid, varid, "straight_vertical_longitude_from_pole", NC_DOUBLE, 1,
			     &ps.svlfp);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_put_att_double(ncid, varid, "latitude_of_projection_origin", NC_DOUBLE, 1,
			     &ps.lopo);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_put_att_double(ncid, varid, "standard_parallel", NC_DOUBLE, 1,
			     &ps.sp);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, varid, "pism_intent", 7, "mapping");
    check_err(stat,__LINE__,__FILE__);

    stat = nc_enddef(ncid);
  }
  return 0;
}

//! Checks if the dimension dim in a NetCDF file is OK.
/*! A dimension is OK if is exists and its length is equal to len. If len < 0,
    then dimension length is ignored.

    On processor 0 returns true if OK, false otherwise. Always returns true on
    processors other than 0.
 */
bool NCTool::check_dimension(const char name[], const int len) {
  int stat, dimid;
  size_t dimlen;

  if (grid->rank == 0) {
    stat = nc_inq_dimid(ncid, name, &dimid);
    if (stat == NC_NOERR) {
      stat = nc_inq_dimlen(ncid, dimid, &dimlen);

      if (stat != NC_NOERR)
	return false;
      if ((len < 0) || ((int)dimlen == len))
	return true;
      if ((int)dimlen != len)
	return false;
    }
  }

  return true;
}

//! Always returns true on processors other than zero.
bool NCTool::check_dimensions() {
  bool t, x, y, z, zb;

  t  = check_dimension("t", -1); // length does not matter
  x  = check_dimension("x", grid->My); // transpose
  y  = check_dimension("y", grid->Mx); // transpose
  z  = check_dimension("z", grid->Mz);
  zb = check_dimension("zb", grid->Mbz);
  
  return (t & x & y & z & zb);
}

//! Create dimensions and coordinate variables.
/*! Assumes that the dataset is in the data mode. */
PetscErrorCode NCTool::create_dimensions() {
  int stat, t, x, y, z, zb, dimid;

  if (grid->rank == 0) {
    // define dimensions and coordinate variables:
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    // t
    stat = nc_def_dim(ncid, "t", NC_UNLIMITED, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_def_var(ncid, "t", NC_DOUBLE, 1, &dimid, &t); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_put_att_text(ncid, t, "long_name", 4, "time"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, t, "units", 17, "years since 1-1-1"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, t, "calendar", 4, "none"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, t, "axis", 1, "T"); check_err(stat,__LINE__,__FILE__);
    // x; note the transpose
    stat = nc_def_dim(ncid, "x", grid->My, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_def_var(ncid, "x", NC_DOUBLE, 1, &dimid, &x); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_put_att_text(ncid, x, "axis", 1, "X"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, x, "long_name", 32, "X-coordinate in Cartesian system"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, x, "standard_name", 23, "projection_x_coordinate"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, x, "units", 1, "m"); check_err(stat,__LINE__,__FILE__);
    // y; note the transpose
    stat = nc_def_dim(ncid, "y", grid->Mx, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_def_var(ncid, "y", NC_DOUBLE, 1, &dimid, &y); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_put_att_text(ncid, y, "axis", 1, "Y"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, y, "long_name", 32, "Y-coordinate in Cartesian system"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, y, "standard_name", 23, "projection_y_coordinate"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, y, "units", 1, "m"); check_err(stat,__LINE__,__FILE__);
    // z
    stat = nc_def_dim(ncid, "z", grid->Mz, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_def_var(ncid, "z", NC_DOUBLE, 1, &dimid, &z); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_put_att_text(ncid, z, "axis", 1, "Z"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, z, "long_name", 32, "z-coordinate in Cartesian system"); check_err(stat,__LINE__,__FILE__);
    // PROPOSED standard_name = projection_z_coordinate
    stat = nc_put_att_text(ncid, z, "units", 1, "m"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, z, "positive", 2, "up"); check_err(stat,__LINE__,__FILE__);
    // zb
    stat = nc_def_dim(ncid, "zb", grid->Mbz, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_def_var(ncid, "zb", NC_DOUBLE, 1, &dimid, &zb); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_put_att_text(ncid, zb, "long_name", 23, "z-coordinate in bedrock"); check_err(stat,__LINE__,__FILE__);
    // PROPOSED standard_name = projection_z_coordinate_in_lithosphere
    stat = nc_put_att_text(ncid, zb, "units", 1, "m"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, zb, "positive", 2, "up"); check_err(stat,__LINE__,__FILE__);
    // 
    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    // set values:
    // Note that the 't' dimension is not modified: it is handled by the append_time method.
    // Also note the transpose.
    stat = put_dimension_regular(x, grid->My,
				 grid->y0 - grid->Ly,
				 grid->y0 + grid->Ly); CHKERRQ(stat);
    stat = put_dimension_regular(y, grid->Mx,
				 grid->x0 - grid->Lx,
				 grid->x0 + grid->Lx); CHKERRQ(stat);
    stat = put_dimension(z, grid->Mz, grid->zlevels); CHKERRQ(stat);
    stat = put_dimension(zb, grid->Mbz, grid->zblevels); CHKERRQ(stat);
  }

  return 0;
}

//! Appends \c time to the t dimension.
PetscErrorCode NCTool::append_time(PetscReal time) {
  int stat, t_id;


  if (grid->rank == 0) {
    size_t t_len;
    stat = nc_inq_dimid(ncid, "t", &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimlen(ncid, t_id, &t_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_inq_varid(ncid, "t", &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_put_var1_double(ncid, t_id, &t_len, &time); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  return 0;
}

//! Opens a NetCDF file for reading.
/*!
  Stops if a file does not exist or could not be opened.
 */
PetscErrorCode NCTool::open_for_reading(const char filename[]) {
  PetscErrorCode ierr;
  int stat = 0;
  if (grid->rank == 0) {
    stat = nc_open(filename, NC_NOWRITE, &ncid);
  }
  ierr = MPI_Bcast(&stat, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);

  if (stat != NC_NOERR) {
    ierr = PetscPrintf(grid->com, "ERROR: Can't open file '%s'!\n",
		       filename); CHKERRQ(ierr);
    PetscEnd();
  }
  
  return 0;
}

//! Closes a NetCDF file.
PetscErrorCode NCTool::close() {
  PetscErrorCode ierr;
  if (grid->rank == 0) {
    ierr = nc_close(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }
  ncid = -1;			// make it invalid
  return 0;
}

//! Opens or creates a NetCDF file for writing the PISM model state.
/*!
  1) Open the file. If the file could not be opened, go to step 3. Otherwise go
  to step 2.

  2) Call check_dimensions(). If dimensions are OK, go to step 5. Otherwise go
  to step 4.

  3) Create the file with NC_CLOBBER (overwrite). Go to step 4.

  4) Define dimensions and create dimension variables. Go to step 5.

  5) Set the NetCDF ID.
 */
PetscErrorCode NCTool::open_for_writing(const char filename[], bool replace) {
  int stat;

  if (grid->rank == 0) {
    bool file_exists = false, dimensions_are_ok = false;
    char tmp[PETSC_MAX_PATH_LEN];

    // Check if the file exists:
    if (FILE *f = fopen(filename, "r")) {
      file_exists = true;
      fclose(f);
    } else {
      file_exists = false;
    }

    if (file_exists) {
      if (replace) {
	strcpy(tmp, filename);
	strcat(tmp, "~");	// tmp <- "foo.nc~"

	stat = rename(filename, tmp);
	if (stat != 0) {
	  stat = verbPrintf(1, grid->com, "PISM ERROR: can't move '%s' to '%s'.\n",
			    filename, tmp);
	  PetscEnd();
	}
	stat = verbPrintf(2, grid->com, 
	   "PISM WARNING: output file '%s' already exists. Moving it to '%s'.\n",
	   filename, tmp);
	file_exists = false;
      } else {
	// Check dimensions
	stat = nc_open(filename, NC_WRITE, &ncid);
	if (stat != NC_NOERR) {
	  stat = verbPrintf(1, grid->com, "PISM ERROR: NetCDF error: %s\n",
			    nc_strerror(stat));
	  PetscEnd();
	}
	dimensions_are_ok = check_dimensions();
      }
    }

    if (!file_exists || !dimensions_are_ok) {
      stat = nc_create(filename, NC_CLOBBER|NC_64BIT_OFFSET, &ncid); 
      CHKERRQ(check_err(stat,__LINE__,__FILE__));
      stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
      stat = create_dimensions(); CHKERRQ(stat);
    }
    stat = nc_set_fill(ncid, NC_NOFILL, NULL); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  } // end of if(grid->rank == 0)

  stat = MPI_Bcast(&ncid, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);
  return 0;
}

//! Writes global attributes to a NetCDF file.
PetscErrorCode NCTool::write_global_attrs(bool have_ssa_velocities, const char conventions[]) {
  int stat, flag = 0;
  char tmp[TEMPORARY_STRING_LENGTH];

  snprintf(tmp, TEMPORARY_STRING_LENGTH, "PISM %s", PISM_Revision);

  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    if (have_ssa_velocities)
      flag = 1;
    stat = nc_put_att_int(ncid, NC_GLOBAL, "ssa_velocities_are_valid", NC_INT, 1, &flag);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", strlen(conventions), conventions); 
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_put_att_text(ncid, NC_GLOBAL, "source", strlen(tmp), tmp); 
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  return 0;
}

//! Finds the length of a dimension. Returns 0 if failed.
PetscErrorCode NCTool::get_dim_length(const char name[], int *len) {
  int stat, dim_id;
  
  if (grid->rank == 0) {
    size_t dim_len;
    stat = nc_inq_dimid(ncid, name, &dim_id);
    if (stat == NC_NOERR) {
      stat = nc_inq_dimlen(ncid, dim_id, &dim_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    } else
      dim_len = 0;

    *len = static_cast<int>(dim_len);
  }

  stat = MPI_Bcast(len, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);

  return 0;
}

//! Gets dimension limits.
/*! Gets dimension limits (i.e min and max values of the coordinate variable).

  Sets min = 0 and max = 0 if the dimension \c name has length 0.

  Set \c min or \c max to NULL to omit the corresponding value.
 */
PetscErrorCode NCTool::get_dim_limits(const char name[], double *min, double *max) {
  PetscErrorCode ierr;
  int len, varid;
  bool variable_exists;
  size_t start = 0, count;
  double range[2] = {0, 0};

  ierr = get_dim_length(name, &len); CHKERRQ(ierr);

  if (len != 0) {
    ierr = find_variable(name, &varid, variable_exists);
    if (!variable_exists) {
      ierr = PetscPrintf(grid->com, "PISM ERROR: coordinate variable '%s' does not exist.\n",
			 name);
      CHKERRQ(ierr);
      PetscEnd();
    }

    if (grid->rank == 0) {
      double *data;
      data = new double[len];

      count = static_cast<size_t>(len);
      ierr = nc_get_vara_double(ncid, varid, &start, &count, data); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

      range[0] = data[0];
      range[1] = data[0];
      for (int j = 1; j < len; j++) {
	range[0] = PetscMin(data[j], range[0]);
	range[1] = PetscMax(data[j], range[1]);
      }
      delete[] data;
    } // end of if(grid->rank == 0)
    ierr = MPI_Bcast(range, 2, MPI_DOUBLE, 0, grid->com); CHKERRQ(ierr);
  } // if (len != 0)

  char internal_units[TEMPORARY_STRING_LENGTH];
  utUnit input, internal;
  bool input_has_units;

  if (strcmp(name, "t") == 0) {
    // Note that this units specification does *not* have a reference date.
    strcpy(internal_units, "seconds");
  } else {
    strcpy(internal_units, "meters");
  }

  if (utScan(internal_units, &internal) != 0) {
    SETERRQ(1, "UDUNITS failed trying to scan internal units.");
  }

  // Get the units information:
  ierr = get_units(varid, input_has_units, input); CHKERRQ(ierr);
  if (!input_has_units) {
    ierr = verbPrintf(3, grid->com,
		      "PISM WARNING: dimensional variable '%s' does not have the units attribute.\n"
		      "     Assuming that it is in '%s'.\n",
		      name, internal_units); CHKERRQ(ierr);
    utCopy(&internal, &input);
  }

  // Find the conversion coefficients:
  double slope, intercept;
  ierr = utConvert(&input, &internal, &slope, &intercept);
  if (ierr != 0) {
    if (ierr == UT_ECONVERT) {
      ierr = PetscPrintf(grid->com,
			 "PISM ERROR: dimensional variable '%s' has the units that are not compatible with '%s'.\n",
			 name, internal_units); CHKERRQ(ierr);
      PetscEnd();
    }
    SETERRQ1(1, "UDUNITS failure: error code = %d", ierr);
  }

  // Change units and return limits:
  if (min != NULL)
    *min = intercept + range[0]*slope;
  if (max != NULL)
    *max = intercept + range[1]*slope;

  return 0;
}

//! Reads a text attribute from a NetCDF file.
/*!
  Missimg and empty attributes are treated the same.
 */
PetscErrorCode NCTool::get_att_text(const int varid, const char name[], string &result) {
  char *str = NULL;
  int ierr, stat, len;

  // Read and broadcast the attribute length:
  if (grid->rank == 0) {
    size_t attlen;
    stat = nc_inq_attlen(ncid, varid, name, &attlen);
    if (stat == NC_NOERR)
      len = static_cast<int>(attlen);
    else
      len = 0;
  }
  ierr = MPI_Bcast(&len, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);

  // Allocate some memory or set result to NULL and return:
  if (len == 0) {
    result = "";
    return 0;
  }
  str = new char[len + 1];
  // Zealously clear the string, so that we don't risk moving unitialized bytes over MPI (because Valgrind can't tell
  // the difference between these harmless bytes and potential memory errors)
  ierr = PetscMemzero(str, len+1);CHKERRQ(ierr);

  // Now read the string and broadcast stat to see if we succeeded:
  if (grid->rank == 0) {
    stat = nc_get_att_text(ncid, varid, name, str);
  }
  ierr = MPI_Bcast(&stat, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  
  // On success, broadcast the string. On error, set str to "".
  if (stat == NC_NOERR) {
    stat = MPI_Bcast(str, len, MPI_CHAR, 0, grid->com); CHKERRQ(stat);
  } else {
    strcpy(str, "");
  }

  result = str;

  delete[] str;
  return 0;
}

//! Reads a scalar attribute from a NetCDF file.
PetscErrorCode NCTool::get_att_double(const int varid, const char name[],
				      vector<double> &result) {
  int ierr, stat, len;

  // Read and broadcast the attribute length:
  if (grid->rank == 0) {
    size_t attlen;
    stat = nc_inq_attlen(ncid, varid, name, &attlen);
    if (stat == NC_NOERR)
      len = static_cast<int>(attlen);
    else
      len = 0;
  }
  ierr = MPI_Bcast(&len, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);

  if (len == 0) {
    result.clear();
    return 0;
  }

  result.reserve(len);
  // Now read the data and broadcast stat to see if we succeeded:
  if (grid->rank == 0) {
    stat = nc_get_att_double(ncid, varid, name, &result[0]);
  }
  ierr = MPI_Bcast(&stat, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  
  // On success, broadcast the data. On error, stop.
  if (stat == NC_NOERR) {
    ierr = MPI_Bcast(&result[0], len, MPI_DOUBLE, 0, grid->com); CHKERRQ(ierr);
  } else {
    SETERRQ(1, "Error reading an attribute.");
  }

  return 0;
}

//! Assembles start and count arrays for a particular variable.
/*!
  Does nothing on processors other than zero.

  Arrays \c start and \c count have to have length 5, to be able to hold
  values for t,x,y,z,zb (in this order).

  Arrays nc_start and nc_count will be allocated *by* this function and have to
  be freed by the user.

  Also note that here X and Y have PISM (internal) meaning.
*/
PetscErrorCode NCTool::compute_start_and_count(const int varid, int *start, int *count,
					       size_t* &nc_start, size_t* &nc_count) {
  int stat, ndims;
  // Indices:
  const int T = 0, X = 1, Y = 2, Z = 3, ZB = 4;
  // IDs of interesting dimensions:
  int t_id, x_id, y_id, z_id, zb_id;
  // IDs of all the dimensions a variables depends on:
  int *dimids;

  // Quit if this is not processor zero:
  if (grid->rank != 0) return 0;

//   for (int j = 0; j < 5; j++)
//     fprintf(stderr, "start[%d] = %d, count[%d] = %d\n", j, start[j], j, count[j]);

  // Get the number of dimensions a variable depends on:
  stat = nc_inq_varndims(ncid, varid, &ndims);
  CHKERRQ(check_err(stat,__LINE__,__FILE__));

  // Allocate all the arrays we need:
  nc_start = new size_t[ndims];
  nc_count = new size_t[ndims];
  dimids   = new int[ndims];

  // Find all the dimensions we care about:
  // t
  stat = nc_inq_dimid(ncid, "t", &t_id);
  if (stat != NC_NOERR) t_id = -1;
  // x
  stat = nc_inq_dimid(ncid, "x", &x_id);
  if (stat != NC_NOERR) x_id = -1;
  // y
  stat = nc_inq_dimid(ncid, "y", &y_id);
  if (stat != NC_NOERR) y_id = -1;
  // z
  stat = nc_inq_dimid(ncid, "z", &z_id);
  if (stat != NC_NOERR) z_id = -1;
  // zb
  stat = nc_inq_dimid(ncid, "zb", &zb_id);
  if (stat != NC_NOERR) zb_id = -1;

  // Get the list of dimensions a variable depends on:
  stat = nc_inq_vardimid(ncid, varid, dimids);
  CHKERRQ(check_err(stat,__LINE__,__FILE__));

  // Assemble nc_start and nc_count:
  for (int j = 0; j < ndims; j++) {
    if (dimids[j] == t_id) {
      nc_start[j] = start[T];
      nc_count[j] = count[T];
    } else if (dimids[j] == x_id) {
      nc_start[j] = start[Y];	// transpose
      nc_count[j] = count[Y];	// transpose
    } else if (dimids[j] == y_id) {
      nc_start[j] = start[X];	// transpose
      nc_count[j] = count[X];	// transpose
    } else if (dimids[j] == z_id) {
      nc_start[j] = start[Z];
      nc_count[j] = count[Z];
    } else if (dimids[j] == zb_id) {
      nc_start[j] = start[ZB];
      nc_count[j] = count[ZB];
    } else {
      nc_start[j] = 0;
      nc_count[j] = 1;
    }
//     fprintf(stderr, "nc_start[%d] = %ld, nc_count[%d] = %ld\n", j, nc_start[j], j, nc_count[j]); 
  }

  delete[] dimids;
  return 0;
}

PetscErrorCode NCTool::transpose(int varid, GridType dim_flag, int *count,
				 double* in, double* out) {
  const int N = 3;		// process at most 3 dimensions
  const int X = 1, Y = 2, Z = 3, ZB = 4;
  int stat, ndims, *dimids;
  int x_id, y_id, z_id = -1;

  if (grid->rank != 0) return 0;

  stat = nc_inq_varndims(ncid, varid, &ndims);
  CHKERRQ(check_err(stat,__LINE__,__FILE__));

  // Find all the dimensions we care about:
  // x
  stat = nc_inq_dimid(ncid, "x", &x_id);
  if (stat != NC_NOERR) x_id = -1;

  // y
  stat = nc_inq_dimid(ncid, "y", &y_id);
  if (stat != NC_NOERR) y_id = -1;

  // z or zb, depending on dim_flag
  if (dim_flag == GRID_3D) {
    stat = nc_inq_dimid(ncid, "z", &z_id);
    if (stat != NC_NOERR) z_id = -1;
  } else if (dim_flag == GRID_3D_BEDROCK) {
    stat = nc_inq_dimid(ncid, "zb", &z_id);
    if (stat != NC_NOERR) z_id = -1;
  }

  // Get the list of dimensions a variable depends on:
  dimids = new int[ndims];
  stat = nc_inq_vardimid(ncid, varid, dimids);
  CHKERRQ(check_err(stat,__LINE__,__FILE__));

  // Create arrays of IDs of spatial dimensions, in PISM and input orders (they
  // are used to create the map):
  int pism_dimids[N] = {y_id, x_id, z_id};
  int input_dimids[N] = {-2, -2, -2}; // invalid and never equal to any of dimids
  int m = 0;
  for (int j = 0; j < ndims; j++) {
    if (dimids[j] == x_id || dimids[j] == y_id || dimids[j] == z_id) {

      if (m == 3) {
	PetscPrintf(grid->com,
		    "PISM ERROR: Can't process a variable because a spatial dimension appears twice in its specification.\n"
		    "     Please see ncdump -h <filename> for details.\n");
	PetscEnd();
      }

      input_dimids[m] = dimids[j];
      m++;
    }
  }
  delete[] dimids;

  // Create the array describing the permutation:
  int Map[N] = {0, 1, 2};
  for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++)
      if (pism_dimids[j] == input_dimids[k])
	Map[j] = k;

//   for (int j = 0; j < N; j++)
//     fprintf(stderr, "Map[%d] = %d\n", j, Map[j]);

  // Create arrays of counts, in PISM and input orders:
  int pism_count[N] = {count[X], count[Y], 1};
  if (dim_flag == GRID_3D)         pism_count[2] = count[Z];
  if (dim_flag == GRID_3D_BEDROCK) pism_count[2] = count[ZB];

  int input_count[N];
  for (int j = 0; j < N; j++)
    input_count[j] = pism_count[Map[j]];

  // Permute dimensions:
  for (int i = 0; i < pism_count[0]; i++) {
    for (int j = 0; j < pism_count[1]; j++) {
      for (int k = 0; k < pism_count[2]; k++) {
	int A[N] = {i,j,k};
	// row-column-level-style indices in the input array:
	int I = A[Map[0]];
	int J = A[Map[1]];
	int K = A[Map[2]];
	// index in the input array:
	int N = (I * input_count[1] + J) * input_count[2] + K;
	// index in the PISM array:
	int M = (i *  pism_count[1] + j) *  pism_count[2] + k;
	
	out[M] = in[N];
      }
    }
  }

  return 0;
}

//! Computes the size of the local block.
int NCTool::compute_block_size(GridType dims, int* count) {
  const int X = 1, Y = 2, Z = 3, ZB = 4;
  switch (dims) {
  case GRID_2D:
    return count[X] * count[Y];
  case GRID_3D:
    return count[X] * count[Y] * count[Z];
  case GRID_3D_BEDROCK:
    return count[X] * count[Y] * count[ZB];
  }
  return 0;
}

//! Initializes the IceGrid object from a NetCDF file.
PetscErrorCode NCTool::get_grid(const char filename[]) {
  PetscErrorCode ierr;
  grid_info gi;
  double *z_levels, *zb_levels;

  ierr = open_for_reading(filename); CHKERRQ(ierr);

  ierr = get_grid_info(gi); CHKERRQ(ierr);
  ierr = get_vertical_dims(z_levels, zb_levels); CHKERRQ(ierr);

  grid->Mx = gi.x_len;
  grid->My = gi.y_len;
  grid->Lx = gi.Lx;
  grid->Ly = gi.Ly;
  grid->periodicity = NONE;
  grid->x0   = gi.x0;
  grid->y0   = gi.y0;
  grid->year = gi.time / secpera;

  ierr = grid->compute_horizontal_spacing(); CHKERRQ(ierr);
  ierr = grid->set_vertical_levels(gi.z_len, gi.zb_len, z_levels, zb_levels); CHKERRQ(ierr);

  // We're ready to call grid->createDA().

  // Cleanup:
  ierr = close(); CHKERRQ(ierr);
  delete[] z_levels;
  delete[] zb_levels;
  return 0;
}

//! Get variable's units information from a NetCDF file.
/*!
  Note that this function intentionally ignores the reference date.
 */
PetscErrorCode NCTool::get_units(int varid, bool &has_units, utUnit &units) {
  PetscErrorCode ierr;
  string units_string;

  // Get the string:
  ierr = get_att_text(varid, "units", units_string); CHKERRQ(ierr);

  // If a variables does not have the units attribute, set the flag and return:
  if (units_string.empty()) {
    has_units = false;
    utClear(&units);
    return 0;
  }

  // This finds the string "since" in the units_string and terminates
  // it on the first 's' of "since", if this sub-string was found. This
  // is done to ignore the reference date in the time units string (the
  // reference date specification always starts with this word).

  int n = units_string.find("since");
  if (n != -1) units_string.resize(n);

  ierr = utScan(units_string.c_str(), &units);
  if (ierr != 0) {
    ierr = PetscPrintf(grid->com, "PISM ERROR: units specification '%s' is unknown or invalid.\n",
		       units_string.c_str());
    PetscEnd();
  }
  
  has_units = true;

  // Cleanup:
  return 0;
}

//! Creates a time-series variable (i.e. one that depends on "t" only).
PetscErrorCode NCTool::create_timeseries(const char name[], const char long_name[],
					 const char units[],
					 nc_type nctype, int *varid) {
  int stat, t_id, var_id;

  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "t", &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_def_var(ncid, name, nctype, 1, &t_id, &var_id);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    if (units != NULL) {
      stat = nc_put_att_text(ncid, var_id, "units", strlen(units), units);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));
    }

    if (long_name != NULL) {
      stat = nc_put_att_text(ncid, var_id, "long_name", strlen(long_name), long_name);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));
    }

    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  stat = MPI_Bcast(&var_id, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);

  if (varid != NULL)
    *varid = var_id;

  return 0;
}

//! Writes \c value to a NetCDF variable \c name(t). Uses the last time slice.
PetscErrorCode NCTool::append_timeseries(const char name[], double value) {
  int stat, varid, t_len;
  bool variable_exists;

  stat = get_dim_length("t", &t_len); CHKERRQ(stat);

  stat = find_variable(name, &varid, variable_exists); CHKERRQ(stat);
  if (!variable_exists) {
    stat = create_timeseries(name, NULL, NULL, NC_FLOAT, &varid); CHKERRQ(stat);
  }

  if (grid->rank == 0) {
    size_t index = static_cast<size_t>(t_len - 1);
    stat = nc_put_var1_double(ncid, varid, &index, &value); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  return 0;
}

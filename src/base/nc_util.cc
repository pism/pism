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



NCTool::NCTool() {
  //FIXME: does there need to be a default MaskInterp?
  myMaskInterp = PETSC_NULL;
  grid = PETSC_NULL;
  ncid = -1;

  // does not return a value
}

NCTool::NCTool(IceGrid *my_grid) {
  //FIXME: does there need to be a default MaskInterp?
  myMaskInterp = PETSC_NULL;
  grid = my_grid;
  ncid = -1;
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
PetscErrorCode  NCTool::find_variable(const char short_name[], const char standard_name[],
				      int *varidp, bool &exists) {
  int ierr;
  size_t attlen;
  int stat, found = 0, my_varid = -1, nvars;
  char attribute[TEMPORARY_STRING_LENGTH];

  // Processor 0 does all the job here.
  if (grid->rank == 0) {
    if (standard_name != NULL) {
      ierr = nc_inq_nvars(ncid, &nvars); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

      for (int j = 0; j < nvars; j++) {
	stat = nc_get_att_text(ncid, j, "standard_name", attribute);
	if (stat != NC_NOERR) {
	  continue;
	}

	// text attributes are not always zero-terminated, so we need to add the
	// trailing zero:
	stat = nc_inq_attlen(ncid, j, "standard_name", &attlen);
	CHKERRQ(check_err(stat,__LINE__,__FILE__));
	attribute[attlen] = 0;

	if (strcmp(attribute, standard_name) == 0) {
	  if (!found) {		// if unique
	    found = 1;
	    my_varid = j;
	  } else {    // if not unique
 	    ierr = PetscPrintf(grid->com,
			       "PISM ERROR: Inconsistency in the input file: "
			       "Variables #%d and #%d have the same standard_name ('%s').\n",
			       my_varid, j, attribute);
	    CHKERRQ(ierr);
	    PetscEnd();
	  }
	}
      }
    } // end of if(standard_name != NULL)

    if (!found) {
      // look for short_name
      stat = nc_inq_varid(ncid, short_name, &my_varid);
      if (stat == NC_NOERR)
	found = 1;
    }
  } // end of if(grid->rank == 0)

  // Broadcast the existence flag and the variable ID.
  ierr = MPI_Bcast(&found, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&my_varid, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);

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

// FIXME: Consider removing this method.
PetscErrorCode NCTool::set_grid(IceGrid *my_grid) {
  grid = my_grid;
  return 0;
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
  ierr = get_dim_length("x",  &g.x_len); CHKERRQ(ierr);
  ierr = get_dim_length("y",  &g.y_len); CHKERRQ(ierr);

  ierr = get_dim_limits("t", NULL, &g.time); CHKERRQ(ierr);
  ierr = get_dim_limits("x", &g.x_min, &g.x_max); CHKERRQ(ierr);
  ierr = get_dim_limits("y", &g.y_min, &g.y_max); CHKERRQ(ierr);

  return 0;
}

//! Read the last value of the time variable t from a NetCDF file.
PetscErrorCode NCTool::get_last_time(double *time) {
  PetscErrorCode ierr;

  ierr = get_dim_limits("t", NULL, time);

  return 0;
}


//! Read in the variables \c z and \c zb from the NetCDF file; <i>do not</i> assume they are equally-spaced.
/*!
Arrays z_read[] and zb_read[] must be pre-allocated arrays of length at least z_len, zb_len, respectively.

Get values of z_len, zb_len by using get_grid_info() before this routine; z_len
= grid_info.z_len and zb_len = grid_info.zb_len.
 */
PetscErrorCode NCTool::get_vertical_dims(const int z_len, const int zb_len,
                                         double z_read[], double zb_read[]) {
  int stat;
  int z_id, zb_id;
  size_t zeroST  = 0,
         zlenST  = (size_t) z_len,
         zblenST = (size_t) zb_len;

  if (grid->rank == 0) {
    stat = nc_inq_varid(ncid, "z", &z_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_varid(ncid, "zb", &zb_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_get_vara_double(ncid, z_id, &zeroST, &zlenST, z_read);
             CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_vara_double(ncid, zb_id, &zeroST, &zblenST, zb_read);
             CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  MPI_Bcast(z_read, z_len, MPI_DOUBLE, 0, grid->com);
  MPI_Bcast(zb_read, zb_len, MPI_DOUBLE, 0, grid->com);
  return 0;
}


//! Put the variable for a dimension in a NetCDF file.  Uses starting value and a spacing for regularly-spaced values.
/*!
Note the variable corresponding to a dimension is always \c double in a PISM NetCDF file.
 */
PetscErrorCode NCTool::put_dimension_regular(int varid, int len, double start, double delta) {
  PetscErrorCode ierr;
  int stat;
  double *v;

  ierr = PetscMalloc(len * sizeof(double), &v); CHKERRQ(ierr);
  for (int i = 0; i < len; i++) {
    v[i] = start + i * delta;
  }
  stat = nc_put_var_double(ncid, varid, v); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  ierr = PetscFree(v); CHKERRQ(ierr);

  return 0;
}


//! Put the variable for a dimension in a NetCDF file.  Makes no assumption about spacing.
/*!
Note the variable corresponding to a dimension is always \c double in a PISM NetCDF file.
 */
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
PetscErrorCode NCTool::get_local_var(const int varid, DA da, Vec v, Vec g,
				     const int *s, const int *c,
				     int dims, void *a_mpi, int a_size) {

  PetscErrorCode ierr;
  ierr = get_global_var(varid, da, g, s, c, dims, a_mpi, a_size);
            CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, g, INSERT_VALUES, v); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, g, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


//! Read a variable in a NetCDF file into a \c DA -managed global \c Vec \c g.  \e In \e parallel.
PetscErrorCode NCTool::get_global_var(const int varid, DA da, Vec g,
				      const int *s, const int *c,
				      int dims, void *a_mpi, int a_size) {

  const int req_tag = 1; // MPI tag for request block
  const int var_tag = 2; // MPI tag for data block
  const int sc_size = 8;
  PetscErrorCode ierr;
  MPI_Status mpi_stat;
  int stat;
  int sc[sc_size]; // buffer to hold both `s' and `c'
  size_t sc_nc[sc_size];
  double *a_double = NULL;

  a_double = (double *)a_mpi;

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

      stat = nc_get_vara_double(ncid, varid, &sc_nc[0], &sc_nc[dims], a_double);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      if (proc != 0) {
        int b_len = 1;
        for (int i = dims; i < 2 * dims; i++) b_len *= sc[i];
	MPI_Send(a_double, b_len, MPI_DOUBLE, proc, var_tag, grid->com);
      }
    }
  } else {
    MPI_Send(sc, 2 * dims, MPI_INT, 0, req_tag, grid->com);
    MPI_Recv(a_double, a_size, MPI_DOUBLE, 0, var_tag, grid->com, &mpi_stat);
  }

  int b_len = 1;
  for (int i = dims; i < 2 * dims; i++) b_len *= sc[i];
  PetscScalar *a_petsc;
  ierr = VecGetArray(g, &a_petsc); CHKERRQ(ierr);
  for (int i = 0; i < b_len; i++) {
    a_petsc[i] = (PetscScalar)a_double[i];
  }

  ierr = VecRestoreArray(g, &a_petsc); CHKERRQ(ierr);
  return 0;
}


//! Put a \c DA -managed local \c Vec \c v into a variable in a NetCDF file; a global \c Vec \c g is used for storage space.
/*!
Just calls put_global_var(), after transfering the local \c Vec called \c vec into the global \c Vec \c g.
 */
PetscErrorCode NCTool::put_local_var(const int varid, DA da, Vec v, Vec g,
				     const int *s, const int *c,
                                     int dims, void *a_mpi, int a_size) {

  PetscErrorCode ierr;
  ierr = DALocalToGlobal(da, v, INSERT_VALUES, g); CHKERRQ(ierr);
  ierr = put_global_var(varid, da, g, s, c, dims, a_mpi, a_size); CHKERRQ(ierr);
  return 0;
}


//! Put a \c DA -managed global \c Vec \c g into a variable in a NetCDF file.  \e In \e parallel.
PetscErrorCode NCTool::put_global_var(const int varid, DA da, Vec g,
				      const int *s, const int *c,
                                      int dims, void *a_mpi, int a_size) {
  const int lim_tag = 1; // MPI tag for limits block
  const int var_tag = 2; // MPI tag for data block
  const int sc_size = 8;
  PetscErrorCode ierr;
  MPI_Status mpi_stat;
  int stat;
  int sc[sc_size];		// buffer to hold both `s' (start) and `c' (count)
  size_t sc_nc[sc_size];
  double *a_double = NULL;

  a_double = (double *)a_mpi;

  for (int i = 0; i < 2 * dims; i++) {
    sc[i] = (i < dims) ? s[i] : c[i - dims];
  }

  int b_len = 1;
  for (int i = 0; i < dims; i++) b_len *= c[i];

  // convert IceModel Vec containing PetscScalar to array of double for NetCDF
  PetscScalar *a_petsc;
  ierr = VecGetArray(g, &a_petsc); CHKERRQ(ierr);
  for (int i = 0; i < b_len; i++) {
    a_double[i] = (double)a_petsc[i];
  }
  ierr = VecRestoreArray(g, &a_petsc); CHKERRQ(ierr);

  if (grid->rank == 0) { // on rank 0 processor, receive messages from every other
                         //    processor, then write it out to the NC file
    for (int proc = 0; proc < grid->size; proc++) { // root will write itself last
      if (proc != 0) {
        MPI_Recv(sc, sc_size, MPI_INT, proc, lim_tag, grid->com, &mpi_stat);
	MPI_Recv(a_double, a_size, MPI_DOUBLE, proc, var_tag, grid->com, &mpi_stat);
      }

//       {
//         printf("[%1d] writing %4d [", proc, varid);
//         for (int i = 0; i < 2 * dims; i++) {
//           if (i == dims) printf("] [");
//           printf("%5d", sc[i]);
//         }
//         printf("]\n");
//       }

      for (int i = 0; i < 2 * dims; i++) sc_nc[i] = (size_t)sc[i]; // we need size_t

      stat = nc_put_vara_double(ncid, varid, &sc_nc[0], &sc_nc[dims], a_double);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));
    }
  } else {  // all other processors send to rank 0 processor
    MPI_Send(sc, 2 * dims, MPI_INT, 0, lim_tag, grid->com);
    MPI_Send(a_double, a_size, MPI_DOUBLE, 0, var_tag, grid->com);
  }
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
                                        DA da, Vec vec, Vec g, 
                                        bool useMaskInterp) {
  PetscErrorCode ierr;
  ierr = regrid_global_var(varid, dim_flag, lic, da, g, useMaskInterp); CHKERRQ(ierr);
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
PetscErrorCode NCTool::regrid_global_var(const int varid, GridType dim_flag, LocalInterpCtx &lic,
                                         DA da, Vec g, bool useMaskInterp) {
  PetscErrorCode ierr;
  const int N = 4, X = 1, Y = 2, Z = 3, ZB = 4; // indices, just for clarity
  const int start_tag = 1; // MPI tag for the start array
  const int count_tag = 2; // MPI tag for the count array
  const int data_tag  = 3; // MPI tag for the data block
  MPI_Status mpi_stat;
  int stat, dims = 0;
  int start[N], count[N];	// enough space for t, x, y, z (or zb)
  bool t_exists;
 
  switch (dim_flag) {
    case GRID_2D:
      dims = 3; // time, x, y
      break;
    case GRID_3D:
      dims = 4; // time, x, y, z
      if (lic.regrid_2d_only) {
        SETERRQ(2, "regrid_2d_only is set, so dim_flag == GRID_3D is not allowed");
      }
      break;
    case GRID_3D_BEDROCK:
      dims = 4; // time, x, y, zb
      if (lic.no_regrid_bedrock) {
        SETERRQ(2, "no_regrid_bedrock is set, so dim_flag == GRID_3D_BEDROCK is not allowed");
      }
  }

  // make local copies of lic.start and lic.count
  for (int i = 0; i < N; i++) {
    start[i] = lic.start[i];
    count[i] = lic.count[i];
  }

  // At this point, start and count are set up correctly for ice 3-D quantities.  In order
  // to avoid duplicating lots of code below, we play some tricks here.
  if (dim_flag == GRID_2D) { // 2-D quantity
    // Treat it as a 3-D quantity with extent 1 in the z-direction.  The netCDF
    // read actually will not use this information (it is passed in, but it
    // won't index these entries in the array, because it knows how many
    // dimensions it there are), but count[3] will be used when node 0 computes how
    // much data it has to send back.  start[3] should never be used in this case.
    start[Z] = 0;
    count[Z] = 1;
  } else if (dim_flag == GRID_3D_BEDROCK) { // Bedrock quantity
    // Just fill in the appropriate values
    start[Z] = lic.start[ZB];
    count[Z] = lic.count[ZB];
  }

  // The next line appears outside the if (grid->rank == 0) since it calls
  // MPI_Bcast.
  ierr = find_dimension("t", NULL, t_exists); CHKERRQ(ierr);

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

      // It is not safe to cast memory.  In particular on amd64 int and size_t are
      // different sizes.  Since netCDF uses size_t for the offsets, we need to here too.
      size_t start_nc[N], count_nc[N];

      for (int i = 0; i < N; i++) {
	start_nc[i] = (size_t) start[i];
	count_nc[i] = (size_t) count[i];
      }

      // Actually read the block into the buffer.
      if (t_exists) {
	stat = nc_get_vara_double(ncid, varid, start_nc, count_nc, lic.a);
	CHKERRQ(check_err(stat,__LINE__,__FILE__));
      } else {
	// This ignores 't' values of start and count:
	stat = nc_get_vara_double(ncid, varid, &start_nc[1], &count_nc[1], lic.a);
	CHKERRQ(check_err(stat,__LINE__,__FILE__));
      }

      /*{
        printf("ncid, varid = %d %d\n", ncid, varid);
        printf("a[] ni {%f %f %f %f}\n", lic.a[50], lic.a[150], lic.a[250], lic.a[350]);

        const int blen = 9;
        float *buf = (float *)malloc(blen * sizeof(float));
        sc_nc[5] = sc_nc[6] = 3; sc_nc[7] = 1;
        stat = nc_get_vara_float(ncid, varid, &sc_nc[0], &sc_nc[4], buf);
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
      for (int i = 0; i < dims; i++) a_len *= count[i];

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
  const int ycount = lic.count[Y];
  int myMz = 1, zcount = 1; // indexing trick so that we don't have to duplicate code for the 2-D case.
  if (dim_flag == GRID_3D) {
    myMz = grid->Mz;
    zcount = lic.count[Z];
  } else if (dim_flag == GRID_3D_BEDROCK) {
    myMz = grid->Mbz;
    zcount = lic.count[ZB];
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

//! Writes \c history to the history global attribute of a NetCDF dataset.
/*!
  Appends if overwrite == false (default).
 */
PetscErrorCode NCTool::write_history(const char history[], bool overwrite) {
  int stat;
  size_t history_len;
  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_attlen(ncid, NC_GLOBAL, "history", &history_len);

    if ((stat == NC_ENOTATT) || overwrite) { // no history attribute present or
					     // asked to overwrite

      stat = nc_put_att_text(ncid, NC_GLOBAL, "history",
			     strlen(history), history);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

    } else {

      CHKERRQ(check_err(stat,__LINE__,__FILE__)); // check for other errors

      history_len += 1;		// add the space for the trailing zero
      // allocate some memory
      char *tmp, *result;
      int length = strlen(history) + history_len;
      result = new char[length];
      tmp = new char[history_len];

      stat = nc_get_att_text(ncid, NC_GLOBAL, "history", tmp);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));
      tmp[history_len - 1] = 0;	// terminate the string

      // concatenate strings
      strcpy(result, history);
      strcat(result, tmp);

      // write
      stat = nc_put_att_text(ncid, NC_GLOBAL, "history",
			     strlen(result), result);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      delete[] tmp;
      delete[] result;
    }
    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  } // if (grid->rank == 0)
  
  return 0;
}

PetscErrorCode NCTool::read_polar_stereographic(double &straight_vertical_longitude_from_pole,
						double &latitude_of_projection_origin,
						double &standard_parallel,
						bool report) {
  PetscErrorCode ierr;
  double lon, lat, par;
  int stat, varid, lon_exists, lat_exists, par_exists;
  bool ps_exists;

  ierr = find_variable("polar_stereographic", NULL, &varid, ps_exists); CHKERRQ(ierr);
  if (!ps_exists && report) {
    ierr = verbPrintf(2,grid->com,
		      "  polar stereo not found, using defaults: svlfp=%6.2f, lopo=%6.2f, sp=%6.2f\n",
		      straight_vertical_longitude_from_pole,
		      latitude_of_projection_origin,
		      standard_parallel); CHKERRQ(ierr);
    return 0;
  }

  if (grid->rank == 0) {
    stat = nc_inq_attid(ncid, varid, "straight_vertical_longitude_from_pole",
			NULL);
    lon_exists = (stat == NC_NOERR);
    stat = nc_inq_attid(ncid,varid,"latitude_of_projection_origin", NULL);
    lat_exists = (stat == NC_NOERR);
    stat = nc_inq_attid(ncid,varid,"standard_parallel", NULL);
    par_exists = (stat == NC_NOERR);
  }
  ierr = MPI_Bcast(&lon_exists, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&lat_exists, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&par_exists, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);

  if (grid->rank == 0) {
    if (lon_exists) {
      stat = nc_get_att_double(ncid, varid, "straight_vertical_longitude_from_pole",
			       &lon); CHKERRQ(nc_check(stat));
    } 
    if (lon_exists) {
      stat = nc_get_att_double(ncid, varid, "latitude_of_projection_origin",
			       &lat); CHKERRQ(nc_check(stat));
    }
    if (par_exists) {
      stat = nc_get_att_double(ncid, varid, "standard_parallel", &par);
      CHKERRQ(nc_check(stat));
    }
  }
  ierr = MPI_Bcast(&lon, 1, MPI_DOUBLE, 0, grid->com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&lat, 1, MPI_DOUBLE, 0, grid->com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&par, 1, MPI_DOUBLE, 0, grid->com); CHKERRQ(ierr);

  if (report) {
    ierr = verbPrintf(2, grid->com,
		      "  polar stereographic var found; attributes present: svlfp=%d, lopo=%d, sp=%d\n"
		      "     values: svlfp = %6.2f, lopo = %6.2f, sp = %6.2f\n",
		      lon_exists, lat_exists, par_exists,
		      lon, lat, par); CHKERRQ(ierr);
  }

  if (lon_exists)
    straight_vertical_longitude_from_pole = lon;
  if (lat_exists)
    latitude_of_projection_origin = lat;
  if (par_exists)
    standard_parallel = par;
  return 0;
}

PetscErrorCode NCTool::write_polar_stereographic(double straight_vertical_longitude_from_pole,
						 double latitude_of_projection_origin,
						 double standard_parallel) {
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
			     &straight_vertical_longitude_from_pole);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_put_att_double(ncid, varid, "latitude_of_projection_origin", NC_DOUBLE, 1,
			     &latitude_of_projection_origin);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_put_att_double(ncid, varid, "standard_parallel", NC_DOUBLE, 1,
			     &standard_parallel);
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
  x  = check_dimension("x", grid->Mx);
  y  = check_dimension("y", grid->My);
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
    stat = nc_put_att_text(ncid, t, "units", 33, "seconds since 2007-01-01 00:00:00"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, t, "calendar", 4, "none"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, t, "axis", 1, "T"); check_err(stat,__LINE__,__FILE__);
    // x
    stat = nc_def_dim(ncid, "x", grid->Mx, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_def_var(ncid, "x", NC_DOUBLE, 1, &dimid, &x); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_put_att_text(ncid, x, "axis", 1, "X"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, x, "long_name", 32, "x-coordinate in Cartesian system"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, x, "standard_name", 23, "projection_x_coordinate"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, x, "units", 1, "m"); check_err(stat,__LINE__,__FILE__);
    // y
    stat = nc_def_dim(ncid, "y", grid->My, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_def_var(ncid, "y", NC_DOUBLE, 1, &dimid, &y); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_put_att_text(ncid, y, "axis", 1, "Y"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, y, "long_name", 32, "y-coordinate in Cartesian system"); check_err(stat,__LINE__,__FILE__);
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
    stat = put_dimension_regular(x, grid->Mx, -grid->Lx, grid->dx); CHKERRQ(stat);
    stat = put_dimension_regular(y, grid->My, -grid->Ly, grid->dy); CHKERRQ(stat);
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
PetscErrorCode NCTool::open_for_reading(const char filename[], bool &exists) {
  PetscErrorCode ierr;
  int stat = 0;
  if (grid->rank == 0) {
    ierr = nc_open(filename, NC_NOWRITE, &ncid);
    if (ierr == NC_NOERR)
      stat = 1;
  }
  ierr = MPI_Bcast(&stat, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  
  exists = (bool) stat;
  
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
	             "\nPISM WARNING: output file '%s' already exists. Moving it to '%s'.\n\n",
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
  }

  stat = MPI_Bcast(&ncid, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);
  return 0;
}

//! Writes global attributes to a NetCDF file.
PetscErrorCode NCTool::write_global_attrs(bool have_ssa_velocities, const char conventions[]) {
  int stat, flag = 0;

  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    if (have_ssa_velocities)
      flag = 1;
    stat = nc_put_att_int(ncid, NC_GLOBAL, "ssa_velocities_are_valid", NC_INT, 1, &flag);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", strlen(conventions), conventions); 
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

    *len = dim_len;
  }

  stat = MPI_Bcast(len, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);

  return 0;
}

//! Gets dimension limits.
/*! Gets dimension limits (i.e min and max values of the coordinate variable).
  Assumes that the coordinate variable corresponding to dimension \c name is
  strictly ascending.

  Sets min = 0 and max = 0 if the dimension \c name has length 0.

  Set \c min or \c max to NULL to omit the corresponding value.
 */
PetscErrorCode NCTool::get_dim_limits(const char name[], double *min, double *max) {
  PetscErrorCode ierr;
  int len, varid;
  bool variable_exists;
  double range[2] = {0, 0};

  ierr = get_dim_length(name, &len); CHKERRQ(ierr);

  if (len != 0) {
    ierr = find_variable(name, NULL, &varid, variable_exists);
    if (!variable_exists) {
      ierr = PetscPrintf(grid->com, "NCTool::get_dim_limits: coordinate variable '%s' does not exist.\n",
			 name);
      CHKERRQ(ierr);
      PetscEnd();
    }

    if ((grid->rank == 0)) {
      size_t index[2] = {0, len - 1};
      ierr = nc_get_var1_double(ncid, varid, &index[0], &range[0]); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
      ierr = nc_get_var1_double(ncid, varid, &index[1], &range[1]); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }
    ierr = MPI_Bcast(range, 2, MPI_DOUBLE, 0, grid->com); CHKERRQ(ierr);
  }

  if (min != NULL)
    *min = range[0];
  if (max != NULL)
    *max = range[1];

  return 0;
}

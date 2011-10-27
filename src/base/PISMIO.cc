#include "PISMIO.hh"
#include "pism_const.hh"

PISMIO::PISMIO(IceGrid *my_grid)
  : NCTool(my_grid->com, my_grid->rank) {
  grid = my_grid;
  myMaskInterp = NULL;
}

//! Read the first and last values, and the lengths, of the x,y,z,zb dimensions from a NetCDF file.  Read the last t.
PetscErrorCode PISMIO::get_grid_info(grid_info &g) const {
  PetscErrorCode ierr;

  ierr = get_grid_info_2d(g); CHKERRQ(ierr);

  ierr = get_dim_length("z",  &g.z_len);  CHKERRQ(ierr);
  ierr = get_dim_length("zb", &g.zb_len); CHKERRQ(ierr);
  ierr = get_dim_limits("zb", &g.zb_min, NULL);     CHKERRQ(ierr);
  ierr = get_dim_limits("z",  NULL,      &g.z_max); CHKERRQ(ierr);
  
  return 0;
}

//! Read the first and last values, and the lengths, of the x,y dimensions from a NetCDF file.  Read the last t.
PetscErrorCode PISMIO::get_grid_info_2d(grid_info &g) const {
  PetscErrorCode ierr;

  ierr = get_dim_length("t",  &g.t_len); CHKERRQ(ierr);
  ierr = get_dim_length("x",  &g.x_len); CHKERRQ(ierr);
  ierr = get_dim_length("y",  &g.y_len); CHKERRQ(ierr);

  ierr = get_dim_limits("t", NULL, &g.time); CHKERRQ(ierr);
  ierr = get_dim_limits("x", &g.x_min, &g.x_max); CHKERRQ(ierr);
  ierr = get_dim_limits("y", &g.y_min, &g.y_max); CHKERRQ(ierr);

  g.x0 = (g.x_max + g.x_min) / 2.0;
  g.y0 = (g.y_max + g.y_min) / 2.0;

  g.Lx = (g.x_max - g.x_min) / 2.0;
  g.Ly = (g.y_max - g.y_min) / 2.0;

  return 0;
}

//! Read a variable in a NetCDF file into a \c DA -managed global \c Vec \c g.  \e In \e parallel.
PetscErrorCode PISMIO::get_var(const int varid, Vec g, GridType dims, int t) const {
  const int N = 5;
  const int start_tag = 1;
  const int count_tag = 2;
  const int data_tag =  3;
  PetscErrorCode ierr;
  MPI_Status mpi_stat;
  int stat;
  double *a_double = NULL;

  if (grid == NULL) SETERRQ(1, "PISMIO::get_global_var(...): grid == NULL");

  int start[N] = {t, grid->xs, grid->ys, 0,        0};
  int count[N] = {1, grid->xm, grid->ym, grid->Mz, grid->Mbz};

  // Find the local size:
  int block_size, buffer_size;
  block_size = compute_block_size(dims, count);
  // And the maximum size of the data block:
  buffer_size = block_size;
  ierr = MPI_Reduce(&block_size, &buffer_size, 1, MPI_INT, MPI_MAX, 0, com); CHKERRQ(ierr);

  // Memory allocation:
  ierr = PetscMalloc(buffer_size * sizeof(double), &a_double); CHKERRQ(ierr);

  if (rank == 0) {
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
        MPI_Recv(start, N, MPI_INT, proc, start_tag, com, &mpi_stat);
        MPI_Recv(count, N, MPI_INT, proc, count_tag, com, &mpi_stat);
      }

      size_t *nc_start, *nc_count;
      ptrdiff_t *imap;
      ierr = compute_start_and_count(varid, start, count, dims, nc_start, nc_count, imap); CHKERRQ(ierr);

      stat = nc_get_varm_double(ncid, varid, nc_start, nc_count, NULL, imap, a_double);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      delete[] nc_start;
      delete[] nc_count;
      delete[] imap;

      if (proc != 0) {
	block_size = compute_block_size(dims, count);
	MPI_Send(a_double, block_size, MPI_DOUBLE, proc, data_tag, com);
      }
    }
  } else {
    MPI_Send(start, N, MPI_INT, 0, start_tag, com);
    MPI_Send(count, N, MPI_INT, 0, count_tag, com);
    MPI_Recv(a_double, buffer_size, MPI_DOUBLE, 0, data_tag, com, &mpi_stat);
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



//! Put a \c DA -managed global \c Vec \c g into a variable in a NetCDF file.  \e In \e parallel.
PetscErrorCode PISMIO::put_var(const int varid, Vec g, GridType dims) const {
  const int start_tag = 1;
  const int count_tag = 2;
  const int data_tag =  3;
  const int N = 5;
  PetscErrorCode ierr;
  MPI_Status mpi_stat;
  int stat;
  double *a_double = NULL;

  if (grid == NULL) SETERRQ(1, "PISMIO::put_global_var(...): grid == NULL");

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
  ierr = MPI_Reduce(&block_size, &buffer_size, 1, MPI_INT, MPI_MAX, 0, com); CHKERRQ(ierr);

  // Memory allocation:
  ierr = PetscMalloc(buffer_size * sizeof(double), &a_double); CHKERRQ(ierr);

  // convert IceModel Vec containing PetscScalar to array of double for NetCDF
  PetscScalar *a_petsc;
  ierr = VecGetArray(g, &a_petsc); CHKERRQ(ierr);
  for (int i = 0; i < block_size; i++) {
    a_double[i] = (double)a_petsc[i];
  }
  ierr = VecRestoreArray(g, &a_petsc); CHKERRQ(ierr);

  if (rank == 0) { // on rank 0 processor, receive messages from every other
                         //    processor, then write it out to the NC file
    for (int proc = 0; proc < grid->size; proc++) { // root will write itself last
      if (proc != 0) {
        MPI_Recv(start, N, MPI_INT, proc, start_tag, com, &mpi_stat);
        MPI_Recv(count, N, MPI_INT, proc, count_tag, com, &mpi_stat);
	MPI_Recv(a_double, buffer_size, MPI_DOUBLE, proc, data_tag, com, &mpi_stat);
      }

      // We initialize the 'stride' array instead of passing a NULL pointer to
      // nc_put_varm_double below to avoid a segfault caused by a bug in NetCDF
      // 4.1.2 and 4.1.3. Please see https://www.unidata.ucar.edu/jira/browse/NCF-88 for details.
      vector<ptrdiff_t> stride(N);

      for (int i = 0; i < N; ++i)
        stride[i] = 1;
      
      size_t *nc_start, *nc_count;
      ptrdiff_t* imap;
      ierr = compute_start_and_count(varid, start, count, dims, nc_start, nc_count, imap); CHKERRQ(ierr);

      stat = nc_put_varm_double(ncid, varid, nc_start, nc_count, &stride[0], imap, a_double);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      delete[] nc_start;
      delete[] nc_count;
      delete[] imap;
    }
  } else {  // all other processors send to rank 0 processor
    MPI_Send(start, N, MPI_INT, 0, start_tag, com);
    MPI_Send(count, N, MPI_INT, 0, count_tag, com);
    MPI_Send(a_double, buffer_size, MPI_DOUBLE, 0, data_tag, com);
  }

  // Cleanup:
  ierr = PetscFree(a_double); CHKERRQ(ierr);
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
PetscErrorCode PISMIO::regrid_var(const int varid, GridType dims, LocalInterpCtx &lic,
				  Vec g, bool useMaskInterp) const {
  PetscErrorCode ierr;
  const int N = 5, X = 1, Y = 2, Z = 3, ZB = 4; // indices, just for clarity
  const int start_tag = 1; // MPI tag for the start array
  const int count_tag = 2; // MPI tag for the count array
  const int data_tag  = 3; // MPI tag for the data block
  MPI_Status mpi_stat;
  int stat, start[N], count[N];	// enough space for t, x, y, z, zb

  if (grid == NULL) SETERRQ(1, "PISMIO::regrid_global_var(...): grid == NULL");

  // make local copies of lic.start and lic.count
  for (int i = 0; i < N; i++) {
    start[i] = lic.start[i];
    count[i] = lic.count[i];
  }

  // This ensures that we don't read too many values (i.e. if varid refers to a
  // 3D variable and 2D-only regridding is requested) and that the size of the
  // buffer (a_len) is calculated correctly.
  switch (dims) {
  case GRID_2D:
    count[Z] = 1;
    count[ZB] = 1;
    break;
  case GRID_3D:
    if (lic.no_regrid_ice) {
      SETERRQ(2, "no_regrid_ice is set, so dims == GRID_3D is not allowed");
    }
    count[ZB] = 1;
    break;
  case GRID_3D_BEDROCK:
    if (lic.no_regrid_bedrock) {
      SETERRQ(2, "no_regrid_bedrock is set, so dims == GRID_3D_BEDROCK is not allowed");
    }
    count[Z] = 1;
    break;
  }

  if (rank == 0) {

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
        MPI_Recv(start, N, MPI_INT, proc, start_tag, com, &mpi_stat);
	MPI_Recv(count, N, MPI_INT, proc, count_tag, com, &mpi_stat);
      }

      // Assemble nc_start and nc_count that are used below in the call to
      // nc_get_vara_double(...)
      size_t *nc_start, *nc_count;
      ptrdiff_t *imap;
      ierr = compute_start_and_count(varid, start, count, dims, nc_start, nc_count, imap); CHKERRQ(ierr);

      stat = nc_get_varm_double(ncid, varid, nc_start, nc_count, NULL, imap, lic.a);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      delete[] nc_start;
      delete[] nc_count;
      delete[] imap;

      // Find out how big the buffer actually is.  Remember that node 0 has a
      // buffer that will only be filled by if the process it is serving has a
      // maximal sized local domain.
      // Also note that at least one of count[Z] and count[ZB] is 1.
      int a_len = count[X] * count[Y] * count[Z] * count[ZB];

      // send the filled buffer
      if (proc != 0) {
        MPI_Send(lic.a, a_len, MPI_DOUBLE, proc, data_tag, com);
      }
    }
  } else { // not process 0:
    MPI_Send(start, N, MPI_INT, 0, start_tag, com);  // send out my start
    MPI_Send(count, N, MPI_INT, 0, count_tag, com);  // send out my count
    MPI_Recv(lic.a, lic.a_len, MPI_DOUBLE, 0, data_tag, com, &mpi_stat); // get back filled buffer
  }

  // At this point, the buffer lic.a[] should contain lic.a_len doubles.  This is the 
  // local processor's part of the source variable.  
  // That is, it should be enough of the source variable so that \e interpolation
  // (not extrapolation) onto the local processor's part of the target grid is possible.
  //ierr = lic.printArray(grid.com); CHKERRQ(ierr);
  
  // indexing parameters
  int myMz = 1, zcount = 1; // indexing trick so that we don't have to duplicate code for the 2-D case.
  switch (dims) {
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
                     z = (dims == GRID_3D) ? grid->zlevels[k] : grid->zblevels[k];

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

        if (dims == GRID_3D || dims == GRID_3D_BEDROCK) {
          // get the index into the source grid, for just below the level z
          const int kc = (dims == GRID_3D) ? lic.kBelowHeight(z)
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
          const double zkc = (dims == GRID_3D) ? lic.zlevs[kc] : lic.zblevs[kc];
          double dz;
          if (dims == GRID_3D) {
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
            SETERRQ(3,"PISMIO::myMaskInterp needed, but not initialized correctly");
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



//! Computes the size of the local block.
int PISMIO::compute_block_size(GridType dims, int* count) const {
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


//! Assembles start, count and imap arrays for a particular variable.
/*!
  Does nothing on processors other than zero.

  Arrays \c start and \c count have to have length 5, to be able to hold
  values for t,x,y,z,zb (in this order).

  Arrays nc_start and nc_count will be allocated *by* this function and have to
  be freed by the user.

  Also note that here X and Y have PISM (internal) meaning.

  Regarding imap: here's the description from section <a
  href="http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c/nc_005fget_005fvarm_005f-type.html">6.27
  Read a Mapped Array of Values</a> of the NetCDF C Interface Guide:

  A vector of integers that speciﬁes the mapping between the dimensions of a
  netCDF variable and the in-memory structure of the internal data array.
  imap[0] gives the distance between elements of the internal array
  corresponding to the most slowly varying dimension of the netCDF variable.
  imap[n-1] (where n is the rank of the netCDF variable) gives the distance
  between elements of the internal array corresponding to the most rapidly
  varying dimension of the netCDF variable. Intervening imap elements
  correspond to other dimensions of the netCDF variable in the obvious way.
  Distances between elements are speciﬁed in type-independent units of elements
*/
PetscErrorCode PISMIO::compute_start_and_count(const int varid, int *start, int *count, GridType dims,
					       size_t* &nc_start, size_t* &nc_count, ptrdiff_t* &imap) const {
  int stat, ndims;
  // Indices:
  const int T = 0, X = 1, Y = 2, Z = 3, ZB = 4;
  // IDs of interesting dimensions:
  int t_id, x_id, y_id, z_id, zb_id;
  // IDs of all the dimensions a variables depends on:
  int *dimids;

  // Quit if this is not processor zero:
  if (rank != 0) return 0;

//   for (int j = 0; j < 5; j++)
//     fprintf(stderr, "start[%d] = %d, count[%d] = %d\n", j, start[j], j, count[j]);

  // Get the number of dimensions a variable depends on:
  stat = nc_inq_varndims(ncid, varid, &ndims);
  CHKERRQ(check_err(stat,__LINE__,__FILE__));

  // Allocate all the arrays we need:
  nc_start = new size_t[ndims];
  nc_count = new size_t[ndims];
  imap     = new ptrdiff_t[ndims];
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

  size_t z_count;
  switch (dims) {
  case GRID_3D:
    {
      z_count = count[Z];
      break;
    }
  case GRID_3D_BEDROCK:
    {
      z_count = count[ZB];
      break;
    }
  default:
    z_count = 1;
  }

  // Assemble nc_start, nc_count and imap:
  for (int j = 0; j < ndims; j++) {
    if (dimids[j] == t_id) {
      nc_start[j] = start[T];
      nc_count[j] = count[T];
      imap[j]     = 1;		// this value does not matter because we never read more than 1 record
    } else if (dimids[j] == y_id) {
      nc_start[j] = start[Y];
      nc_count[j] = count[Y];
      imap[j]     = z_count;
    } else if (dimids[j] == x_id) {
      nc_start[j] = start[X];
      nc_count[j] = count[X];
      imap[j]     = count[Y] * z_count;
    } else if (dimids[j] == z_id) {
      nc_start[j] = start[Z];
      nc_count[j] = count[Z];
      imap[j]     = 1;
    } else if (dimids[j] == zb_id) {
      nc_start[j] = start[ZB];
      nc_count[j] = count[ZB];
      imap[j]     = 1;
    } else {
      nc_start[j] = 0;
      nc_count[j] = 1;
    }
    //    fprintf(stderr, "nc_start[%d] = %ld, nc_count[%d] = %ld, imap[%d] = %d\n", j, nc_start[j], j, nc_count[j], j, imap[j]); 
  }

  delete[] dimids;
  return 0;
}



//! Initializes the IceGrid object from a NetCDF file.
PetscErrorCode PISMIO::get_grid(const char filename[]) {
  PetscErrorCode ierr;
  grid_info gi;
  double *z_levels, *zb_levels;

  if (grid == NULL) SETERRQ(1, "PISMIO::get_grid(...): grid == NULL");

  ierr = open_for_reading(filename); CHKERRQ(ierr);

  ierr = get_grid_info(gi); CHKERRQ(ierr);
  ierr = get_vertical_dims(z_levels, zb_levels); CHKERRQ(ierr);

  grid->Mx = gi.x_len;
  grid->My = gi.y_len;
  grid->Lx = gi.Lx;
  grid->Ly = gi.Ly;
  // NB: We don't try to deduce periodicity from an input file.
  //  grid->periodicity = NONE; 
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


//! Open a NetCDF file for writing.
/*!
  If append == false moves an existing file aside if necessary.

  if check_dims == true, makes sure dimensions are OK.
 */
PetscErrorCode PISMIO::open_for_writing(string filename, bool append,
					bool check_dims) {
  int stat;

  if (ncid >= 0) SETERRQ(1, "PISMIO::open_for_writing(): ncid >= 0 at the beginning of the call");

  if (append == false) {
    // if append == false, we need to check if the file exists and move it
    // before proceeding if it does:
    if (rank == 0) {
      bool file_exists = false;

      // Check if the file exists:
      if (FILE *f = fopen(filename.c_str(), "r")) {
	file_exists = true;
	fclose(f);
      } else {
	file_exists = false;
      }
    
      if (file_exists && !append) {
	string tmp = filename + "~";
      
	stat = rename(filename.c_str(), tmp.c_str());
	if (stat != 0) {
	  stat = verbPrintf(1, com, "PISM ERROR: can't move '%s' to '%s'.\n",
			    filename.c_str(), tmp.c_str());
	  PetscEnd();
	}
	stat = verbPrintf(2, com, 
			  "PISM WARNING: output file '%s' already exists. Moving it to '%s'.\n",
			  filename.c_str(), tmp.c_str());
      }    
    } // end of if (rank == 0)
  }   // end of if (append == false)

  stat = NCTool::open_for_writing(filename.c_str()); CHKERRQ(stat);

  // If we don't need to check dimensions, we're done.
  if (!check_dims)
    return 0;

  if (append == true) {
    bool dimensions_are_ok = check_dimensions();

    if (!dimensions_are_ok) {
      stat = PetscPrintf(com,
			 "PISM ERROR: file '%s' has dimensions incompatible with the current grid. Exiting...\n",
			 filename.c_str()); CHKERRQ(stat);
      PetscEnd();
    }
  } else {
    // the file we just opened is empty, so we need to create dimensions
    stat = create_dimensions(); CHKERRQ(stat);
  }

  return 0;
}


//! Call this before calling regrid_local_var() or regrid_global_var() with argument useMaskInterp = true.
PetscErrorCode PISMIO::set_MaskInterp(MaskInterp *mi_in) {
  myMaskInterp = mi_in;
  return 0;
}


//! Always returns true on processors other than zero.
bool PISMIO::check_dimensions() const {
  bool t, x, y, z, zb;

  if (grid == NULL) SETERRQ(1, "PISMIO::check_dimensions(...): grid == NULL");

  t  = check_dimension("t", -1); // length does not matter
  x  = check_dimension("y", grid->My);
  y  = check_dimension("x", grid->Mx);
  z  = check_dimension("z", grid->Mz);
  zb = check_dimension("zb", grid->Mbz);
  
  return (t & x & y & z & zb);
}


//! Create dimensions and coordinate variables for storing spatial data.
/*! Assumes that the dataset is in the data mode. */
PetscErrorCode PISMIO::create_dimensions() const {
  int stat, t, x, y, z, zb, dimid;

  if (grid == NULL) SETERRQ(1, "PISMIO::create_dimensions(...): grid == NULL");

  if (rank == 0) {
    // define dimensions and coordinate variables:
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    // t
    stat = nc_def_dim(ncid, "t", NC_UNLIMITED, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_def_var(ncid, "t", NC_DOUBLE, 1, &dimid, &t); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_put_att_text(ncid, t, "long_name", 4, "time"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, t, "units", 17, "years since 1-1-1"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, t, "calendar", 4, "none"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, t, "axis", 1, "T"); check_err(stat,__LINE__,__FILE__);
    // x
    stat = nc_def_dim(ncid, "x", grid->Mx, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_def_var(ncid, "x", NC_DOUBLE, 1, &dimid, &x); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_put_att_text(ncid, x, "axis", 1, "X"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, x, "long_name", 32, "X-coordinate in Cartesian system"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, x, "standard_name", 23, "projection_x_coordinate"); check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, x, "units", 1, "m"); check_err(stat,__LINE__,__FILE__);
    // y
    stat = nc_def_dim(ncid, "y", grid->My, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
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
    
    double *x_coords, *y_coords;
    stat = grid->compute_horizontal_coordinates(x_coords, y_coords); CHKERRQ(stat);

    stat = put_dimension(y, grid->My, y_coords); CHKERRQ(stat);
    stat = put_dimension(x, grid->Mx, x_coords); CHKERRQ(stat);
    stat = put_dimension(z, grid->Mz, grid->zlevels); CHKERRQ(stat);
    stat = put_dimension(zb, grid->Mbz, grid->zblevels); CHKERRQ(stat);

    delete[] x_coords;
    delete[] y_coords;
  }

  return 0;
}

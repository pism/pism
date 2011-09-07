// Copyright (C) 2006--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "PISMIO.hh"
#include "pism_const.hh"

PISMIO::PISMIO(IceGrid *my_grid)
  : NCTool(my_grid->com, my_grid->rank) {
  grid = my_grid;

  event_write       = grid->profiler->create("pismio_write", "time spent in PISMIO::put_var()");
  event_write_proc0 = grid->profiler->create("pismio_write_proc0",
                                             "time spent writing data on processor 0");
  event_write_send_and_receive = grid->profiler->create("pismio_write_send_and_receive",
                                                        "time spent sending data to processor 0");
}

//! Read the first and last values, and the lengths, of the x,y,z dimensions from a NetCDF file. Read the last t.
PetscErrorCode PISMIO::get_grid_info(string name, grid_info &g) const {
  PetscErrorCode ierr;

  vector<int> dimids;
  int varid;
  bool exists;

  // try "name" as the standard_name first, then as the short name:
  ierr = find_variable(name, name, &varid, exists); CHKERRQ(ierr);

  if (!exists) {
    return 1;
  }

  ierr = inq_dimids(varid, dimids); CHKERRQ(ierr);

  int ndims = (int)dimids.size();
  for (int i = 0; i < ndims; ++i) {
    int dimid = dimids[i];

    string dimname;
    ierr = inq_dimname(dimid, dimname); CHKERRQ(ierr);

    AxisType dimtype = UNKNOWN_AXIS;
    ierr = inq_dimtype(dimname, dimtype); CHKERRQ(ierr);

    switch (dimtype) {
    case X_AXIS:
      {
        ierr = get_dim_length(dimname, &g.x_len); CHKERRQ(ierr);
        ierr = get_dim_limits(dimname, &g.x_min, &g.x_max); CHKERRQ(ierr);
        break;
      }
    case Y_AXIS:
      {
        ierr = get_dim_length(dimname, &g.y_len); CHKERRQ(ierr);
        ierr = get_dim_limits(dimname, &g.y_min, &g.y_max); CHKERRQ(ierr);
        break;
      }
    case Z_AXIS:
      {
        ierr = get_dim_length(dimname, &g.z_len); CHKERRQ(ierr);
        ierr = get_dim_limits(dimname, &g.z_min, &g.z_max); CHKERRQ(ierr);
        ierr = get_dimension(dimname, g.zlevels); CHKERRQ(ierr);
        break;
      }
    case T_AXIS:
      {
        ierr = get_dim_length(dimname,  &g.t_len); CHKERRQ(ierr);
        ierr = get_dim_limits(dimname, NULL, &g.time); CHKERRQ(ierr);
        break;
      }
    default:
      {
        PetscPrintf(grid->com, "ERROR: Can't figure out which direction dimension '%s' corresponds to.",
                    dimname.c_str());
        PISMEnd();
      }
    } // switch
  }   // for loop

  g.x0 = (g.x_max + g.x_min) / 2.0;
  g.y0 = (g.y_max + g.y_min) / 2.0;
  
  g.Lx = (g.x_max - g.x_min) / 2.0;
  g.Ly = (g.y_max - g.y_min) / 2.0;

  return 0;
}

//! Read a variable in a NetCDF file into a \c DA -managed global \c Vec \c g.  \e In \e parallel.
PetscErrorCode PISMIO::get_var(const int varid, Vec g, int z_count, int t) const {
  const int start_tag = 1;
  const int count_tag = 2;
  const int data_tag =  3;
  const int imap_tag =  4;
  PetscErrorCode ierr;
  MPI_Status mpi_stat;
  int stat;
  double *a_double = NULL;

  ierr = data_mode(); CHKERRQ(ierr); 

  vector<int> start, count, imap;
  ierr = compute_start_and_count(varid,
                                     t,
                                     grid->xs, grid->xm,
                                     grid->ys, grid->ym,
                                     0, z_count,
                                     start, count, imap); CHKERRQ(ierr);

  int N = (int)start.size();

  // Find the local size:
  int block_size, buffer_size;
  block_size = compute_block_size(count);
  // And the maximum size of the data block:
  buffer_size = block_size;
  ierr = MPI_Reduce(&block_size, &buffer_size, 1, MPI_INT, MPI_MAX, 0, com); CHKERRQ(ierr);

  // Memory allocation:
  ierr = PetscMalloc(buffer_size * sizeof(double), &a_double); CHKERRQ(ierr);

  if (rank == 0) {
    vector<int> start0(N), count0(N), imap0(N);
    for (int j = 0; j < N; j++) {
      // root needs to save its range
      start0[j] = start[j];
      count0[j] = count[j];
      imap0[j]  = imap[j];
    }
    vector<size_t> nc_start(N), nc_count(N);
    vector<ptrdiff_t> nc_imap(N);

    for (int proc = grid->size - 1; proc >= 0; proc--) { // root will read itself last

      if (proc == 0) {
        for (int j = 0; j < N; j++) {
	  start[j] = start0[j];
	  count[j] = count0[j];
          imap[j]  = imap0[j];
	}
      } else {
        MPI_Recv(&start[0], N, MPI_INT, proc, start_tag, com, &mpi_stat);
        MPI_Recv(&count[0], N, MPI_INT, proc, count_tag, com, &mpi_stat);
        MPI_Recv(&imap[0],  N, MPI_INT, proc, imap_tag,  com, &mpi_stat);
      }

      // Convert start, count and imap to types NetCDF requires:
      for (int i = 0; i < N; ++i) {
        nc_start[i] = start[i];
        nc_count[i] = count[i];
        nc_imap[i]  = imap[i];
      }

      stat = nc_get_varm_double(ncid, varid, &nc_start[0], &nc_count[0], NULL, &nc_imap[0], a_double);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      if (proc != 0) {
	block_size = compute_block_size(count);
	MPI_Send(a_double, block_size, MPI_DOUBLE, proc, data_tag, com);
      }
    }
  } else {                      // if (rank == 0)
    MPI_Send(&start[0], N, MPI_INT, 0, start_tag, com);
    MPI_Send(&count[0], N, MPI_INT, 0, count_tag, com);
    MPI_Send(&imap[0],  N, MPI_INT, 0, imap_tag,  com);
    MPI_Recv(a_double, buffer_size, MPI_DOUBLE, 0, data_tag, com, &mpi_stat);
  }

  block_size = compute_block_size(count);
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
PetscErrorCode PISMIO::put_var(const int varid, Vec g, int z_count) const {
  const int start_tag = 1;
  const int count_tag = 2;
  const int data_tag =  3;
  const int imap_tag =  4;
  PetscErrorCode ierr;
  MPI_Status mpi_stat;
  int stat;
  double *a_double = NULL;

  ierr = data_mode(); CHKERRQ(ierr); 

  grid->profiler->begin(event_write);

  // Fill start and count arrays:
  unsigned int t;
  ierr = get_dim_length("t", &t); CHKERRQ(ierr);
  t = t - 1;

  vector<int> start, count, imap;
  ierr = compute_start_and_count(varid,
                                     t,
                                     grid->xs, grid->xm,
                                     grid->ys, grid->ym,
                                     0, z_count,
                                     start, count, imap); CHKERRQ(ierr);

  int N = (int)start.size();

  // Find the local size:
  int block_size, buffer_size;
  block_size = compute_block_size(count);
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

  if (rank == 0) {
    // on rank 0 processor, receive messages from every other processor, then
    // write it out to the NC file
    vector<size_t> nc_start(N), nc_count(N);
    vector<ptrdiff_t> nc_imap(N);

    for (int proc = 0; proc < grid->size; proc++) { // root will write itself last
      if (proc != 0) {
        grid->profiler->begin(event_write_send_and_receive);
        MPI_Recv(&start[0], N, MPI_INT, proc, start_tag, com, &mpi_stat);
        MPI_Recv(&count[0], N, MPI_INT, proc, count_tag, com, &mpi_stat);
        MPI_Recv(&imap[0],  N, MPI_INT, proc, imap_tag,  com, &mpi_stat);
	MPI_Recv(a_double, buffer_size, MPI_DOUBLE, proc, data_tag, com, &mpi_stat);
        grid->profiler->end(event_write_send_and_receive);
      }
      
      grid->profiler->begin(event_write_proc0);

      // Convert start, count and imap to types NetCDF requires:

      // We initialize the 'stride' array instead of passing a NULL pointer to
      // nc_put_varm_double below to avoid a segfault caused by a bug in NetCDF
      // 4.1.2 and 4.1.3. Please see https://www.unidata.ucar.edu/jira/browse/NCF-88 for details.
      vector<ptrdiff_t> stride(N);

      for (int i = 0; i < N; ++i) {
        nc_start[i] = start[i];
        nc_count[i] = count[i];
        nc_imap[i]  = imap[i];
        stride[i]   = 1;
      }

      stat = nc_put_varm_double(ncid, varid, &nc_start[0], &nc_count[0],
                                &stride[0], &nc_imap[0], a_double);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      grid->profiler->end(event_write_proc0);
    }
  } else {  // all other processors send to rank 0 processor
    grid->profiler->begin(event_write_send_and_receive);
    MPI_Send(&start[0], N, MPI_INT, 0, start_tag, com);
    MPI_Send(&count[0], N, MPI_INT, 0, count_tag, com);
    MPI_Send(&imap[0],  N, MPI_INT, 0, imap_tag,  com);
    MPI_Send(a_double, buffer_size, MPI_DOUBLE, 0, data_tag, com);
    grid->profiler->end(event_write_send_and_receive);
  }

  // Cleanup:
  ierr = PetscFree(a_double); CHKERRQ(ierr);

  grid->profiler->end(event_write);
  
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
PetscErrorCode PISMIO::regrid_var(const int varid, const vector<double> &zlevels_out, LocalInterpCtx *lic,
				  Vec g) const {
  PetscErrorCode ierr;
  const int T = 0, X = 1, Y = 2, Z = 3; // indices, just for clarity
  const int start_tag = 1; // MPI tag for the start array
  const int count_tag = 2; // MPI tag for the count array
  const int data_tag  = 3; // MPI tag for the data block
  const int imap_tag  = 4;
  MPI_Status mpi_stat;
  int stat;
  vector<int> start, count, imap;

  ierr = data_mode(); CHKERRQ(ierr);

  int t_start = lic->start[T],
    x_start = lic->start[X],
    y_start = lic->start[Y],
    x_count = lic->count[X],
    y_count = lic->count[Y],
    z_start = lic->start[Z],
    z_count = lic->count[Z];

  // FIXME: this will need to go
  vector<double> &zlevels_in = lic->zlevels;
  double *buffer = lic->a;
  int buffer_length = lic->a_len;

  int nlevels = (int)zlevels_out.size();

  ierr = compute_start_and_count(varid, t_start,
                                     x_start, x_count, y_start, y_count, z_start, z_count,
                                     start, count, imap); CHKERRQ(ierr);

  int N = (int)start.size();

  if (rank == 0) {

    // Node 0 will service all the other nodes before itself.  We need to save
    // start[] and count[] so that it knows how to get its block at the end.
    vector<int> start0(N), count0(N), imap0(N);

    for (int i = 0; i < N; i++) {
      start0[i] = start[i];
      count0[i] = count[i];
      imap0[i]  = imap[i];
    }

    vector<size_t> nc_start(N), nc_count(N);
    vector<ptrdiff_t> nc_imap(N);

    for (int proc = grid->size - 1; proc >= 0; proc--) {
      if (proc == 0) {// Get the bounds.
	for (int i = 0; i < N; i++) {
	  start[i] = start0[i];
	  count[i] = count0[i];
          imap[i]  = imap0[i];
	}
      } else {
        MPI_Recv(&start[0], N, MPI_INT, proc, start_tag, com, &mpi_stat);
	MPI_Recv(&count[0], N, MPI_INT, proc, count_tag, com, &mpi_stat);
	MPI_Recv(&imap[0],  N, MPI_INT, proc, imap_tag,  com, &mpi_stat);
      }

      // Convert start, count and imap to types NetCDF requires and compute how
      // much data we're going to read (for communication, below):
      int a_len = 1;
      for (int i = 0; i < N; ++i) {
        nc_start[i] = start[i];
        nc_count[i] = count[i];
        nc_imap[i]  = imap[i];
        a_len = a_len * count[i];
      }

      stat = nc_get_varm_double(ncid, varid, &nc_start[0], &nc_count[0], NULL, &nc_imap[0], buffer);
      CHKERRQ(check_err(stat,__LINE__,__FILE__));

      // send the filled buffer
      if (proc != 0) {
        MPI_Send(buffer, a_len, MPI_DOUBLE, proc, data_tag, com);
      }
    }
  } else { // not process 0:
    MPI_Send(&start[0], N, MPI_INT, 0, start_tag, com);  // send out my start
    MPI_Send(&count[0], N, MPI_INT, 0, count_tag, com);  // send out my count
    MPI_Send(&imap[0],  N, MPI_INT, 0, imap_tag,  com);  // send out my imap
    MPI_Recv(buffer, buffer_length, MPI_DOUBLE, 0, data_tag, com, &mpi_stat); // get back filled buffer
  }

  // At this point, the buffer buffer[] should contain buffer_length doubles.  This is the 
  // local processor's part of the source variable.  
  // That is, it should be enough of the source variable so that \e interpolation
  // (not extrapolation) onto the local processor's part of the target grid is possible.
  //ierr = lic->printArray(grid.com); CHKERRQ(ierr);
  
  // We'll work with the raw storage here so that the array we are filling is
  // indexed the same way as the buffer we are pulling from (buffer)
  PetscScalar *vec_a;
  ierr = VecGetArray(g, &vec_a); CHKERRQ(ierr);

  for (int i = grid->xs; i < grid->xs + grid->xm; i++) {
    for (int j = grid->ys; j < grid->ys + grid->ym; j++) {

      for (int k = 0; k < nlevels; k++) {
        // location (x,y,z) is in target computational domain
        const double
          x = grid->x[i] - grid->x0,
          y = grid->y[j] - grid->y0,
          z = zlevels_out[k];

        // We need to know how the point (x,y,z) sits within the local block we
        // pulled from the netCDF file.  This part is special to a regular
        // grid.  In particular floor(ic) is the index of the 'left neighbor'
        // and ceil(ic) is the index of the 'right neighbor'.
        double ic = (x - lic->fstart[X]) / lic->delta[X];
        double jc = (y - lic->fstart[Y]) / lic->delta[Y];

        // bounds checking on ic and jc; cases where this is essential *have* been observed
        if ((int)floor(ic) < 0) {
          ic = 0.0;
          //SETERRQ2(101,"(int)floor(ic) < 0      [%d < 0; ic = %16.15f]",(int)floor(ic),ic);
        }
        if ((int)floor(jc) < 0) {
          jc = 0.0;
          //SETERRQ2(102,"(int)floor(jc) < 0      [%d < 0; jc = %16.15f]",(int)floor(jc),jc);
        }
        if ((int)ceil(ic) > x_count-1) {
          ic = (double)(x_count-1);
          //SETERRQ3(103,"(int)ceil(ic) > lic->count[1]-1      [%d > %d; ic = %16.15f]",
          //     (int)ceil(ic), lic->count[1]-1, ic);
        }
        if ((int)ceil(jc) > y_count-1) {
          jc = (double)(y_count-1);
          //SETERRQ3(104,"(int)ceil(jc) > lic->count[2]-1      [%d > %d; jc = %16.15f]",
          //     (int)ceil(jc), lic->count[2]-1, jc);
        }

        double a_mm, a_mp, a_pm, a_pp;  // filled differently in 2d and 3d cases

        const int Im = (int)floor(ic), Ip = (int)ceil(ic),
          Jm = (int)floor(jc), Jp = (int)ceil(jc);

        if (nlevels > 1) {
          // get the index into the source grid, for just below the level z
          const int kc = k_below(z, zlevels_in);

          // We pretend that there are always 8 neighbors (4 in the map plane,
          // 2 vertical levels). And compute the indices into the buffer for
          // those neighbors.
          // Note that floor(ic) + 1 = ceil(ic) does not hold when ic is an
          // integer.  Computation of the domain (in constructor of LocalInterpCtx; note
          // that lic->count uses ceil) must be done in a compatible way,
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
          const int mmm = (Im * y_count + Jm) * z_count + kc;
          const int mmp = (Im * y_count + Jm) * z_count + kc + 1;
          const int mpm = (Im * y_count + Jp) * z_count + kc;
          const int mpp = (Im * y_count + Jp) * z_count + kc + 1;
          const int pmm = (Ip * y_count + Jm) * z_count + kc;
          const int pmp = (Ip * y_count + Jm) * z_count + kc + 1;
          const int ppm = (Ip * y_count + Jp) * z_count + kc;
          const int ppp = (Ip * y_count + Jp) * z_count + kc + 1;

          // We know how to index the neighbors, but we don't yet know where the
          // point lies within this box.  This is represented by kk in [0,1].
          // For the irregular case, with left index km and right index kp in the source grid, we
          // would have
          //   kk = (km == kp) ? 0.0 : (z - Z(km)) / (Z(kp) - Z(km))
          // where Z(.) are the physical coordinates on the source grid.  Note
          // that any value in [0,1] would be okay when km == kp.
          const double zkc = zlevels_in[kc];
          double dz;
          if (kc == z_count - 1) {
            dz = zlevels_in[kc] - zlevels_in[kc-1];
          } else {
            dz = zlevels_in[kc+1] - zlevels_in[kc];
          }
          const double kk = (z - zkc) / dz;

          // linear interpolation in the z-direction
          a_mm = buffer[mmm] * (1.0 - kk) + buffer[mmp] * kk;
          a_mp = buffer[mpm] * (1.0 - kk) + buffer[mpp] * kk;
          a_pm = buffer[pmm] * (1.0 - kk) + buffer[pmp] * kk;
          a_pp = buffer[ppm] * (1.0 - kk) + buffer[ppp] * kk;
        } else {
          // we don't need to interpolate vertically for the 2-D case
          a_mm = buffer[Im * y_count + Jm];
          a_mp = buffer[Im * y_count + Jp];
          a_pm = buffer[Ip * y_count + Jm];
          a_pp = buffer[Ip * y_count + Jp];
        }

        const double jj = jc - floor(jc);

        // interpolate in y direction
        const double a_m = a_mm * (1.0 - jj) + a_mp * jj;
        const double a_p = a_pm * (1.0 - jj) + a_pp * jj;

        const double ii = ic - floor(ic);
        int index = ((i - grid->xs) * grid->ym + (j - grid->ys)) * nlevels + k;

        // index into the new array and interpolate in x direction
        vec_a[index] = a_m * (1.0 - ii) + a_p * ii;
        // done with the point at (x,y,z)
      }
    }
  }

  ierr = VecRestoreArray(g, &vec_a); CHKERRQ(ierr);

  return 0;
}

int PISMIO::k_below(double z, const vector<double> &zlevels) const {
  double z_min = zlevels.front(), z_max = zlevels.back();
  PetscInt mcurr = 0;

  if (z < z_min - 1.0e-6 || z > z_max + 1.0e-6) {
    PetscPrintf(com, 
                "PISMIO::k_below(): z = %5.4f is outside the allowed range.\n", z);
    PISMEnd();
  }

  while (zlevels[mcurr+1] < z)
    mcurr++;

  return mcurr;
}


//! Computes the size of the local block.
int PISMIO::compute_block_size(vector<int> count) const {
  int N = (int)count.size(), result = 1;

  for (int i = 0; i < N; ++i)
    result = result * count[i];

  return result;
}




//! Initializes the IceGrid object from a NetCDF file.
PetscErrorCode PISMIO::get_grid(string filename, string var_name) {
  PetscErrorCode ierr;
  grid_info gi;
  vector<double> z_levels, zb_levels;

  if (grid == NULL) SETERRQ(1, "PISMIO::get_grid(...): grid == NULL");

  ierr = open_for_reading(filename); CHKERRQ(ierr);

  ierr = get_grid_info(var_name, gi);
  // Close the file and return 1 if the variable could not be found. We don't
  // use CHKERRQ(ierr) to let the caller handle this error.
  if (ierr != 0) {
    close();
    return 1;
  }

  // if we have no vertical grid information, create a fake 2-level vertical grid.
  if (gi.zlevels.size() < 2) {
    double Lz = grid->config.get("grid_Lz");
    ierr = verbPrintf(3, com,
                      "WARNING: Can't determine vertical grid information using '%s' in %s'\n"
                      "         Using 2 levels and Lz of %3.3fm\n",
                      var_name.c_str(), filename.c_str(), Lz); CHKERRQ(ierr);

    gi.zlevels.push_back(0);
    gi.zlevels.push_back(Lz);
  }

  grid->Mx = gi.x_len;
  grid->My = gi.y_len;
  grid->Lx = gi.Lx;
  grid->Ly = gi.Ly;
  // NB: We don't try to deduce periodicity from an input file.
  //  grid->periodicity = NONE; 
  grid->x0   = gi.x0;
  grid->y0   = gi.y0;
  grid->year = convert(gi.time, "seconds", "years");

  ierr = grid->compute_horizontal_spacing(); CHKERRQ(ierr);
  ierr = grid->set_vertical_levels(gi.zlevels); CHKERRQ(ierr);

  // We're ready to call grid->createDA().

  // Cleanup:
  ierr = close(); CHKERRQ(ierr);
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
	  PISMEnd();
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
      PISMEnd();
    }
  } else {
    // the file we just opened is empty, so we need to create dimensions
    stat = create_dimensions(); CHKERRQ(stat);
  }

  return 0;
}

//! Always returns true on processors other than zero.
bool PISMIO::check_dimensions() const {

  if (grid == NULL) SETERRQ(1, "PISMIO::check_dimensions(...): grid == NULL");

  return check_dimension("t", -1); // length does not matter
}


//! Create dimensions and coordinate variables for storing spatial data.
/*! Assumes that the dataset is in the data mode. */
PetscErrorCode PISMIO::create_dimensions() const {
  int ierr, dimid, varid;
  map<string,string> attrs;

  if (grid == NULL) SETERRQ(1, "PISMIO::create_dimensions(...): grid == NULL");
  
  if (grid->rank != 0) return 0;

  // define dimensions and coordinate variables:

  // t
  attrs["long_name"] = "time";
  attrs["calendar"]  = "365_day";
  attrs["units"]     = "years since " + grid->config.get_string("reference_date");
  attrs["axis"]      = "T";
  ierr = create_dimension("t", NC_UNLIMITED, attrs, dimid, varid); CHKERRQ(ierr);

  return 0;
}

//! Assembles start, count and imap arrays for a particular variable.
/*!
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

  The default I/O code (sending data to/from processor 0) only calls this on
  processor 0.
*/
PetscErrorCode PISMIO::compute_start_and_count(int varid,
                                               int t_start,
                                               int x_start, int x_count,
                                               int y_start, int y_count,
                                               int z_start, int z_count,
                                               vector<int> &start,
                                               vector<int> &count,
                                               vector<int> &imap) const {
  PetscErrorCode ierr;
  vector<int> dimids;

  ierr = inq_dimids(varid, dimids); CHKERRQ(ierr);
  int ndims = (int)dimids.size();

  // Resize output vectors:
  start.resize(ndims);
  count.resize(ndims);
  imap.resize(ndims);

  // Assemble start, count and imap:
  for (int j = 0; j < ndims; j++) {
    string dimname;
    ierr = inq_dimname(dimids[j], dimname); CHKERRQ(ierr);
    
    AxisType dimtype;
    ierr = inq_dimtype(dimname, dimtype); CHKERRQ(ierr);

    switch (dimtype) {
    case T_AXIS:
      start[j] = t_start;
      count[j] = 1;             // t_count is always 1
      imap[j]  = 1; // this value does not matter because we never read more than 1 record
      break;
    case X_AXIS:
      start[j] = x_start;
      count[j] = x_count;
      imap[j]  = y_count * z_count;
      break;
    case Y_AXIS:
      start[j] = y_start;
      count[j] = y_count;
      imap[j]  = z_count;
      break;
    case Z_AXIS:
      start[j] = z_start;
      count[j] = z_count;
      imap[j]  = 1;
      break;
    default:
      start[j] = 0;
      count[j] = 1;
      imap[j]  = 1;             // is this right?
    }

    //    fprintf(stderr, "start[%d] = %ld, count[%d] = %ld, imap[%d] = %d\n", j, start[j], j, count[j], j, imap[j]); 
  }

  return 0;
}


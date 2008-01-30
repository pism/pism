// Copyright (C) 2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __nc_util_hh
#define __nc_util_hh

#include <netcdf.h>
#include <petscmat.h>
#include "grid.hh"


//! The "local interpolation context" describes the processor's part of the source NetCDF file (for regridding).
/*!
The local interpolation context contains the details of how the processor's block
of the new computational domain fits into the domain of the netCDF file.  Note that each vertical column 
of the grid is owned by exactly one processor.

For any particular dimension, we have a new computational domain \f$[a,b]\f$ with
spacing \f$h\f$ so there are \f$n = (b - a) / h\f$ interior cells, indexed by \f$\{i_0, \dots, i_n\}\f$.
The local processor owns a range \f$\{i_m, \dots, i_{m'}\}\f$.  Suppose the netCDF file has
domain \f$[A,B]\f$, spacing \f$H\f$, and \f$N = (B - A) / H\f$ cells.  In order to interpolate 
onto these points, we need the indices \f$\{I_m, \dots, I_{m'}\}\f$ of the netCDF file so that

  \f[  [x(i_m), x(i_{m'})] \quad \text{is a subset of} \quad  [x(I_m), x(I_{m'})]  \f]

In the code, \c xbdy\f$[2] = \{x(i_m), x(i_{m'})\}\f$.  We have obtained the netCDF bounds 
\f$x(I_0)\f$ and \f$x(I_N)\f$ in the \c bdy array and the
number of elements \f$N+1\f$ in the \c dim array.

This is a summary of the fields of a \c LocalInterpCtx, and a translation between scalars in the previous
paragraphs and the fields:
 - \f$I_m = \operatorname{floor}((x(i_m) - A) / H)\f$ is called \c start
 - \f$I_{m'} - I_m = \operatorname{ceil}((x(i_{m'}) - X(I_m)) / H\f$ is \c count - 1
 - \f$X(I_m)\f$ is called \c fstart
 - \f$H\f$ is called \c delta

\c delta and \c fstart are not used in the vertical dimension because the spacing is not 
generally equal.  Rather, \c zlevs and \c zblevs contain the needed information.

The arrays \c start and \c count have 5 integer entries, corresponding to the dimensions
\f$t, x, y, z, zb\f$.
 */
class LocalInterpCtx {
public:
//  double fstart[3], delta[4];
  double fstart[3], delta[3];
  int start[5], count[5];    // Indices in netCDF file.
  float *a;
  int a_len;
  int ncid;
  int nz, nzb;
  double *zlevs, *zblevs;

public:
  LocalInterpCtx(int ncid, const size_t dim[], const double bdy[],
                 const double zlevsIN[], const double zblevsIN[], IceGrid &grid);
  ~LocalInterpCtx();
  int kBelowHeight(const double height, MPI_Comm com);
  int kbBelowHeight(const double elevation, MPI_Comm com);
  PetscErrorCode printGrid(MPI_Comm com);
  PetscErrorCode printArray(MPI_Comm com);
};


struct MaskInterp {
  int number_allowed;
  int allowed_levels[50]; // must be strictly increasing
};

int nc_check(int stat);
int check_err(const int stat, const int line, const char *file);


//! Collects together parallel NetCDF methods used by IceModel and IceModelVec.
class NCTool {

public:
NCTool();

PetscErrorCode put_dimension(int ncid, int v_id, int len, PetscScalar *vals);
PetscErrorCode put_dimension_regular(int ncid, int v_id, int len, double start, double delta);

PetscErrorCode get_dims_limits_lengths(int ncid, size_t dim[], double bdy[], MPI_Comm com);
PetscErrorCode get_ends_1d_var(int ncid, int vid, PetscScalar *gfirst, PetscScalar *glast, MPI_Comm com);

PetscErrorCode get_vertical_dims(int ncid, int z_len, int zb_len, 
                                 double z_read[], double zb_read[], MPI_Comm com);

PetscErrorCode put_local_var(const IceGrid *grid, int ncid, const int var_id, nc_type type,
                             DA da, Vec v, Vec g, const int *s, const int *c,
                             int dims, void *a_mpi, int a_size);
PetscErrorCode put_global_var(const IceGrid *grid, int ncid, const int var_id, nc_type type,
                              DA da, Vec g, const int *s, const int *c,
                              int dims, void *a_mpi, int a_size);

PetscErrorCode get_local_var(const IceGrid *grid, int ncid, const char *name, nc_type type,
                             DA da, Vec v, Vec g, const int *s, const int *c,
                             int dims, void *a_mpi, int a_size);
PetscErrorCode get_global_var(const IceGrid *grid, int ncid, const char *name, nc_type type,
                              DA da, Vec g, const int *s, const int *c,
                              int dims, void *a_mpi, int a_size);

PetscErrorCode var_to_da_vec(IceGrid &grid, int ncid, int vid, DA da, Vec vecl,
                             Vec vecg, Vec vindzero);
PetscErrorCode var_to_da_vec(IceGrid &grid, int ncid, int vid, DA da, Vec vecl,
                             Vec vecg, Vec vindzero, MaskInterp masktool);

PetscErrorCode regrid_local_var(const char *vars, char c, const char *name, int dim_flag,
                                LocalInterpCtx &lic, IceGrid &grid, DA da, Vec vec, Vec g);
PetscErrorCode regrid_global_var(const char *vars, char c, const char *name, int dim_flag,
                                 LocalInterpCtx &lic, IceGrid &grid, DA da, Vec g);
PetscErrorCode form_LocalInterpCtx(int ncid, const size_t dim[], const double bdy[],
                                   const double zlev[], const double zblev[],
                                   LocalInterpCtx &lic, IceGrid &grid);
};

#endif // __nc_util_hh


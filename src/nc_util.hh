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

struct LocalInterpCtx {
  float fstart[3], delta[4];
  double timestart;
  int start[5], count[5];    // Indices in netCDF file.
  float *a;
  int a_len;
  int ncid;
};

int check_err(const int stat, const int line, const char *file);
PetscErrorCode put_local_var(const IceGrid *grid, int ncid, const int var_id, nc_type type,
                             DA da, Vec v, Vec g, const int *s, const int *c,
                             int dims, void *a_mpi, int a_size);
PetscErrorCode put_dimension_regular(int ncid, int v_id, int len, float start, float delta);
PetscErrorCode get_dimensions(int ncid, size_t dim[], float bdy[], double *bdy_time, MPI_Comm com);
PetscErrorCode get_local_var(const IceGrid *grid, int ncid, const char *name, nc_type type,
                             DA da, Vec v, Vec g, const int *s, const int *c,
                             int dims, void *a_mpi, int a_size);
PetscErrorCode regrid_local_var(const char *vars, char c, const char *name, int dim_flag,
                                LocalInterpCtx &lic, IceGrid &grid, DA da, Vec vec, Vec g);
PetscErrorCode get_LocalInterpCtx(int ncid, const size_t dim[], const float bdy[], const double bdy_time,
                                  LocalInterpCtx &lic, IceGrid &grid);

#endif // __nc_util_hh

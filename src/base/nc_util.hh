// Copyright (C) 2007--2009 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "LocalInterpCtx.hh"
#include "grid.hh"
#include "../udunits/udunits.h"

struct MaskInterp {
  int number_allowed;
  int allowed_levels[50]; // must be strictly increasing
};

typedef enum {GRID_2D = 2, GRID_3D, GRID_3D_BEDROCK} GridType;

int nc_check(int stat);
int check_err(const int stat, const int line, const char *file);

//! Collects together parallel NetCDF methods used by IceModel and IceModelVec.
class NCTool {
public:
  int ncid;

public:
  NCTool(IceGrid *my_grid);
  PetscErrorCode open_for_reading(const char filename[], bool &exists);
  PetscErrorCode open_for_writing(const char filename[], bool replace = true);
  PetscErrorCode close();
  PetscErrorCode find_variable(const char short_name[], const char standard_name[],
			       int *varid, bool &exists);
  PetscErrorCode find_dimension(const char short_name[], int *dimid, bool &exists);
  PetscErrorCode create_dimensions();
  PetscErrorCode append_time(PetscReal time);
  PetscErrorCode write_global_attrs(bool have_ssa_velocities, const char conventions[]);
  PetscErrorCode write_history(const char history[], bool overwrite = false);
  PetscErrorCode get_last_time(double *time);
  PetscErrorCode get_dim_length(const char name[], int *len);
  PetscErrorCode get_dim_limits(const char name[], double *min, double *max);
  PetscErrorCode get_grid(const char filename[]);
  PetscErrorCode get_grid_info(grid_info &g);
  PetscErrorCode get_grid_info_2d(grid_info &g);
  PetscErrorCode get_vertical_dims(double* &z_levels, double* &zb_levels);

  PetscErrorCode put_dimension(int varid, int len, PetscScalar *vals);
  PetscErrorCode put_dimension_regular(int varid, int len, double start, double delta);

  PetscErrorCode read_polar_stereographic(double &straight_vertical_longitude_from_pole,
					  double &latitude_of_projection_origin,
					  double &standard_parallel,
					  bool report = false);
  PetscErrorCode get_att_text(int varid, const char name[], int *length, char **result);
  PetscErrorCode get_att_double(int varid, const char name[],
				int length, double *result);
  PetscErrorCode get_units(int varid, bool &has_units, utUnit &units);
  PetscErrorCode write_polar_stereographic(double straight_vertical_longitude_from_pole,
					   double latitude_of_projection_origin,
					   double standard_parallel);

  PetscErrorCode get_local_var(int varid, DA da, Vec v, GridType dims, int t);
  PetscErrorCode get_global_var(int varid, Vec g, GridType dims, int t);

  PetscErrorCode put_local_var(int varid, DA da, Vec v, GridType dims);
  PetscErrorCode put_global_var(int varid, Vec g, GridType dims);

  PetscErrorCode set_MaskInterp(MaskInterp *mi_in);
  PetscErrorCode regrid_local_var(int varid, GridType dim_flag,
				  LocalInterpCtx &lic, DA da, Vec vec,
				  bool useMaskInterp);
  PetscErrorCode regrid_global_var(int varid, GridType dim_flag,
				   LocalInterpCtx &lic, Vec g,
				   bool useMaskInterp);

private:
  int compute_block_size(GridType dims, int* count);
  PetscErrorCode compute_start_and_count(int varid, int *pism_start, int *pism_count,
					 size_t* &nc_start, size_t* &nc_count);
  PetscErrorCode transpose(int varid, GridType dim_flag, int *count,
			   double* in, double* out);
  bool check_dimension(const char dim[], int len);
  bool check_dimensions();
  MaskInterp  *myMaskInterp;
  IceGrid* grid;
};

#endif // __nc_util_hh


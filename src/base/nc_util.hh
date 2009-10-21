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

#include <string>
#include <vector>
#include <netcdf.h>
#include <petscmat.h>
#include "LocalInterpCtx.hh"
#include "grid.hh"
#include "../udunits/udunits.h"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

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
  NCTool(IceGrid *my_grid);
  NCTool(MPI_Comm c, PetscMPIInt r);
  ~NCTool();
  PetscErrorCode open_for_reading(const char filename[]);
  PetscErrorCode open_for_writing(const char filename[], bool append,
				  bool check_dims = false);
  PetscErrorCode set_MaskInterp(MaskInterp *mi_in);
  PetscErrorCode get_grid(const char filename[]);
  PetscErrorCode close();

  int get_ncid() const;
  PetscErrorCode find_variable(string short_name, string standard_name,
			       int *varid, bool &exists, bool &found_by_standard_name) const;
  PetscErrorCode find_variable(string short_name, string standard_name,
			       int *varid, bool &exists) const;
  PetscErrorCode find_variable(string short_name, int *varid, bool &exists) const;
  PetscErrorCode find_dimension(const char short_name[], int *dimid, bool &exists) const;
  PetscErrorCode create_dimensions() const;
  PetscErrorCode append_time(PetscReal time) const;
  PetscErrorCode write_history(const char history[], bool overwrite = false) const;
  PetscErrorCode get_last_time(double *time) const;
  PetscErrorCode get_dim_length(const char name[], int *len) const;
  PetscErrorCode get_dim_limits(const char name[], double *min, double *max) const;
  PetscErrorCode get_grid_info(grid_info &g) const;
  PetscErrorCode get_grid_info_2d(grid_info &g) const;
  PetscErrorCode get_vertical_dims(double* &z_levels, double* &zb_levels) const;

  PetscErrorCode put_dimension(int varid, int len, PetscScalar *vals) const;
  PetscErrorCode put_dimension_regular(int varid, int len, double start, double delta) const;

  PetscErrorCode inq_unlimdim(int &unlimdimid) const;
  PetscErrorCode inq_dimname(int dimid, string &name) const;
  PetscErrorCode inq_dimids(int varid, vector<int> &dimids) const;
  PetscErrorCode inq_nattrs(int varid, int &N) const;
  PetscErrorCode inq_att_name(int varid, int n, string &name) const;
  PetscErrorCode inq_att_type(int varid, const char name[], nc_type &typep) const;
  PetscErrorCode get_att_text(int varid, const char name[], string &result) const;
  PetscErrorCode get_att_double(int varid, const char name[], vector<double> &result) const;
  PetscErrorCode get_units(int varid, bool &has_units, utUnit &units) const;

  PetscErrorCode create_timeseries(const char name[], const char long_name[],
				   const char units[], nc_type nctype, int *varid) const;
  PetscErrorCode append_timeseries(const char name[], double value) const;

  PetscErrorCode get_local_var(int varid, DA da, Vec v, GridType dims, int t) const;
  PetscErrorCode get_global_var(int varid, Vec g, GridType dims, int t) const;

  PetscErrorCode put_local_var(int varid, DA da, Vec v, GridType dims) const;
  PetscErrorCode put_global_var(int varid, Vec g, GridType dims) const;

  PetscErrorCode regrid_local_var(int varid, GridType dim_flag,
				  LocalInterpCtx &lic, DA da, Vec vec,
				  bool useMaskInterp) const;
  PetscErrorCode regrid_global_var(int varid, GridType dim_flag,
				   LocalInterpCtx &lic, Vec g,
				   bool useMaskInterp) const;

private:
  PetscErrorCode open_for_writing(const char filename[]);
  int compute_block_size(GridType dims, int* count) const;
  PetscErrorCode compute_start_and_count(int varid, int *pism_start, int *pism_count, GridType dims,
					 size_t* &nc_start, size_t* &nc_count, ptrdiff_t* &imap) const;
  bool check_dimension(const char dim[], int len) const;
  bool check_dimensions() const;
  int ncid;
  MaskInterp  *myMaskInterp;
  IceGrid* grid;
  MPI_Comm com;
  PetscMPIInt rank;
};

#endif // __nc_util_hh


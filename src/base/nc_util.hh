// Copyright (C) 2007--2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petsc.h>
#include <netcdf.h>		// nc_type
#include "../udunits/udunits.h"	// utUnit
#include <string>
#include <vector>

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
  NCTool(MPI_Comm c, PetscMPIInt r);
  virtual ~NCTool();
  virtual PetscErrorCode open_for_reading(const char filename[]);
  virtual PetscErrorCode open_for_writing(const char filename[]);
  virtual PetscErrorCode close();

  virtual PetscErrorCode find_variable(string short_name, string standard_name,
			       int *varid, bool &exists, bool &found_by_standard_name) const;
  virtual PetscErrorCode find_variable(string short_name, string standard_name,
			       int *varid, bool &exists) const;
  virtual PetscErrorCode find_variable(string short_name, int *varid, bool &exists) const;
  virtual PetscErrorCode find_dimension(const char short_name[], int *dimid, bool &exists) const;
  virtual PetscErrorCode append_time(PetscReal time) const;
  virtual PetscErrorCode write_history(const char history[], bool overwrite = false) const;
  virtual PetscErrorCode get_last_time(double *time) const;
  virtual PetscErrorCode get_vertical_dims(double* &z_levels, double* &zb_levels) const;
  virtual bool check_dimension(const char dim[], int len) const;

  virtual PetscErrorCode get_dim_length(const char name[], int *len) const;
  virtual PetscErrorCode get_dim_limits(const char name[], double *min, double *max) const;

  virtual PetscErrorCode put_dimension(int varid, int len, PetscScalar *vals) const;
  virtual PetscErrorCode put_dimension_regular(int varid, int len, double start, double delta) const;

  virtual PetscErrorCode inq_unlimdim(int &unlimdimid) const;
  virtual PetscErrorCode inq_dimname(int dimid, string &name) const;
  virtual PetscErrorCode inq_dimids(int varid, vector<int> &dimids) const;
  virtual PetscErrorCode inq_nattrs(int varid, int &N) const;
  virtual PetscErrorCode inq_att_name(int varid, int n, string &name) const;
  virtual PetscErrorCode inq_att_type(int varid, const char name[], nc_type &typep) const;
  virtual PetscErrorCode get_att_text(int varid, const char name[], string &result) const;
  virtual PetscErrorCode get_att_double(int varid, const char name[], vector<double> &result) const;
  virtual PetscErrorCode get_units(int varid, bool &has_units, utUnit &units) const;

  virtual int get_ncid() const;
protected:
  int ncid;
  MPI_Comm com;
  PetscMPIInt rank;
};

#endif // __nc_util_hh


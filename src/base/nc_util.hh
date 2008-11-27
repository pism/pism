// Copyright (C) 2007, 2008 Jed Brown, Ed Bueler and Constantine Khroulev
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
  double fstart[3], delta[3];
  int start[5], count[5];    // Indices in netCDF file.
  double *a;
  int a_len;
  int nz, nzb;
  double *zlevs, *zblevs;
  bool regrid_2d_only, no_regrid_bedrock;

public:
  LocalInterpCtx(const size_t dim[], const double bdy[],
                 const double zlevsIN[], const double zblevsIN[], IceGrid &grid);
  ~LocalInterpCtx();
  int kBelowHeight(const double height, MPI_Comm com);
  int kbBelowHeight(const double elevation, MPI_Comm com);
  PetscErrorCode printGrid(MPI_Comm com);
  PetscErrorCode printArray(MPI_Comm com);
protected:
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
  int ncid;

public:
  NCTool(IceGrid *my_grid);
  NCTool();
  PetscErrorCode set_grid(IceGrid *my_grid);
  PetscErrorCode open_for_reading(const char filename[], bool &exists);
  PetscErrorCode open_for_writing(const char filename[]);
  PetscErrorCode close();
  PetscErrorCode find_variable(const char short_name[], const char standard_name[],
			       int *varid, bool &exists);
  PetscErrorCode find_dimension(const char short_name[], int *dimid, bool &exists);
  PetscErrorCode create_dimensions();
  PetscErrorCode append_time(PetscReal time);
  PetscErrorCode write_global_attrs(bool have_ssa_velocities, const char conventions[]);
  PetscErrorCode write_history(const char history[]);
  PetscErrorCode get_last_time(double *time);
  PetscErrorCode get_dims_limits_lengths(size_t dim[], double bdy[]);
  PetscErrorCode get_dims_limits_lengths_2d(size_t dim[], double bdy[]);
  PetscErrorCode get_vertical_dims(const int z_len, const int zb_len, 
				   double z_read[], double zb_read[]);

  PetscErrorCode put_dimension(int varid, int len, PetscScalar *vals);
  PetscErrorCode put_dimension_regular(int varid, int len, double start, double delta);

  PetscErrorCode read_polar_stereographic(double &straight_vertical_longitude_from_pole,
					  double &latitude_of_projection_origin,
					  double &standard_parallel);
  PetscErrorCode write_polar_stereographic(double straight_vertical_longitude_from_pole,
					   double latitude_of_projection_origin,
					   double standard_parallel);

  PetscErrorCode get_local_var(const int varid, DA da, Vec v, Vec g,
			       const int *s, const int *c,
			       int dims, void *a_mpi, int a_size);
  PetscErrorCode get_global_var(const int varid, DA da, Vec g,
				const int *s, const int *c,
				int dims, void *a_mpi, int a_size);

  PetscErrorCode put_local_var(const int varid, DA da, Vec v, Vec g,
			       const int *s, const int *c,
			       int dims, void *a_mpi, int a_size);
  PetscErrorCode put_global_var(const int varid, DA da, Vec g,
				const int *s, const int *c,
				int dims, void *a_mpi, int a_size);

  PetscErrorCode set_MaskInterp(MaskInterp *mi_in);
  PetscErrorCode regrid_local_var(const int varid, int dim_flag,
				  LocalInterpCtx &lic, DA da, Vec vec, Vec g,
				  bool useMaskInterp);
  PetscErrorCode regrid_global_var(const int varid, int dim_flag,
				   LocalInterpCtx &lic, DA da, Vec g,
				   bool useMaskInterp);

private:
  bool check_dimension(const char dim[], const int len);
  bool check_dimensions();
  MaskInterp  *myMaskInterp;
  IceGrid* grid;
};

#endif // __nc_util_hh


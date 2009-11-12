// Copyright (C) 2009 Constantine Khroulev
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

#ifndef __IceModelVec2T_hh
#define __IceModelVec2T_hh

#include <petsc.h>
#include <string>
#include <vector>
#include "NCVariable.hh"
#include "iceModelVec.hh"
#include "nc_util.hh"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

//! A class for storing and accessing 2D time-series (for climate forcing)
/*! This class was created to read time-dependent and spatially-varying climate
  forcing data, in particular snow temperatures and precipitation.

  It allocates a number (given as an argument to the create() method) of
  records and reads them from a file if necessary.

  If requests (calls to update()) go in sequence, every records should be read
  only once.

  Note that this class is optimized for use with a PDD scheme -- it stores
  records so that data corresponding to a grid point are stored in adjacent
  memory locations.

  IceModelVec2T is always global (%i.e. has no ghosts).

  Both versions of interp() use linear interpolation and extrapolate (by a
  constant) outside the available range.

  Usage example:
  \code
  // initialization:
  char filename[] = "climate_inputs.nc";
  IceModelVec2T v;
  ierr = v.create(grid, "snowtemp", config.get("climate_forcing_buffer_size")); CHKERRQ(ierr);
  ierr = v.set_attrs("climate_forcing", "snow temperature", "K", ""); CHKERRQ(ierr);
  ierr = v.init(filename); CHKERRQ(ierr);

  // actual use:

  // ask it how long a time-step is possible:
  double t = 0, max_dt = 1e16;
  ierr = v.max_timestep(t, max_dt); CHKERRQ(ierr);

  // update (this might read more data from the file)
  ierr = v.update(t, max_dt); CHKERRQ(ierr);

  // get a 2D field (such as precipitation):
  ierr = v.interp(t + max_dt / 2.0); CHKERRQ(ierr);
  // at this point v "looks" almost like an IceModelVec2 with data we need

  // get a time-series for every grid location:
  int N = 21;
  vector<double> ts(N), values(N);
  double dt = max_dt / (N - 1);

  for (int j = 0; j < N; ++j) {
    ts[j] = t + dt * j;
  }

  // is is OK to call update() again, it will not re-read data if at all possible
  ierr = v.update(t, max_dt); CHKERRQ(ierr);

  ierr = v.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      ierr = v.interp(i, j, N, &ts[0], &values[0]); CHKERRQ(ierr);
      // do more
    }
  ierr = v.end_access(); CHKERRQ(ierr);
  
  \endcode
 */
class IceModelVec2T : public IceModelVec2 {
public:
  IceModelVec2T();
  IceModelVec2T(const IceModelVec2T &other);
  virtual PetscErrorCode create(IceGrid &mygrid, const char my_short_name[], int N);
  virtual PetscErrorCode destroy();
  virtual PetscErrorCode get_arrays(PetscScalar** &a2, PetscScalar*** &a3);
  virtual PetscErrorCode init(string filename);
  virtual PetscErrorCode update(double t_years, double dt_years);
  virtual PetscErrorCode set_record(int n);
  virtual PetscErrorCode get_record(int n);
  virtual PetscErrorCode max_timestep(double t_years, double &dt_years);
  virtual PetscErrorCode interp(double t_years);
  virtual PetscErrorCode interp(int i, int j, int N,
				PetscScalar *times, PetscScalar *results);
  virtual PetscErrorCode write(string filename, double t_years, nc_type nctype);
  virtual PetscErrorCode begin_access();
  virtual PetscErrorCode end_access();
private:
  vector<double> times,		//!< all the times available in filename
    T;				//!< times stored in memory
  string filename;		//!< file to read (regrid) from
  LocalInterpCtx *lic;		//!< We store the context because we might need
				//!< to regrid many times
  DA da3;
  Vec v3;			//!< a 3D Vec used to store records
  void ***array3;
  int n_records,		//!< maximum number of records to store in memory
    first;			//!< in-file index of the first record stored in memory
  
  virtual PetscErrorCode update(int start);
  virtual PetscErrorCode discard(int N);
};


#endif // __IceModelVec2T_hh

// Copyright (C) 2009--2012 Constantine Khroulev
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

#include "iceModelVec.hh"

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
  v.set_n_records(config.get("climate_forcing_buffer_size"));
  ierr = v.create(grid, "snowtemp", false); CHKERRQ(ierr);
  ierr = v.set_attrs("climate_forcing", "snow surface temperature", "K", ""); CHKERRQ(ierr);
  ierr = v.init(filename); CHKERRQ(ierr);

  // actual use:

  // ask it how long a time-step is possible:
  double t = 0, max_dt = 1e16;
  ierr = v.max_timestep(t, max_dt); CHKERRQ(ierr);

  // update (this might read more data from the file)
  ierr = v.update(t, max_dt); CHKERRQ(ierr);

  // get a 2D field (such as precipitation):
  ierr = v.interp(t + max_dt / 2.0); CHKERRQ(ierr);
  // at this point v "looks" almost like an IceModelVec2S with data we need (and
  // can be type-cast to IceModelVec2S)

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

  // compute an average value over a time interval at a certain grid location:

  double average;
  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = v.average(grid.xs, grid.ys, t, max_dt, average); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);

  \endcode
 */
class IceModelVec2T : public IceModelVec2S {
public:
  IceModelVec2T();
  IceModelVec2T(const IceModelVec2T &other);
  virtual ~IceModelVec2T();

  virtual void set_n_records(unsigned int N);
  virtual unsigned int get_n_records();
  using IceModelVec2S::create;
  virtual PetscErrorCode create(IceGrid &mygrid, string my_short_name,
                                bool local, int width = 1);
  virtual PetscErrorCode init(string filename);
  virtual PetscErrorCode update(double my_t, double my_dt);
  virtual PetscErrorCode set_record(int n);
  virtual PetscErrorCode get_record(int n);
  virtual double         max_timestep(double my_t);

  virtual PetscErrorCode at_time(double my_t);

  virtual PetscErrorCode interp(double my_t);
  virtual PetscErrorCode interp(int i, int j, int N,
				PetscScalar *times, PetscScalar *results);
  virtual PetscErrorCode average(double my_t, double my_dt);
  virtual PetscErrorCode average(int i, int j, double my_t, double my_dt,
				 double &result);
  virtual PetscErrorCode begin_access();
  virtual PetscErrorCode end_access();

protected:
  vector<double> time,		//!< all the times available in filename
    time_bounds;		//!< time bounds
  string filename;		//!< file to read (regrid) from
  DM da3;
  Vec v3;			//!< a 3D Vec used to store records
  void ***array3;
  int n_records,		//!< maximum number of records to store in memory
    first,			//!< in-file index of the first record stored in memory
    N;                   //!< number of records kept in memory
  LocalInterpCtx *lic;

  virtual PetscErrorCode destroy();
  virtual PetscErrorCode get_array3(PetscScalar*** &a3);
  virtual PetscErrorCode update(int start);
  virtual PetscErrorCode discard(int N);
};


#endif // __IceModelVec2T_hh

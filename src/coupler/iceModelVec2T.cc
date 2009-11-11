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

#include <petsc.h>
#include <algorithm>
#include "iceModelVec2T.hh"
#include "../base/nc_util.hh"
#include "../base/pism_const.hh"

IceModelVec2T::IceModelVec2T() : IceModelVec() {
  localp = false;
  da3 = PETSC_NULL;
  v3 = PETSC_NULL;
  array3 = NULL;
  lic = NULL;
  first = -1;
}

IceModelVec2T::IceModelVec2T(const IceModelVec2T &other) : IceModelVec(other) {
  shallow_copy = true;
  localp = false;

  dimension = other.dimension;
  times = other.times;
  T = other.T;
  filename = other.filename;
  lic = other.lic;
  da3 = other.da3;
  v3 = other.v3;
  array3 = other.array3;
  n_records = other.n_records;
  first = other.first;
}

PetscErrorCode IceModelVec2T::create(IceGrid &my_grid, const char my_short_name[], int my_n_records) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }
  if (v != PETSC_NULL) {
    SETERRQ1(2,"IceModelVec2T with name='%s' already allocated\n", my_short_name);
  }

  grid = &my_grid;
  dims = GRID_2D;
  n_records = my_n_records;
  T.resize(1);			// we assume that T[0] is always available

  // create the 2D DA:
  PetscInt       M, N, m, n;
  PetscErrorCode ierr;
  ierr = DAGetInfo(my_grid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate2d(my_grid.com, DA_XYPERIODIC, DA_STENCIL_BOX,
		    N, M,
		    n, m,
		    1, 1,	// dof and stencil width
                    PETSC_NULL, PETSC_NULL, &da); CHKERRQ(ierr);

  // allocate the 2D Vec:
  ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);

  // create the 3D DA:
  ierr = DACreate3d(my_grid.com, DA_YZPERIODIC, DA_STENCIL_STAR,
		    n_records, N, M,
		    1,         n, m,
		    1, 1,	// dof and stencil width
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da3); CHKERRQ(ierr);

  // allocate the 3D Vec:
  ierr = DACreateGlobalVector(da3, &v3); CHKERRQ(ierr);

  name = my_short_name;

  var1.init(name, my_grid, GRID_2D);

  return 0;
}

PetscErrorCode IceModelVec2T::destroy() {
  PetscErrorCode ierr;

  ierr = IceModelVec::destroy(); CHKERRQ(ierr);

  if (v3 != PETSC_NULL) {
    ierr = VecDestroy(v3); CHKERRQ(ierr);
    v3 = PETSC_NULL;
  }
  if (da3 != PETSC_NULL) {
    ierr = DADestroy(da3); CHKERRQ(ierr);
    da3 = PETSC_NULL;
  }

  delete lic;

  return 0;
}

PetscErrorCode IceModelVec2T::get_arrays(PetscScalar** &a2, PetscScalar*** &a3) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a2 = (PetscScalar**) array;
  a3 = (PetscScalar***) array3;
  return 0;
}

PetscErrorCode IceModelVec2T::begin_access() {
  PetscErrorCode ierr = IceModelVec::begin_access(); CHKERRQ(ierr);
  if (array3 == NULL) {
    ierr = DAVecGetArray(da3, v3, &array3); CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode IceModelVec2T::end_access() {
  PetscErrorCode ierr = IceModelVec::end_access(); CHKERRQ(ierr);
  if (array3 != NULL) {
    ierr = DAVecRestoreArray(da3, v3, &array3); CHKERRQ(ierr);
    array3 = PETSC_NULL;
  }
  return 0;
}

PetscErrorCode IceModelVec2T::init(string fname, string dim_name) {
  PetscErrorCode ierr;
  NCTool nc(grid);
  grid_info gi;

  filename = fname;
  
  dimension.init(dim_name, dim_name, grid->com, grid->rank);
  ierr = dimension.set_units("years"); CHKERRQ(ierr);
  ierr = dimension.read(filename.c_str(), times); CHKERRQ(ierr);
  
  bool is_increasing = true;
  for (unsigned int j = 1; j < times.size(); ++j) {
    if (times[j] - times[j-1] < 1e-16) {
      is_increasing = false;
      break;
    }
  }
  if (!is_increasing) {
    ierr = PetscPrintf(grid->com, "PISM ERROR: times '%s' have to be strictly increasing (read from '%s').\n",
		       dimension.short_name.c_str(), filename.c_str());
    PetscEnd();
  }

  ierr = nc.open_for_reading(filename.c_str()); CHKERRQ(ierr);
  ierr = nc.get_grid_info(gi); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  lic = new LocalInterpCtx(gi, NULL, NULL, *grid);
  lic->report_range = false;

  double dt = 1e16;		// a very big number
  ierr = max_timestep(grid->year, dt); CHKERRQ(ierr);
  ierr = update(grid->year, dt); CHKERRQ(ierr);

  return 0;
}

//! Read some data to make sure that the interval (t_years, t_years + dt_years) is covered.
PetscErrorCode IceModelVec2T::update(double t_years, double dt_years) {
  PetscErrorCode ierr;
  int j, start, N = times.size();
  vector<double>::iterator begin = times.begin();

  // just return if we have all the data we need:
  if ((t_years >= T[0]) && (t_years + dt_years <= T.back()))
    return 0;

  j = upper_bound(begin, times.end(), t_years) - begin - 1;

  if ((n_records >= N) || (j < 0)) {
    start = 0;
  } else {
    start = j;
    
    if ((start < N - n_records) &&
	(t_years + dt_years > times[start + n_records - 1]))
      SETERRQ(1, "IceModelVec2T::update(): timestep is too big");
  }

  ierr = update(start); CHKERRQ(ierr);

  return 0;
}

//! Update by reading at most n_records records from the file.
PetscErrorCode IceModelVec2T::update(int start) {
  PetscErrorCode ierr;
  int N = times.size();

  if ((start < 0) || (start >= N))
    SETERRQ1(1, "IceModelVec2T::update(int start): start = %d is invalid", start);

  int missing = PetscMin(n_records, N - start);
  
  int kept = 0;
  if (first >= 0) {
    int last = first + T.size() - 1;
    if ((start >= first) && (start <= last)) {
      int discarded = start - first; 
      kept = last - start + 1;
      ierr = discard(discarded); CHKERRQ(ierr);
      missing -= kept;
      start += kept;
      first += discarded;
    }
  } else {
    first = start;
  }

  if (missing <= 0) return 0;
  
  T.resize(kept + missing);
  string long_name = string_attr("long_name");
  ierr = verbPrintf(2, grid->com, "  reading \"%s\" (short_name = %s): %d records, years %3.3f through %3.3f...\n",
		    long_name.c_str(), name.c_str(), missing,
		    times[start], times[start + missing - 1]);

  for (int j = 0; j < missing; ++j) {
    lic->start[0] = start + j;
    ierr = regrid(filename.c_str(), *lic, true); CHKERRQ(ierr);
    ierr = verbPrintf(5, grid->com, " %s: reading entry #%02d, year %f...\n",
		      name.c_str(), start + j, times[start + j]);
    ierr = set_record(kept + j); CHKERRQ(ierr);
    T[kept + j] = times[start + j];
  }

  return 0;
}

//! Discard the first N records, shifting the rest of them towards the "beginning".
PetscErrorCode IceModelVec2T::discard(int N) {
  PetscErrorCode ierr;
  PetscScalar **a2, ***a3;
  PetscInt M = T.size() - N;

  for (PetscInt k = 0; k < M; ++k)
    T[k] = T[k + N];
  T.resize(M);

  ierr = get_arrays(a2, a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      for (PetscInt k = 0; k < M; ++k)
	a3[i][j][k] = a3[i][j][k + N];
  ierr = end_access(); CHKERRQ(ierr);
  
  return 0;
}

//! Sets the record number n to the contents of the Vec v.
PetscErrorCode IceModelVec2T::set_record(int n) {
  PetscErrorCode ierr;
  PetscScalar **a2, ***a3;

  ierr = get_arrays(a2, a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      a3[i][j][n] = a2[i][j];
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! Sets the Vec v to the contents of the nth record.
PetscErrorCode IceModelVec2T::get_record(int n) {
  PetscErrorCode ierr;
  PetscScalar **a2, ***a3;

  ierr = get_arrays(a2, a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      a2[i][j] = a3[i][j][n];
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! Given the time t_years and the current selected time-step dt_years,
//! determines the maximum possible time-step this IceModelVec2T allows.
PetscErrorCode IceModelVec2T::max_timestep(double t_years, double &dt_years) {
  int j, k, N = times.size();
  vector<double>::iterator begin = times.begin();

  // no restriction if all the records are in memory:
  if (n_records >= N) {
    return 0;
  }

  // no restriction if t_years is outside the interval of available records (on the right)
  if (t_years >= times[N - 1]) {
    return 0;
  }

  // if t_years is outside the interval on the left, we can go to times[n_records - 2].
  if (t_years < times[0]) {
    dt_years = PetscMin(dt_years, times[n_records - 2] - t_years);
    return 0;
  }

  j = upper_bound(begin, times.end(), t_years) - begin - 1;
  k = j + n_records;

  if (k < N) {
    dt_years = PetscMin(dt_years, times[k - 2] - t_years);
    return 0;
  }

  return 0;
}

//! Writes the snapshot corresponding to t_years to a file.
PetscErrorCode IceModelVec2T::write(string filename, double t_years, nc_type nctype) {
  PetscErrorCode ierr;

  ierr = interp(t_years); CHKERRQ(ierr);
  ierr = IceModelVec::write(filename.c_str(), nctype); CHKERRQ(ierr);

  return 0;
}

//! Extract data corresponding to t_years using linear interpolation.
PetscErrorCode IceModelVec2T::interp(double t_years) {
  PetscErrorCode ierr;
  vector<double>::iterator end = T.end(), j;
  
  j = upper_bound(T.begin(), end, t_years); // binary search

  if (j == end) {
    ierr = get_record(T.size() - 1); CHKERRQ(ierr);
    return 0;
  }
    
  int index = j - T.begin();

  if (index == 0) {
    ierr = get_record(0); CHKERRQ(ierr);
    return 0;
  }

  double lambda = (t_years - T[index]) / (T[index + 1] - T[index]);
  
  PetscScalar **a2, ***a3;

  ierr = get_arrays(a2, a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      a2[i][j] = a3[i][j][index] * (1 - lambda) + a3[i][j][index + 1] * lambda;
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! Gets an interpolated time-series out. Has to be surrounded with begin_access() and end_access().
PetscErrorCode IceModelVec2T::interp(int i, int j, int N,
				    PetscScalar *ts, PetscScalar *values) {
  int mcurr = 0;
  PetscScalar ***a3 = (PetscScalar***) array3;

  for (int k = 0; k < N; ++k) {
    // extrapolate on the left:
    if (ts[k] <= T[0]) {
      values[k] = a3[i][j][0];
      continue;
    }
    // extrapolate on the right:
    if (ts[k] >= T.back()) {
      values[k] = a3[i][j][T.size()-1];
      continue;
    }

    while (T[mcurr+1] < ts[k]) {
      mcurr++;
    }

    const PetscScalar incr = (ts[k] - T[mcurr]) / (T[mcurr+1] - T[mcurr]);
    const PetscScalar valm = a3[i][j][mcurr];
    values[k] = valm + incr * (a3[i][j][mcurr+1] - valm);
  }

  return 0;
}

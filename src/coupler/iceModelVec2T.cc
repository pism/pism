// Copyright (C) 2009--2011 Constantine Khroulev
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
#include "PISMIO.hh"
#include "pism_const.hh"

IceModelVec2T::IceModelVec2T() : IceModelVec2S() {
  localp = false;
  da3 = PETSC_NULL;
  v3 = PETSC_NULL;
  array3 = NULL;
  first = -1;
  n_records = 50;		// just a default
  T.resize(1);			// so that T[0] is always available
  report_range = false;
  lic = NULL;
  strict_timestep_limit = false;
}

IceModelVec2T::IceModelVec2T(const IceModelVec2T &other) : IceModelVec2S(other) {
  shallow_copy = true;

  T = other.T;
  array3 = other.array3;
  da3 = other.da3;
  filename = other.filename;
  first = other.first;
  lic = other.lic;
  localp = other.localp;
  n_records = other.n_records;
  times = other.times;
  v3 = other.v3;
}

IceModelVec2T::~IceModelVec2T() {
  if (!shallow_copy) {
    delete lic;
  }
}


//! Sets the number of records to store in memory. Call it before calling create().
void IceModelVec2T::set_n_records(unsigned int N) {
  n_records = N;
}

PetscErrorCode IceModelVec2T::create(IceGrid &my_grid, const char my_short_name[],
				     bool local, int width) {
  PetscErrorCode ierr;

  if (local) {
    SETERRQ(1, "IceModelVec2T cannot be 'local'");
  }

  ierr = IceModelVec2S::create(my_grid, my_short_name, false, width); CHKERRQ(ierr);

  // create the 3D DA:
  ierr = DACreate3d(my_grid.com, DA_YZPERIODIC, DA_STENCIL_STAR,
		    n_records, grid->My, grid->Mx,
		    1,         grid->Ny, grid->Nx,
		    1, 1,	// dof and stencil width
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da3); CHKERRQ(ierr);

  // allocate the 3D Vec:
  ierr = DACreateGlobalVector(da3, &v3); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2T::destroy() {
  PetscErrorCode ierr;

  ierr = IceModelVec2S::destroy(); CHKERRQ(ierr);

  if (v3 != PETSC_NULL) {
    ierr = VecDestroy(v3); CHKERRQ(ierr);
    v3 = PETSC_NULL;
  }
  if (da3 != PETSC_NULL) {
    ierr = DADestroy(da3); CHKERRQ(ierr);
    da3 = PETSC_NULL;
  }

  return 0;
}

PetscErrorCode IceModelVec2T::get_array3(PetscScalar*** &a3) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a3 = (PetscScalar***) array3;
  return 0;
}

PetscErrorCode IceModelVec2T::begin_access() {
  PetscErrorCode ierr = IceModelVec2S::begin_access(); CHKERRQ(ierr);
  if (array3 == NULL) {
    ierr = DAVecGetArray(da3, v3, &array3); CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode IceModelVec2T::end_access() {
  PetscErrorCode ierr = IceModelVec2S::end_access(); CHKERRQ(ierr);
  if (array3 != NULL) {
    ierr = DAVecRestoreArray(da3, v3, &array3); CHKERRQ(ierr);
    array3 = PETSC_NULL;
  }
  return 0;
}

PetscErrorCode IceModelVec2T::init(string fname) {
  PetscErrorCode ierr;
  NCTimeseries time_dimension;

  filename = fname;

  NCTool nc(grid->com, grid->rank);
  bool t_exists, time_exists;

  ierr = nc.open_for_reading(filename.c_str()); CHKERRQ(ierr);
  ierr = nc.find_variable("t", NULL, t_exists); CHKERRQ(ierr);
  ierr = nc.find_variable("time", NULL, time_exists); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  if (t_exists) {
    time_dimension.init("t", "t", grid->com, grid->rank);
  } else if (time_exists) {
    time_dimension.init("time", "time", grid->com, grid->rank);
  }

  ierr = time_dimension.set_units("years"); CHKERRQ(ierr);
  ierr = time_dimension.read(filename.c_str(), times); CHKERRQ(ierr);
  
  bool is_increasing = true;
  for (unsigned int j = 1; j < times.size(); ++j) {
    if (times[j] - times[j-1] < 1e-16) {
      is_increasing = false;
      break;
    }
  }
  if (!is_increasing) {
    ierr = PetscPrintf(grid->com, "PISM ERROR: times have to be strictly increasing (read from '%s').\n",
		       filename.c_str());
    PISMEnd();
  }

  ierr = get_interp_context(filename, lic); CHKERRQ(ierr);

  return 0;
}

//! Read some data to make sure that the interval (t_years, t_years + dt_years) is covered.
PetscErrorCode IceModelVec2T::update(double t_years, double dt_years) {
  PetscErrorCode ierr;
  int j, start, N = (int)times.size();
  vector<double>::iterator begin = times.begin();

  // just return if we have all the data we need:
  if ((t_years >= T[0]) && (t_years + dt_years <= T.back()))
    return 0;

  j = (int)(upper_bound(begin, times.end(), t_years) - begin - 1);

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
  int N = (int)times.size();

  if ((start < 0) || (start >= N))
    SETERRQ1(1, "IceModelVec2T::update(int start): start = %d is invalid", start);

  int missing = PetscMin(n_records, N - start);
  
  int kept = 0;
  if (first >= 0) {
    int last = first + (int)T.size() - 1;
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
  ierr = verbPrintf(2, grid->com, "  reading \"%s\" records into buffer\n    (short_name = %s): %d records, years %3.3f through %3.3f...\n",
		    long_name.c_str(), name.c_str(), missing,
		    times[start], times[start + missing - 1]);

  for (int j = 0; j < missing; ++j) {
    lic->start[0] = start + j;
    lic->report_range = false;

    ierr = vars[0].regrid(filename.c_str(), *lic, true, false, 0.0, v); CHKERRQ(ierr);

    // ierr = vars[0].read(filename.c_str(), start + j, v); CHKERRQ(ierr);

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
  PetscInt M = (int)T.size() - N;

  for (PetscInt k = 0; k < M; ++k)
    T[k] = T[k + N];
  T.resize(M);

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      for (PetscInt k = 0; k < M; ++k)
	a3[i][j][k] = a3[i][j][k + N];
  ierr = end_access(); CHKERRQ(ierr);
  
  return 0;
}

//! Sets the record number n to the contents of the (internal) Vec v.
PetscErrorCode IceModelVec2T::set_record(int n) {
  PetscErrorCode ierr;
  PetscScalar **a2, ***a3;

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      a3[i][j][n] = a2[i][j];
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! Sets the (internal) Vec v to the contents of the nth record.
PetscErrorCode IceModelVec2T::get_record(int n) {
  PetscErrorCode ierr;
  PetscScalar **a2, ***a3;

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      a2[i][j] = a3[i][j][n];
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Given the time t_years and the current selected time-step dt_years,
//! determines the maximum possible time-step this IceModelVec2T allows.
/*!
  Returns -1 if any time step is OK at t_years.
 */
double IceModelVec2T::max_timestep(double t_years) {
  int j, k, N = (int)times.size();
  vector<double>::iterator begin = times.begin();

  // only allow going to till the next record if strict time-step restriction
  // is requested
  if (strict_timestep_limit) {
    vector<double>::iterator l = upper_bound(begin, times.end(), t_years);
    if (l != times.end()) {
      PetscReal tmp = *l - t_years;

      if (tmp > 1e-8)
        return tmp;
      else if (l + 1 != times.end())
        return *(l + 1) - *l;
      else
        return -1;
    } else
      return -1;
  }

  // no restriction if all the records are in memory:
  if (n_records >= N) {
    return -1;
  }

  // no restriction if t_years is outside the interval of available records (on the right)
  if (t_years >= times.back()) {
    return -1;
  }

  // if t_years is outside the interval on the left, we can go to times[n_records - 2].
  if (t_years < times[0]) {
    return times[n_records - 2] - t_years;
  }

  j = (int)(upper_bound(begin, times.end(), t_years) - begin - 1);
  k = j + n_records;

  if (k < N) {
    return times[k - 2] - t_years;
  }

  return -1;
}

//! \brief Get the record that is just before t_years.
PetscErrorCode IceModelVec2T::get_record_years(double t_years) {
  PetscErrorCode ierr;

  vector<double>::iterator end = T.end(), k;
  
  k = upper_bound(T.begin(), end, t_years); // binary search

  if (k == end) {
    ierr = get_record((int)T.size() - 1); CHKERRQ(ierr);
    return 0;
  }
    
  int index = (int)(k - T.begin() - 1);

  if (index < 0) {
    ierr = get_record(0); CHKERRQ(ierr);
    return 0;
  }

  ierr = get_record(index); CHKERRQ(ierr);

  return 0;
}


//! Extract data corresponding to t_years using linear interpolation.
/*!
  Note: this method does not check if an update() call is necessary!
 */
PetscErrorCode IceModelVec2T::interp(double t_years) {
  PetscErrorCode ierr;
  vector<double>::iterator end = T.end(), k;
  
  k = upper_bound(T.begin(), end, t_years); // binary search

  if (k == end) {
    ierr = get_record((int)T.size() - 1); CHKERRQ(ierr);
    return 0;
  }
    
  int index = (int)(k - T.begin() - 1);

  if (index < 0) {
    ierr = get_record(0); CHKERRQ(ierr);
    return 0;
  }

  double lambda = (t_years - T[index]) / (T[index + 1] - T[index]);
  
  PetscScalar **a2, ***a3;

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      a2[i][j] = a3[i][j][index] * (1 - lambda) + a3[i][j][index + 1] * lambda;
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! Gets an interpolated time-series out. Has to be surrounded with begin_access() and end_access().
/*!
  Note: this method does not check ownership and does not check if an update() call is necessary!
 */
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

//! \brief Finds the average value at i,j over the interval (t_years, t_years +
//! dt_years) using trapezoidal rule.
/*!
  Can (and should) be optimized. Later, though.
 */
PetscErrorCode IceModelVec2T::average(int i, int j, double t_years, double dt_years,
				      double &result) {
  PetscErrorCode ierr;

  // Determine the number of small time-steps to use for averaging:
  int N = (int) ceil(52 * (dt_years) + 1); // (52 weeks in a year)
  if (N < 2) N = 2;

  vector<double> ts(N), values(N);
  double dt = dt_years / (N - 1);
  for (int k = 0; k < N; k++)
    ts[k] = t_years + k * dt;
  
  ierr = interp(i, j, N, &ts[0], &values[0]); CHKERRQ(ierr);

  // trapezoidal rule; uses the fact that all 'small' time intervals used here
  // are the same:
  result = 0;
  for (int k = 0; k < N - 1; ++k)
    result += values[k] + values[k + 1];
  result /= 2*(N - 1);

  return 0;
}

PetscErrorCode IceModelVec2T::average(double t_years, double dt_years) {
  PetscErrorCode ierr;
  PetscScalar **a2;

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = get_array(a2);
  for (PetscInt   i = grid->xs; i < grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j < grid->ys+grid->ym; ++j) {
      ierr = average(i, j, t_years, dt_years, a2[i][j]); CHKERRQ(ierr);
    }
  }

  ierr = end_access(); CHKERRQ(ierr);
  return 0;
}


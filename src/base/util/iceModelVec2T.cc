// Copyright (C) 2009--2014 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include "PIO.hh"
#include "pism_const.hh"
#include "PISMTime.hh"
#include "LocalInterpCtx.hh"
#include "IceGrid.hh"

namespace pism {

IceModelVec2T::IceModelVec2T() : IceModelVec2S() {
  m_has_ghosts           = false;
  m_v3                     = NULL;
  array3                 = NULL;
  first                  = -1;
  N                      = 0;
  n_records              = 50;  // just a default
  m_report_range         = false;
  m_period               = 0;
  m_reference_time       = 0.0;
  n_evaluations_per_year = 53;

  m_da3.reset();
}

IceModelVec2T::~IceModelVec2T() {
  destroy();
}


//! Sets the number of records to store in memory. Call it before calling create().
void IceModelVec2T::set_n_records(unsigned int my_N) {
  n_records = my_N;
}

void IceModelVec2T::set_n_evaluations_per_year(unsigned int M) {
  n_evaluations_per_year = M;
}

unsigned int IceModelVec2T::get_n_records() {
  return n_records;
}

PetscErrorCode IceModelVec2T::create(IceGrid &my_grid, const std::string &my_short_name,
                                     bool local, int width) {
  PetscErrorCode ierr;

  if (local) {
    SETERRQ(grid->com, 1, "IceModelVec2T cannot be 'local'");
  }

  ierr = IceModelVec2S::create(my_grid, my_short_name, WITHOUT_GHOSTS, width); CHKERRQ(ierr);

  // initialize the m_da3 member:
  ierr = grid->get_dm(this->n_records, this->m_da_stencil_width, m_da3); CHKERRQ(ierr);

  // allocate the 3D Vec:
  ierr = DMCreateGlobalVector(m_da3->get(), &m_v3); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2T::destroy() {
  PetscErrorCode ierr;

  ierr = IceModelVec2S::destroy(); CHKERRQ(ierr);

  if (m_v3 != NULL) {
    ierr = VecDestroy(&m_v3); CHKERRQ(ierr);
    m_v3 = NULL;
  }

  return 0;
}

PetscErrorCode IceModelVec2T::get_array3(double*** &a3) {
  PetscErrorCode ierr = begin_access(); CHKERRQ(ierr);
  a3 = (double***) array3;
  return 0;
}

PetscErrorCode IceModelVec2T::begin_access() const {
  PetscErrorCode ierr;
  if (m_access_counter == 0) {
    ierr = DMDAVecGetArrayDOF(m_da3->get(), m_v3, &array3); CHKERRQ(ierr);
  }

  // this call will increment the m_access_counter
  ierr = IceModelVec2S::begin_access(); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode IceModelVec2T::end_access() const {
  // this call will decrement the m_access_counter
  PetscErrorCode ierr = IceModelVec2S::end_access(); CHKERRQ(ierr);

  if (m_access_counter == 0) {
    ierr = DMDAVecRestoreArrayDOF(m_da3->get(), m_v3, &array3); CHKERRQ(ierr);
    array3 = NULL;
  }

  return 0;
}

PetscErrorCode IceModelVec2T::init(const std::string &fname, unsigned int period, double reference_time) {
  PetscErrorCode ierr;

  filename         = fname;
  m_period         = period;
  m_reference_time = reference_time;

  // We find the variable in the input file and
  // try to find the corresponding time dimension.

  PIO nc(*grid, "guess_mode");
  std::string name_found;
  bool exists, found_by_standard_name;
  ierr = nc.open(filename, PISM_READONLY); CHKERRQ(ierr);
  ierr = nc.inq_var(m_metadata[0].get_name(), m_metadata[0].get_string("standard_name"),
                    exists, name_found, found_by_standard_name); CHKERRQ(ierr);
  if (exists == false) {
    PetscPrintf(grid->com, "PISM ERROR: can't find %s (%s) in %s.\n",
                m_metadata[0].get_string("long_name").c_str(), m_metadata[0].get_name().c_str(),
                filename.c_str());
    PISMEnd();
  }

  // find the time dimension:
  std::vector<std::string> dims;
  ierr = nc.inq_vardims(name_found, dims); CHKERRQ(ierr);
  
  std::string dimname = "";
  bool time_found = false;
  for (unsigned int i = 0; i < dims.size(); ++i) {
    AxisType dimtype;
    dimname = dims[i];

    ierr = nc.inq_dimtype(dimname, dimtype); CHKERRQ(ierr);

    if (dimtype == T_AXIS) {
      time_found = true;
      break;
    }
  }

  if (time_found) {
    // we're found the time dimension
    NCTimeseries time_dimension(dimname, dimname, grid->get_unit_system());

    ierr = time_dimension.set_units(grid->time->units_string()); CHKERRQ(ierr);
    ierr = nc.read_timeseries(time_dimension, grid->time, time); CHKERRQ(ierr);

    std::string bounds_name;
    ierr = nc.get_att_text(dimname, "bounds", bounds_name); CHKERRQ(ierr);

    if (time.size() > 1) {
      if (bounds_name.empty() == false) {
        // read time bounds data from a file
        NCTimeBounds tb(bounds_name, dimname, grid->get_unit_system());
        ierr = tb.set_units(time_dimension.get_string("units")); CHKERRQ(ierr);

        ierr = nc.read_time_bounds(tb, grid->time, time_bounds); CHKERRQ(ierr);

        // time bounds data overrides the time variable: we make t[j] be the
        // right end-point of the j-th interval
        for (unsigned int k = 0; k < time.size(); ++k)
          time[k] = time_bounds[2*k + 1];
      } else {
        // no time bounds attribute
        PetscPrintf(grid->com,
                    "PISM ERROR: Variable '%s' does not have the time_bounds attribute.\n"
                    "  Cannot use time-dependent forcing data '%s' (%s) without time bounds.\n",
                    dimname.c_str(),  m_metadata[0].get_string("long_name").c_str(), m_metadata[0].get_name().c_str());
        PISMEnd();
      }
    } else {
      // only one time record; set fake time bounds:
      time_bounds.resize(2);
      time_bounds[0] = time[0] - 1;
      time_bounds[1] = time[0] + 1;
    }

  } else {
    // no time dimension; assume that we have only one record and set the time
    // to 0
    time.resize(1);
    time[0] = 0;

    // set fake time bounds:
    time_bounds.resize(2);
    time_bounds[0] = -1;
    time_bounds[1] =  1;
  }

  if (is_increasing(time) == false) {
    ierr = PetscPrintf(grid->com, "PISM ERROR: times have to be strictly increasing (read from '%s').\n",
                       filename.c_str());
    PISMEnd();
  }

  ierr = nc.close(); CHKERRQ(ierr);

  if (m_period != 0) {
    if ((size_t)n_records < time.size())
      SETERRQ(grid->com, 1, "buffer has to be big enough to hold all records of periodic data");

    // read periodic data right away (we need to hold it all in memory anyway)
    ierr = update(0); CHKERRQ(ierr);
  }

  return 0;
}

//! Initialize as constant in time and space
PetscErrorCode IceModelVec2T::init_constant(double value) {
  PetscErrorCode ierr;

  // set constant value everywhere
  ierr = set(value); CHKERRQ(ierr);

  // set the time to zero
  time.resize(1);
  time[0] = 0;
  //N = 1 ;

  // set fake time bounds:
  time_bounds.resize(2);
  time_bounds[0] = -1;
  time_bounds[1] =  1;

  return 0;
}

//! Read some data to make sure that the interval (my_t, my_t + my_dt) is covered.
PetscErrorCode IceModelVec2T::update(double my_t, double my_dt) {
  PetscErrorCode ierr;
  std::vector<double>::iterator i, j;
  unsigned int m, n, last;

  if (time_bounds.size() == 0) {
    ierr = update(0); CHKERRQ(ierr);
    return 0;
  }

  if (m_period != 0) {
    // we read all data in IceModelVec2T::init() (see above)
    return 0;
  }

  if (N > 0) {
    last = first + (N - 1);

    // find the interval covered by data held in memory:
    double t0 = time_bounds[first * 2],
      t1 = time_bounds[last * 2 + 1];

    // just return if we have all the data we need:
    if (my_t >= t0 && my_t + my_dt <= t1)
      return 0;
  }

  // i will point to the smallest t so that t >= my_t 
  i = lower_bound(time_bounds.begin(), time_bounds.end(), my_t);

  // j will point to the smallest t so that t >= my_t + my_dt
  j = lower_bound(time_bounds.begin(), time_bounds.end(), my_t + my_dt);

  // some of the the interval (my_t, my_t + my_dt) is outside of the
  // time interval available in the file (on the right)
  // 
  if (i == time_bounds.end()) {
    m = (int)(time.size() - 1);
  } else {
    m = (int)((i - time_bounds.begin() - 1) / 2);
  }

  if (j == time_bounds.end()) {
    n = (int)(time.size() - 1);
  } else {
    n = (int)((j - time_bounds.begin() - 1) / 2);
  }

  // check if all the records necessary to cover this interval fit in the
  // buffer:
  if (n - m + 1 > n_records)
    SETERRQ(grid->com, 1, "IceModelVec2T::update(): timestep is too big");

  ierr = update(m); CHKERRQ(ierr);

  return 0;
}

//! Update by reading at most n_records records from the file.
PetscErrorCode IceModelVec2T::update(unsigned int start) {
  PetscErrorCode ierr;
  unsigned int time_size = (int)time.size();

  if (start >= time_size)
    SETERRQ1(grid->com, 1, "IceModelVec2T::update(int start): start = %d is invalid", start);

  unsigned int missing = PetscMin(n_records, time_size - start);

  if (start == static_cast<unsigned int>(first))
    return 0;                   // nothing to do

  int kept = 0;
  if (first >= 0) {
    unsigned int last = first + (N - 1);
    if ((N > 0) && (start >= (unsigned int)first) && (start <= last)) {
      int discarded = start - first;
      kept = last - start + 1;
      ierr = discard(discarded); CHKERRQ(ierr);
      missing -= kept;
      start += kept;
      first += discarded;
    } else {
      first = start;
    }
  } else {
    first = start;
  }

  if (missing <= 0) return 0;
  
  N = kept + missing;

  if (this->get_n_records() > 1 || getVerbosityLevel() > 4) {
    ierr = verbPrintf(2, grid->com,
                      "  reading \"%s\" into buffer\n"
                      "          (short_name = %s): %d records, time intervals (%s, %s) through (%s, %s)...\n",
                      metadata().get_string("long_name").c_str(), m_name.c_str(), missing,
                      grid->time->date(time_bounds[start*2]).c_str(),
                      grid->time->date(time_bounds[start*2 + 1]).c_str(),
                      grid->time->date(time_bounds[(start + missing - 1)*2]).c_str(),
                      grid->time->date(time_bounds[(start + missing - 1)*2 + 1]).c_str()); CHKERRQ(ierr);
    m_report_range = false;
  } else {
    m_report_range = true;
  }

  PIO nc(*grid, "guess_mode");
  ierr = nc.open(filename, PISM_READONLY); CHKERRQ(ierr);

  for (unsigned int j = 0; j < missing; ++j) {
    ierr = m_metadata[0].regrid(nc, start + j,
                                CRITICAL, m_report_range, 0.0, m_v); CHKERRQ(ierr);

    ierr = verbPrintf(5, grid->com, " %s: reading entry #%02d, year %s...\n",
                      m_name.c_str(),
                      start + j,
                      grid->time->date(time[start + j]).c_str());
    ierr = set_record(kept + j); CHKERRQ(ierr);
  }

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Discard the first N records, shifting the rest of them towards the "beginning".
PetscErrorCode IceModelVec2T::discard(int number) {
  PetscErrorCode ierr;
  double **a2, ***a3;

  if (number == 0)
    return 0;

  N -= number;

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (unsigned int k = 0; k < N; ++k) {
      a3[i][j][k] = a3[i][j][k + number];
    }
  }
  ierr = end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);
  
  return 0;
}

//! Sets the record number n to the contents of the (internal) Vec v.
PetscErrorCode IceModelVec2T::set_record(int n) {
  PetscErrorCode ierr;
  double **a2, ***a3;

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    a3[i][j][n] = a2[i][j];
  }
  ierr = end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! Sets the (internal) Vec v to the contents of the nth record.
PetscErrorCode IceModelVec2T::get_record(int n) {
  PetscErrorCode ierr;
  double **a2, ***a3;

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    a2[i][j] = a3[i][j][n];
  }
  ierr = end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Given the time my_t and the current selected time-step my_dt,
//! determines the maximum possible time-step this IceModelVec2T allows.
/*!
  Returns -1 if any time step is OK at my_t.
 */
double IceModelVec2T::max_timestep(double my_t) {
  // only allow going to the next record
  std::vector<double>::iterator l = upper_bound(time_bounds.begin(),
                                           time_bounds.end(), my_t);
  if (l != time_bounds.end()) {
    double tmp = *l - my_t;

    if (tmp > 1)                // never take time-steps shorter than 1 second
      return tmp;
    else if ((l + 1) != time_bounds.end() && (l + 2) != time_bounds.end())
      return *(l + 2) - *l;
    else
      return -1;
  } else
    return -1;

}

/*
 * \brief Use piecewise-constant interpolation to initialize
 * IceModelVec2T with the value at time `my_t`.
 *
 * \note This method does not check if an update() call is necessary!
 *
 * @param[in] my_t requested time
 *
 * @return 0 on success
 */
PetscErrorCode IceModelVec2T::interp(double my_t) {
  PetscErrorCode ierr;

  std::vector<double> t_vector(1);
  t_vector[0] = my_t;
  ierr = init_interpolation(t_vector); CHKERRQ(ierr);

  ierr = get_record(m_interp_indices[0]); CHKERRQ(ierr);

  return 0;
}


/** 
 * Compute the average value over the time interval `[my_t, my_t + my_dt]`.
 *
 * @param my_t  start of the time interval, in seconds
 * @param my_dt length of the time interval, in seconds
 *
 * @return 0 on success
 */
PetscErrorCode IceModelVec2T::average(double my_t, double my_dt) {
  PetscErrorCode ierr;
  double **a2;
  double dt_years = grid->convert(my_dt, "seconds", "years"); // *not* time->year(my_dt)

  // if only one record, nothing to do
  if (time.size() == 1) {
    return 0;
  }

  // Determine the number of small time-steps to use for averaging:
  int M = (int) ceil(n_evaluations_per_year * (dt_years));
  if (M < 1)
    M = 1;

  std::vector<double> ts(M);
  double dt = my_dt / M;
  for (int k = 0; k < M; k++)
    ts[k] = my_t + k * dt;

  ierr = init_interpolation(ts); CHKERRQ(ierr);

  ierr = get_array(a2);         // calls begin_access()
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    ierr = average(i, j, a2[i][j]); CHKERRQ(ierr); // NB! order
  }
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

/**
 * \brief Compute weights for the piecewise-constant interpolation.
 * This is used *both* for time-series and "snapshots".
 *
 * @param ts requested times, in seconds
 *
 * @return 0 on success
 */
PetscErrorCode IceModelVec2T::init_interpolation(const std::vector<double> &ts) {
  unsigned int index = 0,
    last = first + N - 1;

  size_t ts_length = ts.size();

  // Compute "periodized" times if necessary.
  std::vector<double> times_requested(ts_length);
  if (m_period != 0) {
    for (unsigned int k = 0; k < ts_length; ++k)
      times_requested[k] = grid->time->mod(ts[k] - m_reference_time, m_period);
  } else {
    for (unsigned int k = 0; k < ts_length; ++k)
      times_requested[k] = ts[k];
  }

  m_interp_indices.resize(ts_length);

  if (time_bounds.size() == 0) {
    for (unsigned int k = 0; k < ts_length; ++k) {
      m_interp_indices[k] = 0;
    }
  }

  for (unsigned int k = 0; k < ts_length; ++k) {

    if (k > 0 && times_requested[k] < times_requested[k-1])
      index = 0; // reset the index: times_requested are not increasing!

    // extrapolate on the left:
    if (times_requested[k] < time_bounds[2*first]) {
      m_interp_indices[k] = 0;
      continue;
    }

    // extrapolate on the right:
    if (times_requested[k] >= time_bounds[2*last + 1]) {
      m_interp_indices[k] = N - 1;
      continue;
    }

    while (index < N) {
      if (time_bounds[2*(first + index) + 0] <= times_requested[k] &&
          time_bounds[2*(first + index) + 1] >  times_requested[k])
        break;

      index++;
    }

    m_interp_indices[k] = index;
  }

  return 0;
}

/** 
 * \brief Compute values of the time-series using precomputed indices
 * (and piecewise-constant interpolation).
 *
 * @param i,j map-plane grid point
 * @param result pointer to an allocated array of `weights.size()` `double`
 *
 * @return 0 on success
 */
PetscErrorCode IceModelVec2T::interp(int i, int j, std::vector<double> &result) {
  double ***a3 = (double***) array3;
  unsigned int ts_length = m_interp_indices.size();

  for (unsigned int k = 0; k < ts_length; ++k) {
    result[k] = a3[i][j][m_interp_indices[k]];
  }
  
  return 0;
}

//! \brief Finds the average value at i,j over the interval (my_t, my_t +
//! my_dt) using the rectangle rule.
/*!
  Can (and should) be optimized. Later, though.
 */
PetscErrorCode IceModelVec2T::average(int i, int j, double &result) {
  PetscErrorCode ierr;
  unsigned int M = m_interp_indices.size();

  if (N == 1) {
    double ***a3 = (double***) array3;
    result = a3[i][j][0];
  } else {
    std::vector<double> values(M);

    ierr = interp(i, j, values); CHKERRQ(ierr);

    // rectangular rule (uses the fact that points are equally-spaces
    // in time)
    result = 0;
    for (unsigned int k = 0; k < M; ++k)
      result += values[k];
    result /= (double)M;
  }
  return 0;
}



} // end of namespace pism

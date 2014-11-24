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

#include "error_handling.hh"

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

void IceModelVec2T::create(IceGrid &my_grid, const std::string &my_short_name,
                                     bool local, int width) {

  if (local) {
    throw RuntimeError("IceModelVec2T cannot be 'local'");
  }

  IceModelVec2S::create(my_grid, my_short_name, WITHOUT_GHOSTS, width);

  // initialize the m_da3 member:
  m_da3 = grid->get_dm(this->n_records, this->m_da_stencil_width);

  // allocate the 3D Vec:
  DMCreateGlobalVector(*m_da3, &m_v3);
}

void IceModelVec2T::destroy() {
  PetscErrorCode ierr;

  IceModelVec2S::destroy();

  if (m_v3 != NULL) {
    ierr = VecDestroy(&m_v3);
    PISM_PETSC_CHK(ierr, "VecDestroy");
    m_v3 = NULL;
  }
}

void IceModelVec2T::get_array3(double*** &a3) {
  begin_access();
  a3 = (double***) array3;
}

void IceModelVec2T::begin_access() const {
  if (m_access_counter == 0) {
    DMDAVecGetArrayDOF(*m_da3, m_v3, &array3);
  }

  // this call will increment the m_access_counter
  IceModelVec2S::begin_access();
  
}

void IceModelVec2T::end_access() const {
  // this call will decrement the m_access_counter
  IceModelVec2S::end_access();

  if (m_access_counter == 0) {
    DMDAVecRestoreArrayDOF(*m_da3, m_v3, &array3);
    array3 = NULL;
  }
}

void IceModelVec2T::init(const std::string &fname, unsigned int period, double reference_time) {

  filename         = fname;
  m_period         = period;
  m_reference_time = reference_time;

  // We find the variable in the input file and
  // try to find the corresponding time dimension.

  PIO nc(*grid, "guess_mode");
  std::string name_found;
  bool exists, found_by_standard_name;
  nc.open(filename, PISM_READONLY);
  nc.inq_var(m_metadata[0].get_name(), m_metadata[0].get_string("standard_name"),
             exists, name_found, found_by_standard_name);
  if (exists == false) {
    throw RuntimeError::formatted("can't find %s (%s) in %s.",
                                  m_metadata[0].get_string("long_name").c_str(), m_metadata[0].get_name().c_str(),
                                  filename.c_str());
  }

  // find the time dimension:
  std::vector<std::string> dims;
  dims = nc.inq_vardims(name_found);
  
  std::string dimname = "";
  bool time_found = false;
  for (unsigned int i = 0; i < dims.size(); ++i) {
    dimname = dims[i];

    AxisType dimtype = nc.inq_dimtype(dimname);

    if (dimtype == T_AXIS) {
      time_found = true;
      break;
    }
  }

  if (time_found) {
    // we're found the time dimension
    NCTimeseries time_dimension(dimname, dimname, grid->get_unit_system());

    time_dimension.set_units(grid->time->units_string());
    nc.read_timeseries(time_dimension, grid->time, time);

    std::string bounds_name = nc.get_att_text(dimname, "bounds");

    if (time.size() > 1) {
      if (bounds_name.empty() == false) {
        // read time bounds data from a file
        NCTimeBounds tb(bounds_name, dimname, grid->get_unit_system());
        tb.set_units(time_dimension.get_string("units"));

        nc.read_time_bounds(tb, grid->time, time_bounds);

        // time bounds data overrides the time variable: we make t[j] be the
        // right end-point of the j-th interval
        for (unsigned int k = 0; k < time.size(); ++k) {
          time[k] = time_bounds[2*k + 1];
        }
      } else {
        // no time bounds attribute
        throw RuntimeError::formatted("Variable '%s' does not have the time_bounds attribute.\n"
                                      "Cannot use time-dependent forcing data '%s' (%s) without time bounds.",
                                      dimname.c_str(),  m_metadata[0].get_string("long_name").c_str(),
                                      m_metadata[0].get_name().c_str());
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
    throw RuntimeError::formatted("times have to be strictly increasing (read from '%s').",
                                  filename.c_str());
  }

  nc.close();

  if (m_period != 0) {
    if ((size_t)n_records < time.size()) {
      throw RuntimeError("buffer has to be big enough to hold all records of periodic data");
    }

    // read periodic data right away (we need to hold it all in memory anyway)
    update(0);
  }
}

//! Initialize as constant in time and space
void IceModelVec2T::init_constant(double value) {

  // set constant value everywhere
  set(value);

  // set the time to zero
  time.resize(1);
  time[0] = 0;
  //N = 1 ;

  // set fake time bounds:
  time_bounds.resize(2);
  time_bounds[0] = -1;
  time_bounds[1] =  1;
}

//! Read some data to make sure that the interval (my_t, my_t + my_dt) is covered.
void IceModelVec2T::update(double my_t, double my_dt) {
  std::vector<double>::iterator i, j;
  unsigned int m, n, last;

  if (time_bounds.size() == 0) {
    update(0);
    return;
  }

  if (m_period != 0) {
    // we read all data in IceModelVec2T::init() (see above)
    return;
  }

  if (N > 0) {
    last = first + (N - 1);

    // find the interval covered by data held in memory:
    double t0 = time_bounds[first * 2],
      t1 = time_bounds[last * 2 + 1];

    // just return if we have all the data we need:
    if (my_t >= t0 && my_t + my_dt <= t1) {
      return;
    }
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
  if (n - m + 1 > n_records) {
    throw RuntimeError("IceModelVec2T::update(): timestep is too big");
  }

  update(m);
}

//! Update by reading at most n_records records from the file.
void IceModelVec2T::update(unsigned int start) {

  unsigned int time_size = (int)time.size();

  if (start >= time_size) {
    throw RuntimeError::formatted("IceModelVec2T::update(int start): start = %d is invalid", start);
  }

  unsigned int missing = PetscMin(n_records, time_size - start);

  if (start == static_cast<unsigned int>(first)) {
    // nothing to do
    return;
  }

  int kept = 0;
  if (first >= 0) {
    unsigned int last = first + (N - 1);
    if ((N > 0) && (start >= (unsigned int)first) && (start <= last)) {
      int discarded = start - first;
      kept = last - start + 1;
      discard(discarded);
      missing -= kept;
      start += kept;
      first += discarded;
    } else {
      first = start;
    }
  } else {
    first = start;
  }

  if (missing <= 0) {
    return;
  }
  
  N = kept + missing;

  if (this->get_n_records() > 1 || getVerbosityLevel() > 4) {
    verbPrintf(2, grid->com,
               "  reading \"%s\" into buffer\n"
               "          (short_name = %s): %d records, time intervals (%s, %s) through (%s, %s)...\n",
               metadata().get_string("long_name").c_str(), m_name.c_str(), missing,
               grid->time->date(time_bounds[start*2]).c_str(),
               grid->time->date(time_bounds[start*2 + 1]).c_str(),
               grid->time->date(time_bounds[(start + missing - 1)*2]).c_str(),
               grid->time->date(time_bounds[(start + missing - 1)*2 + 1]).c_str());
    m_report_range = false;
  } else {
    m_report_range = true;
  }

  PIO nc(*grid, "guess_mode");
  nc.open(filename, PISM_READONLY);

  for (unsigned int j = 0; j < missing; ++j) {
    m_metadata[0].regrid(nc, start + j,
                         CRITICAL, m_report_range, 0.0, m_v);

    verbPrintf(5, grid->com, " %s: reading entry #%02d, year %s...\n",
                      m_name.c_str(),
                      start + j,
                      grid->time->date(time[start + j]).c_str());
    set_record(kept + j);
  }

  nc.close();
}

//! Discard the first N records, shifting the rest of them towards the "beginning".
void IceModelVec2T::discard(int number) {
  double **a2, ***a3;

  if (number == 0) {
    return;
  }

  N -= number;

  get_array(a2);
  get_array3(a3);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (unsigned int k = 0; k < N; ++k) {
      a3[i][j][k] = a3[i][j][k + number];
    }
  }
  end_access();
  end_access();
  
}

//! Sets the record number n to the contents of the (internal) Vec v.
void IceModelVec2T::set_record(int n) {
  double **a2, ***a3;

  get_array(a2);
  get_array3(a3);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    a3[i][j][n] = a2[i][j];
  }
  end_access();
  end_access();
}

//! Sets the (internal) Vec v to the contents of the nth record.
void IceModelVec2T::get_record(int n) {
  double **a2, ***a3;

  get_array(a2);
  get_array3(a3);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    a2[i][j] = a3[i][j][n];
  }
  end_access();
  end_access();
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

    if (tmp > 1) {                // never take time-steps shorter than 1 second
      return tmp;
    } else if ((l + 1) != time_bounds.end() && (l + 2) != time_bounds.end()) {
      return *(l + 2) - *l;
    } else {
      return -1;
    }
  } else {
    return -1;
  }

}

/*
 * \brief Use piecewise-constant interpolation to initialize
 * IceModelVec2T with the value at time `my_t`.
 *
 * \note This method does not check if an update() call is necessary!
 *
 * @param[in] my_t requested time
 *
 */
void IceModelVec2T::interp(double my_t) {

  std::vector<double> t_vector(1);
  t_vector[0] = my_t;
  init_interpolation(t_vector);

  get_record(m_interp_indices[0]);
}


/** 
 * Compute the average value over the time interval `[my_t, my_t + my_dt]`.
 *
 * @param my_t  start of the time interval, in seconds
 * @param my_dt length of the time interval, in seconds
 *
 */
void IceModelVec2T::average(double my_t, double my_dt) {
  double **a2;
  double dt_years = grid->convert(my_dt, "seconds", "years"); // *not* time->year(my_dt)

  // if only one record, nothing to do
  if (time.size() == 1) {
    return;
  }

  // Determine the number of small time-steps to use for averaging:
  int M = (int) ceil(n_evaluations_per_year * (dt_years));
  if (M < 1) {
    M = 1;
  }

  std::vector<double> ts(M);
  double dt = my_dt / M;
  for (int k = 0; k < M; k++) {
    ts[k] = my_t + k * dt;
  }

  init_interpolation(ts);

  get_array(a2);         // calls begin_access()
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    average(i, j, a2[i][j]); // NB! order
  }
  end_access();
}

/**
 * \brief Compute weights for the piecewise-constant interpolation.
 * This is used *both* for time-series and "snapshots".
 *
 * @param ts requested times, in seconds
 *
 */
void IceModelVec2T::init_interpolation(const std::vector<double> &ts) {
  unsigned int index = 0,
    last = first + N - 1;

  size_t ts_length = ts.size();

  // Compute "periodized" times if necessary.
  std::vector<double> times_requested(ts_length);
  if (m_period != 0) {
    for (unsigned int k = 0; k < ts_length; ++k) {
      times_requested[k] = grid->time->mod(ts[k] - m_reference_time, m_period);
    }
  } else {
    for (unsigned int k = 0; k < ts_length; ++k) {
      times_requested[k] = ts[k];
    }
  }

  m_interp_indices.resize(ts_length);

  if (time_bounds.size() == 0) {
    for (unsigned int k = 0; k < ts_length; ++k) {
      m_interp_indices[k] = 0;
    }
  }

  for (unsigned int k = 0; k < ts_length; ++k) {

    if (k > 0 && times_requested[k] < times_requested[k-1]) {
      // reset the index: times_requested are not increasing!
      index = 0;
    }

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
          time_bounds[2*(first + index) + 1] >  times_requested[k]) {
        break;
      }

      index++;
    }

    m_interp_indices[k] = index;
  }
}

/** 
 * \brief Compute values of the time-series using precomputed indices
 * (and piecewise-constant interpolation).
 *
 * @param i,j map-plane grid point
 * @param result pointer to an allocated array of `weights.size()` `double`
 *
 */
void IceModelVec2T::interp(int i, int j, std::vector<double> &result) {
  double ***a3 = (double***) array3;
  unsigned int ts_length = m_interp_indices.size();

  for (unsigned int k = 0; k < ts_length; ++k) {
    result[k] = a3[i][j][m_interp_indices[k]];
  }
}

//! \brief Finds the average value at i,j over the interval (my_t, my_t +
//! my_dt) using the rectangle rule.
/*!
  Can (and should) be optimized. Later, though.
 */
void IceModelVec2T::average(int i, int j, double &result) {
  unsigned int M = m_interp_indices.size();

  if (N == 1) {
    double ***a3 = (double***) array3;
    result = a3[i][j][0];
  } else {
    std::vector<double> values(M);

    interp(i, j, values);

    // rectangular rule (uses the fact that points are equally-spaces
    // in time)
    result = 0;
    for (unsigned int k = 0; k < M; ++k) {
      result += values[k];
    }
    result /= (double)M;
  }
}



} // end of namespace pism

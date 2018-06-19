// Copyright (C) 2009--2018 Constantine Khroulev
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
#include <cassert>

#include "iceModelVec2T.hh"
#include "pism/util/io/PIO.hh"
#include "pism_utilities.hh"
#include "Time.hh"
#include "IceGrid.hh"
#include "ConfigInterface.hh"

#include "error_handling.hh"
#include "io/io_helpers.hh"
#include "pism/util/Logger.hh"

namespace pism {


/*!
 * Allocate an instance that will be used to load and use a forcing field from a file.
 *
 * Checks the number of records in a file and allocates storage accordingly.
 *
 * If `periodic` is true, allocate enough storage to hold all the records, otherwise
 * allocate storage for at most `max_buffer_size` records.
 *
 * @param[in] grid computational grid
 * @param[in] file input file
 * @param[in] short_name variable name in `file`
 * @param[in] standard_name standard name (if available); leave blank to ignore
 * @param[in] max_buffer_size maximum buffer size for non-periodic fields
 * @param[in] evaluations_per_year number of evaluations per year to use when averaging
 * @param[in] periodic true if this forcing field should be interpreted as periodic
 */
IceModelVec2T::Ptr IceModelVec2T::ForcingField(IceGrid::ConstPtr grid,
                                               const PIO &file,
                                               const std::string &short_name,
                                               const std::string &standard_name,
                                               int max_buffer_size,
                                               int evaluations_per_year,
                                               bool periodic) {

  int n_records = file.inq_nrecords(short_name, standard_name,
                                    grid->ctx()->unit_system());

  if (not periodic) {
    n_records = std::min(n_records, max_buffer_size);
  }
  // In the periodic case we try to keep all the records in RAM.

  // Allocate storage for one record if the variable was not found. This is needed to be
  // able to cheaply allocate and then discard an "-atmosphere given" model
  // (atmosphere::Given) when "-surface given" (Given) is selected.
  n_records = std::max(n_records, 1);

  return IceModelVec2T::Ptr(new IceModelVec2T(grid, short_name, n_records,
                                              evaluations_per_year));
}


IceModelVec2T::IceModelVec2T(IceGrid::ConstPtr grid, const std::string &short_name,
                             unsigned int n_records,
                             unsigned int n_evaluations_per_year)
  : IceModelVec2S() {
  m_has_ghosts           = false;
  m_array3                 = NULL;
  m_first                  = -1;
  m_N                      = 0;
  m_report_range         = false;
  m_period               = 0;
  m_reference_time       = 0.0;
  m_n_evaluations_per_year = 53;

  m_n_evaluations_per_year = n_evaluations_per_year;

  const unsigned int width = 1;

  IceModelVec2S::create(grid, short_name, WITHOUT_GHOSTS, width);

  // initialize the m_da3 member:
  m_da3 = m_grid->get_dm(n_records, this->m_da_stencil_width);
  m_n_records = n_records;

  // allocate the 3D Vec:
  PetscErrorCode ierr = DMCreateGlobalVector(*m_da3, m_v3.rawptr());
  PISM_CHK(ierr, "DMCreateGlobalVector");
}

IceModelVec2T::~IceModelVec2T() {
  // empty
}

unsigned int IceModelVec2T::n_records() {
  return m_n_records;
}

double*** IceModelVec2T::get_array3() {
  begin_access();
  return reinterpret_cast<double***>(m_array3);
}

void IceModelVec2T::begin_access() const {
  if (m_access_counter == 0) {
    PetscErrorCode ierr = DMDAVecGetArrayDOF(*m_da3, m_v3, &m_array3);
    PISM_CHK(ierr, "DMDAVecGetArrayDOF");
  }

  // this call will increment the m_access_counter
  IceModelVec2S::begin_access();

}

void IceModelVec2T::end_access() const {
  // this call will decrement the m_access_counter
  IceModelVec2S::end_access();

  if (m_access_counter == 0) {
    PetscErrorCode ierr = DMDAVecRestoreArrayDOF(*m_da3, m_v3, &m_array3);
    PISM_CHK(ierr, "DMDAVecRestoreArrayDOF");
    m_array3 = NULL;
  }
}

void IceModelVec2T::init(const std::string &fname, unsigned int period, double reference_time) {

  const Logger &log = *m_grid->ctx()->log();

  m_filename       = fname;
  m_period         = period;
  m_reference_time = reference_time;

  // We find the variable in the input file and
  // try to find the corresponding time dimension.

  PIO nc(m_grid->com, "guess_mode", m_filename, PISM_READONLY);
  std::string name_found;
  bool exists, found_by_standard_name;
  nc.inq_var(m_metadata[0].get_name(), m_metadata[0].get_string("standard_name"),
             exists, name_found, found_by_standard_name);
  if (not exists) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "can't find %s (%s) in %s.",
                                  m_metadata[0].get_string("long_name").c_str(),
                                  m_metadata[0].get_name().c_str(),
                                  m_filename.c_str());
  }

  // find the time dimension:
  std::vector<std::string> dims;
  dims = nc.inq_vardims(name_found);
  
  std::string dimname = "";
  bool time_found = false;
  for (unsigned int i = 0; i < dims.size(); ++i) {
    dimname = dims[i];

    AxisType dimtype = nc.inq_dimtype(dimname, m_grid->ctx()->unit_system());

    if (dimtype == T_AXIS) {
      time_found = true;
      break;
    }
  }

  if (time_found) {
    // we're found the time dimension
    TimeseriesMetadata time_dimension(dimname, dimname, m_grid->ctx()->unit_system());

    time_dimension.set_string("units", m_grid->ctx()->time()->units_string());
    io::read_timeseries(nc, time_dimension,
                        *m_grid->ctx()->time(), log, m_time);

    std::string bounds_name = nc.get_att_text(dimname, "bounds");

    if (m_time.size() > 1) {
      if (not bounds_name.empty()) {
        // read time bounds data from a file
        TimeBoundsMetadata tb(bounds_name, dimname, m_grid->ctx()->unit_system());
        tb.set_string("units", time_dimension.get_string("units"));

        io::read_time_bounds(nc, tb, *m_grid->ctx()->time(),
                             log, m_time_bounds);

        // time bounds data overrides the time variable: we make t[j] be the
        // right end-point of the j-th interval
        for (unsigned int k = 0; k < m_time.size(); ++k) {
          m_time[k] = m_time_bounds[2*k + 1];
        }
      } else {
        // no time bounds attribute
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Variable '%s' does not have the time_bounds attribute.\n"
                                      "Cannot use time-dependent forcing data '%s' (%s) without time bounds.",
                                      dimname.c_str(),  m_metadata[0].get_string("long_name").c_str(),
                                      m_metadata[0].get_name().c_str());
      }
    } else {
      // only one time record; set fake time bounds:
      m_time_bounds.resize(2);
      m_time_bounds[0] = m_time[0] - 1;
      m_time_bounds[1] = m_time[0] + 1;
    }

  } else {
    // no time dimension; assume that we have only one record and set the time
    // to 0
    m_time.resize(1);
    m_time[0] = 0;

    // set fake time bounds:
    m_time_bounds.resize(2);
    m_time_bounds[0] = -1;
    m_time_bounds[1] =  1;
  }

  if (not is_increasing(m_time)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "times have to be strictly increasing (read from '%s').",
                                  m_filename.c_str());
  }

  if (m_period != 0) {
    if ((size_t)m_n_records < m_time.size()) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "buffer has to be big enough to hold all records of periodic data");
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
  m_time.resize(1);
  m_time[0] = 0;
  //N = 1 ;

  // set fake time bounds:
  m_time_bounds.resize(2);
  m_time_bounds[0] = -1;
  m_time_bounds[1] =  1;
}

//! Read some data to make sure that the interval (t, t + dt) is covered.
void IceModelVec2T::update(double t, double dt) {
  unsigned int m, n, last;

  if (m_filename.empty()) {
    // We are not reading data from a file.
    return;
  }

  if (m_time_bounds.size() == 0) {
    update(0);
    return;
  }

  if (m_period != 0) {
    // we read all data in IceModelVec2T::init() (see above)
    return;
  }

  if (m_N > 0) {
    last = m_first + (m_N - 1);

    // find the interval covered by data held in memory:
    double t0 = m_time_bounds[m_first * 2],
      t1 = m_time_bounds[last * 2 + 1];

    // just return if we have all the data we need:
    if (t >= t0 && t + dt <= t1) {
      return;
    }
  }

  // i will point to the smallest t so that t >= t 
  auto i = lower_bound(m_time_bounds.begin(), m_time_bounds.end(), t);

  // j will point to the smallest t so that t >= t + dt
  auto j = lower_bound(m_time_bounds.begin(), m_time_bounds.end(), t + dt);

  // some of the the interval (t, t + dt) is outside of the
  // time interval available in the file (on the right)
  // 
  if (i == m_time_bounds.end()) {
    m = (int)(m_time.size() - 1);
  } else {
    m = (int)((i - m_time_bounds.begin() - 1) / 2);
  }

  if (j == m_time_bounds.end()) {
    n = (int)(m_time.size() - 1);
  } else {
    n = (int)((j - m_time_bounds.begin() - 1) / 2);
  }

  // check if all the records necessary to cover this interval fit in the
  // buffer:
  if (n - m + 1 > m_n_records) {
    throw RuntimeError(PISM_ERROR_LOCATION, "IceModelVec2T::update(): timestep is too big");
  }

  update(m);
}

//! Update by reading at most n_records records from the file.
void IceModelVec2T::update(unsigned int start) {

  unsigned int time_size = (int)m_time.size();

  if (start >= time_size) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "IceModelVec2T::update(int start): start = %d is invalid", start);
  }

  unsigned int missing = std::min(m_n_records, time_size - start);

  if (start == static_cast<unsigned int>(m_first)) {
    // nothing to do
    return;
  }

  int kept = 0;
  if (m_first >= 0) {
    unsigned int last = m_first + (m_N - 1);
    if ((m_N > 0) && (start >= (unsigned int)m_first) && (start <= last)) {
      int discarded = start - m_first;
      kept = last - start + 1;
      discard(discarded);
      missing -= kept;
      start += kept;
      m_first += discarded;
    } else {
      m_first = start;
    }
  } else {
    m_first = start;
  }

  if (missing <= 0) {
    return;
  }
  
  m_N = kept + missing;

  Time::ConstPtr t = m_grid->ctx()->time();

  Logger::ConstPtr log = m_grid->ctx()->log();
  if (this->n_records() > 1) {
    log->message(4,
               "  reading \"%s\" into buffer\n"
               "          (short_name = %s): %d records, time intervals (%s, %s) through (%s, %s)...\n",
               metadata().get_string("long_name").c_str(), m_name.c_str(), missing,
               t->date(m_time_bounds[start*2]).c_str(),
               t->date(m_time_bounds[start*2 + 1]).c_str(),
               t->date(m_time_bounds[(start + missing - 1)*2]).c_str(),
               t->date(m_time_bounds[(start + missing - 1)*2 + 1]).c_str());
    m_report_range = false;
  } else {
    m_report_range = true;
  }

  PIO nc(m_grid->com, "guess_mode", m_filename, PISM_READONLY);

  const bool allow_extrapolation = m_grid->ctx()->config()->get_boolean("grid.allow_extrapolation");

  for (unsigned int j = 0; j < missing; ++j) {
    {
      petsc::VecArray tmp_array(m_v);
      io::regrid_spatial_variable(m_metadata[0], *m_grid, nc, start + j, CRITICAL,
                                  m_report_range, allow_extrapolation,
                                  0.0, m_interpolation_type, tmp_array.get());
    }

    m_grid->ctx()->log()->message(5, " %s: reading entry #%02d, year %s...\n",
                                  m_name.c_str(),
                                  start + j,
                                  t->date(m_time[start + j]).c_str());

    set_record(kept + j);
  }
}

//! Discard the first N records, shifting the rest of them towards the "beginning".
void IceModelVec2T::discard(int number) {

  if (number == 0) {
    return;
  }

  m_N -= number;

  double ***a3 = get_array3();
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (unsigned int k = 0; k < m_N; ++k) {
      a3[j][i][k] = a3[j][i][k + number];
    }
  }
  end_access();
}

//! Sets the record number n to the contents of the (internal) Vec v.
void IceModelVec2T::set_record(int n) {

  double  **a2 = get_array();
  double ***a3 = get_array3();
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    a3[j][i][n] = a2[j][i];
  }
  end_access();
  end_access();
}

//! Sets the (internal) Vec v to the contents of the nth record.
void IceModelVec2T::get_record(int n) {

  double  **a2 = get_array();
  double ***a3 = get_array3();
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    a2[j][i] = a3[j][i][n];
  }
  end_access();
  end_access();
}

//! \brief Given the time t and the current selected time-step dt,
//! determines the maximum possible time-step this IceModelVec2T allows.
/*!
  Returns -1 if any time step is OK at t.
 */
MaxTimestep IceModelVec2T::max_timestep(double t) const {
  // only allow going to the next record
  auto l = upper_bound(m_time_bounds.begin(), m_time_bounds.end(), t);
  if (l != m_time_bounds.end()) {
    double tmp = *l - t;

    if (tmp > 1.0) {                // never take time-steps shorter than 1 second
      return MaxTimestep(tmp);
    } else if ((l + 1) != m_time_bounds.end() and
               (l + 2) != m_time_bounds.end()) {
      return MaxTimestep(*(l + 2) - *l);
    } else {
      return MaxTimestep();
    }
  } else {
    return MaxTimestep();
  }

}

/*
 * \brief Use piecewise-constant interpolation to initialize
 * IceModelVec2T with the value at time `t`.
 *
 * \note This method does not check if an update() call is necessary!
 *
 * @param[in] t requested time
 *
 */
void IceModelVec2T::interp(double t) {

  std::vector<double> t_vector(1);
  t_vector[0] = t;
  init_interpolation(t_vector);

  get_record(m_interp_indices[0]);
}


/** 
 * Compute the average value over the time interval `[t, t + dt]`.
 *
 * @param t  start of the time interval, in seconds
 * @param dt length of the time interval, in seconds
 *
 */
void IceModelVec2T::average(double t, double dt) {

  double dt_years = units::convert(m_grid->ctx()->unit_system(),
                                   dt, "seconds", "years"); // *not* time->year(dt)

  // if only one record, nothing to do
  if (m_time.size() == 1) {
    return;
  }

  // Determine the number of small time-steps to use for averaging:
  int M = (int) ceil(m_n_evaluations_per_year * (dt_years));
  if (M < 1) {
    M = 1;
  }

  std::vector<double> ts(M);
  double ts_dt = dt / M;
  for (int k = 0; k < M; k++) {
    ts[k] = t + k * ts_dt;
  }

  init_interpolation(ts);

  double **a2 = get_array();         // calls begin_access()
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    a2[j][i] = average(i, j);
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

  assert(m_first >= 0);

  unsigned int
    index = 0,
    last  = m_first + m_N - 1;

  size_t ts_length = ts.size();

  // Compute "periodized" times if necessary.
  std::vector<double> times_requested(ts_length);
  if (m_period != 0) {
    for (unsigned int k = 0; k < ts_length; ++k) {
      times_requested[k] = m_grid->ctx()->time()->mod(ts[k] - m_reference_time, m_period);
    }
  } else {
    for (unsigned int k = 0; k < ts_length; ++k) {
      times_requested[k] = ts[k];
    }
  }

  m_interp_indices.resize(ts_length);

  if (m_time_bounds.size() == 0) {
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
    if (times_requested[k] < m_time_bounds[2*m_first]) {
      m_interp_indices[k] = 0;
      continue;
    }

    // extrapolate on the right:
    if (times_requested[k] >= m_time_bounds[2*last + 1]) {
      m_interp_indices[k] = m_N - 1;
      continue;
    }

    while (index < m_N) {
      if (m_time_bounds[2*(m_first + index) + 0] <= times_requested[k] &&
          m_time_bounds[2*(m_first + index) + 1] >  times_requested[k]) {
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
  double ***a3 = (double***) m_array3;
  unsigned int ts_length = m_interp_indices.size();

  for (unsigned int k = 0; k < ts_length; ++k) {
    result[k] = a3[j][i][m_interp_indices[k]];
  }
}

//! \brief Finds the average value at i,j over the interval (t, t +
//! dt) using the rectangle rule.
/*!
  Can (and should) be optimized. Later, though.
 */
double IceModelVec2T::average(int i, int j) {
  unsigned int M = m_interp_indices.size();
  double result = 0.0;

  if (m_N == 1) {
    double ***a3 = (double***) m_array3;
    result = a3[j][i][0];
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
  return result;
}



} // end of namespace pism

// Copyright (C) 2009--2024 Constantine Khroulev
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
#include <cassert>
#include <cmath>                // std::floor
#include <array>

#include "pism/util/array/Forcing.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Time.hh"
#include "pism/util/Grid.hh"
#include "pism/util/ConfigInterface.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Logger.hh"
#include "pism/util/interpolation.hh"
#include "pism/util/Context.hh"
#include "pism/util/array/Array_impl.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {
namespace array {

struct Forcing::Data {
  Data()
    : array(nullptr),
      first(-1),
      n_records(0),
      period(0.0),
      period_start(0.0) {
    // empty
  }
  //! all the times available in filename
  std::vector<double> time;

  //! the range of times covered by data in `filename`
  std::array<double,2> time_range;

  //! name of the file to read (regrid) from
  std::string filename;

  //! DM with dof equal to buffer_size
  std::shared_ptr<petsc::DM> da;

  //! a 3D Vec used to store records
  petsc::Vec v;

  //! pointer used to access records stored in memory
  double ***array;

  //! maximum number of records stored in memory
  unsigned int buffer_size;

  //! in-file index of the first record stored in memory (a signed `int` to allow
  //! `first==-1` as an *invalid* first value)
  int first;

  //! number of records currently kept in memory
  unsigned int n_records;

  //! temporal interpolation type
  InterpolationType interp_type;

  //! temporal interpolation code
  std::shared_ptr<Interpolation> interp;

  //! forcing period, in seconds
  double period;

  //! start of the period, in seconds
  double period_start;

  //! minimum time step length in max_timestep(), in seconds
  double dt_min;
};

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
 * @param[in] periodic true if this forcing field should be interpreted as periodic
 * @param[in] interpolation_type type of temporal interpolation (LINEAR or PIECEWISE_CONSTANT)
 */
Forcing::Forcing(std::shared_ptr<const Grid> grid,
                 const File &file,
                 const std::string &short_name,
                 const std::string &standard_name,
                 unsigned int max_buffer_size,
                 bool periodic,
                 InterpolationType interpolation_type)
  : array::Scalar(grid, short_name, 0),
    m_data(new Data()) {

  unsigned int n_records = file.nrecords(short_name, standard_name,
                                         grid->ctx()->unit_system());

  unsigned int buffer_size = 0;
  if (periodic) {
    // In the periodic case we try to keep all the records in RAM.
    buffer_size = n_records;

    if (interpolation_type == LINEAR) {
      // add two more records at the beginning and the end of the period to simplify
      // interpolation in time
      buffer_size += 2;
    }
  } else {
    buffer_size = std::min(n_records, max_buffer_size);
  }

  // Allocate storage for one record if the variable was not found. This is needed to be
  // able to cheaply allocate and then discard an "-atmosphere given" model
  // (atmosphere::Given) when "-surface given" (Given) is selected.
  buffer_size = std::max(buffer_size, 1U);

  allocate(buffer_size, interpolation_type);
}

std::shared_ptr<Forcing> Forcing::Constant(std::shared_ptr<const Grid> grid,
                                           const std::string &short_name,
                                           double value) {
  // note: cannot use std::make_shared because of a private constructor
  std::shared_ptr<Forcing> result(new Forcing(grid, short_name, 1,
                                              PIECEWISE_CONSTANT));

  // set constant value everywhere
  result->set(value);
  result->set_record(0);

  // set the time to zero
  result->m_data->time = {0.0};
  result->m_data->n_records = 1;
  result->m_data->first = 0;

  // set fake time bounds:
  double eps = 0.5 * result->m_data->dt_min;
  result->m_data->time_range = {-eps, eps};

  return result;
}

Forcing::Forcing(std::shared_ptr<const Grid> grid, const std::string &short_name,
                 unsigned int buffer_size,
                 InterpolationType interpolation_type)
  : array::Scalar(grid, short_name, 0),
    m_data(new Data()) {
  allocate(buffer_size, interpolation_type);
}

void Forcing::allocate(unsigned int buffer_size, InterpolationType interpolation_type) {
  m_impl->report_range = false;

  m_data->interp_type = interpolation_type;
  m_data->buffer_size = buffer_size;

  auto config = m_impl->grid->ctx()->config();

  m_data->dt_min = config->get_number("time_stepping.resolution");

  if (not (m_data->interp_type == PIECEWISE_CONSTANT or
           m_data->interp_type == LINEAR)) {
    throw RuntimeError(PISM_ERROR_LOCATION, "unsupported interpolation type");
  }

  // initialize the m_data->da member:
  m_data->da = m_impl->grid->get_dm(buffer_size, this->m_impl->da_stencil_width);

  // allocate the 3D Vec:
  PetscErrorCode ierr = DMCreateGlobalVector(*m_data->da, m_data->v.rawptr());
  PISM_CHK(ierr, "DMCreateGlobalVector");
}

Forcing::~Forcing() {
  delete m_data;
}

unsigned int Forcing::buffer_size() {
  return m_data->buffer_size;
}

double*** Forcing::array3() {
  return m_data->array;
}

void Forcing::begin_access() const {
  if (m_impl->access_counter == 0) {
    PetscErrorCode ierr = DMDAVecGetArrayDOF(*m_data->da, m_data->v, &m_data->array);
    PISM_CHK(ierr, "DMDAVecGetArrayDOF");
  }

  // this call will increment the m_access_counter
  array::Scalar::begin_access();
}

void Forcing::end_access() const {
  // this call will decrement the m_access_counter
  array::Scalar::end_access();

  if (m_impl->access_counter == 0) {
    PetscErrorCode ierr = DMDAVecRestoreArrayDOF(*m_data->da, m_data->v, &m_data->array);
    PISM_CHK(ierr, "DMDAVecRestoreArrayDOF");
    m_data->array = nullptr;
  }
}

void Forcing::init(const std::string &filename, bool periodic) {
  try {
    auto ctx = m_impl->grid->ctx();
    auto time = ctx->time();

    m_data->filename = filename;

    File file(m_impl->grid->com, m_data->filename, io::PISM_GUESS, io::PISM_READONLY);
    auto var = file.find_variable(m_impl->metadata[0].get_name(),
                                  m_impl->metadata[0]["standard_name"]);
    if (not var.exists) {
      throw RuntimeError(PISM_ERROR_LOCATION, "variable not found");
    }

    auto time_name = io::time_dimension(ctx->unit_system(), file, var.name);

    // dimension_length() will return 0 if a dimension is missing
    bool one_record = file.dimension_length(time_name) < 2;

    if (not one_record) {
      std::vector<double> times{};
      std::vector<double> bounds{};
      io::read_time_info(*ctx->log(), ctx->unit_system(),
                         file, time_name, time->units_string(),
                         times, bounds);

      if (periodic) {
        m_data->period_start = bounds.front();
        m_data->period = bounds.back() - bounds.front();
      }

      if (m_data->interp_type == PIECEWISE_CONSTANT) {
        // Time bounds data overrides the time variable: we make t[j] be the
        // left end-point of the j-th interval
        for (unsigned int k = 0; k < times.size(); ++k) {
          times[k] = bounds[2*k + 0];
        }
      }

      m_data->time       = times;
      m_data->time_range = {bounds.front(), bounds.back()};

      bool extrapolate = ctx->config()->get_flag("input.forcing.time_extrapolation");
      if (not (extrapolate or periodic)) {
        check_forcing_duration(*time, bounds.front(), bounds.back());
      }

    } else {
      // Only one time record or no time dimension at all: set fake time bounds assuming
      // that the user wants to use constant-in-time forcing for the whole simulation

      // this value does not matter
      m_data->time = {0.0};

      // set fake time bounds:
      double eps = 0.5 * m_data->dt_min;
      m_data->time_range = {-eps, eps};

      // note that in this case all data is periodic and constant in time
      m_data->period       = 0.0;
      m_data->period_start = 0.0;
    }

    if (periodic) {
      // read periodic data right away (we need to hold it all in memory anyway)
      init_periodic_data(file);
    }
  } catch (RuntimeError &e) {
    e.add_context("reading %s (%s) from '%s'",
                  m_impl->metadata[0].get_string("long_name").c_str(),
                  m_impl->metadata[0].get_name().c_str(),
                  m_data->filename.c_str());
    throw;
  }
}

/*!
 * Read all periodic data from the file and add two more records to simplify
 * interpolation.
 */
void Forcing::init_periodic_data(const File &file) {

  auto ctx = grid()->ctx();

  auto name = get_name();
  auto n_records = file.nrecords(name, metadata()["standard_name"],
                                 ctx->unit_system());

  auto buffer_required = n_records + 2 * static_cast<int>(m_data->interp_type == LINEAR);

  if (m_data->buffer_size < buffer_required) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "the buffer is too small to contain periodic data");
  }

  int offset = m_data->interp_type == LINEAR ? 1 : 0;

  // Read all the records and store them. The index offset leaves room for an extra record
  // needed to simplify interpolation
  auto variable = m_impl->metadata[0];
  auto V = file.find_variable(variable.get_name(), variable["standard_name"]);

  grid::InputGridInfo input_grid(file, V.name, variable.unit_system(), grid()->registration());

  LocalInterpCtx lic(input_grid, *grid(), levels(), m_impl->interpolation_type);

  for (unsigned int j = 0; j < n_records; ++j) {
    {
      lic.start[T_AXIS] = (int)j;
      lic.count[T_AXIS] = 1;

      petsc::VecArray tmp_array(vec());
      io::regrid_spatial_variable(variable, *grid(), lic, file, tmp_array.get());
    }

    auto time = ctx->time();
    auto log  = ctx->log();
    log->message(5, " %s: reading entry #%02d, time %s...\n", name.c_str(), j,
                 time->date(m_data->time[j]).c_str());

    set_record(offset + j);
  }

  m_data->n_records = buffer_required;
  m_data->first     = 0;

  if (m_data->interp_type == PIECEWISE_CONSTANT) {
    return;
  }

  double t0 = m_data->time_range[0];
  double t1 = m_data->time_range[1];

  // compute the interpolation factor used to find the value at the beginning of the
  // period
  double alpha = 0.0;
  {
    double dt1 = m_data->time.front() - t0;
    double dt2 = t1 - m_data->time.back();

    alpha = dt1 + dt2 > 0 ? dt2 / (dt1 + dt2) : 0.0;
  }

  // indexes used to access the first and the last entry in the buffer
  int first = 1;
  int last  = buffer_required - 2;

  array::AccessScope list{ this };

  // compute values at the beginning (and so at the end) of the period
  double **a2  = array();
  double ***a3 = array3();
  for (auto p = grid()->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    a2[j][i] = (1.0 - alpha) * a3[j][i][last] + alpha * a3[j][i][first];
  }

  set_record(0);
  set_record(m_data->buffer_size - 1);

  // add two more records to m_data->time
  {
    if (m_data->time[0] - t0 > 0.0) {
      m_data->time.insert(m_data->time.begin(), t0);
    } else {
      // The first time record is at the beginning of the time interval covered by
      // forcing. This means that the first record we added (set_record(0) above) is the
      // same as the second one (at index 1). Only one of them is needed.
      //
      // Here we use an arbitrary time *before* the beginning of the forcing time interval
      // to ensure that m_data->time is strictly increasing.
      //
      // Note: this time will not be used.
      const double dt = 1.0; // arbitrary; could be any positive number
      m_data->time.insert(m_data->time.begin(), t0 - dt);
    }
    if (t1 - m_data->time.back() > 0.0) {
      m_data->time.push_back(t1);
    } else {
      // The last time record is at the end of the time interval covered by forcing. This
      // means that the last record we added (set_record(m_data->buffer_size - 1) above)
      // is the same as the one before it. Only one of them is needed.
      //
      // Here we use an arbitrary time *after* the end of the forcing time interval to
      // ensure that m_data->time is strictly increasing.
      //
      // Note: this time will not be used.
      const double dt = 1.0; // arbitrary; could be any positive number
      m_data->time.push_back(t1 + dt);
    }
  }
}

//! Read some data to make sure that the interval (t, t + dt) is covered.
void Forcing::update(double t, double dt) {

  if (m_data->filename.empty()) {
    // We are not reading data from a file.
    return;
  }

  if (m_data->period > 0.0) {
    // we read all data in Forcing::init() (see above)
    return;
  }

  if (m_data->n_records > 0) {
    // in-file index of the last record currently in memory
    unsigned int last = m_data->first + (m_data->n_records - 1);

    double t0{}, t1{};
    if (m_data->interp_type == LINEAR) {
      t0 = m_data->time[m_data->first];
      t1 = m_data->time[last];
    } else {
      // piece-wise constant
      t0 = m_data->time[m_data->first];
      if (last + 1 < m_data->time.size()) {
        t1 = m_data->time[last + 1];
      } else {
        t1 = m_data->time_range[1];
      }
    }

    // just return if we have all the data we need:
    if (t >= t0 and t + dt <= t1) {
      return;
    }
  }

  Interpolation I(m_data->interp_type, m_data->time, { t, t + dt });

  unsigned int first = I.left(0), last = I.right(1), N = last - first + 1;

  // check if all the records necessary to cover this interval fit in the
  // buffer:
  if (N > m_data->buffer_size) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot read %d records of %s (buffer size: %d)", N,
                                  m_impl->name.c_str(), m_data->buffer_size);
  }

  update(first);
}

//! Update by reading at most buffer_size records from the file.
void Forcing::update(unsigned int start) {

  unsigned int time_size = (int)m_data->time.size();

  if (start >= time_size) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Forcing::update(int start): start = %d is invalid", start);
  }

  unsigned int missing = std::min(m_data->buffer_size, time_size - start);

  if (start == static_cast<unsigned int>(m_data->first)) {
    // nothing to do
    return;
  }

  int kept = 0;
  if (m_data->first >= 0) {
    unsigned int last = m_data->first + (m_data->n_records - 1);
    if ((m_data->n_records > 0) && (start >= (unsigned int)m_data->first) && (start <= last)) {
      int discarded = start - m_data->first;
      kept          = last - start + 1;
      discard(discarded);
      missing -= kept;
      start += kept;
      m_data->first += discarded;
    } else {
      m_data->first = start;
    }
  } else {
    m_data->first = start;
  }

  if (missing <= 0) {
    return;
  }

  m_data->n_records = kept + missing;

  auto t = m_impl->grid->ctx()->time();

  auto log = m_impl->grid->ctx()->log();
  if (this->buffer_size() > 1) {
    log->message(4,
                 "  reading \"%s\" into buffer\n"
                 "          (short_name = %s): %d records, time %s through %s...\n",
                 metadata().get_string("long_name").c_str(), m_impl->name.c_str(), missing,
                 t->date(m_data->time[start]).c_str(),
                 t->date(m_data->time[start + missing - 1]).c_str());
    m_impl->report_range = false;
  } else {
    m_impl->report_range = true;
  }

  File file(m_impl->grid->com, m_data->filename, io::PISM_GUESS, io::PISM_READONLY);

  auto variable = m_impl->metadata[0];

  try {
    auto V = file.find_variable(variable.get_name(), variable["standard_name"]);
    grid::InputGridInfo input_grid(file, V.name, variable.unit_system(), grid()->registration());

    LocalInterpCtx lic(input_grid, *grid(), levels(), m_impl->interpolation_type);

    for (unsigned int j = 0; j < missing; ++j) {
      lic.start[T_AXIS] = (int)(start + j);
      lic.count[T_AXIS] = 1;

      petsc::VecArray tmp_array(vec());
      io::regrid_spatial_variable(variable, *grid(), lic, file, tmp_array.get());

      log->message(5, " %s: reading entry #%02d, year %s...\n", m_impl->name.c_str(), start + j,
                   t->date(m_data->time[start + j]).c_str());

      set_record(kept + j);
    }
  } catch (RuntimeError &e) {
    e.add_context("regridding '%s' from '%s'", this->get_name().c_str(), m_data->filename.c_str());
    throw;
  }
}

//! Discard the first N records, shifting the rest of them towards the "beginning".
void Forcing::discard(int number) {

  if (number == 0) {
    return;
  }

  m_data->n_records -= number;

  array::AccessScope l{this};

  double ***a3 = array3();
  for (auto p = m_impl->grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (unsigned int k = 0; k < m_data->n_records; ++k) {
      a3[j][i][k] = a3[j][i][k + number];
    }
  }
}

//! Sets the record number n to the contents of the (internal) Vec v.
void Forcing::set_record(int n) {

  array::AccessScope l{this};

  double  **a2 = array();
  double ***a3 = array3();
  for (auto p = m_impl->grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();
    a3[j][i][n] = a2[j][i];
  }
}

//! @brief Given the time t determines the maximum possible time-step this Forcing
//! allows.
MaxTimestep Forcing::max_timestep(double t) const {
  auto time_size = m_data->time.size();

  if (m_data->period > 0.0) {
    // all periodic data is stored in RAM and there is no time step restriction
    return {};
  }

  if (t >= m_data->time.back()) {
    // Reached the end of forcing: no time step restriction. We will need only one record
    // to use constant extrapolation and it will surely fit in the buffer.
    return {};
  }

  // find the index k such that m_data->time[k] <= x < m_data->time[k + 1]
  // Note: `L` below will be strictly less than `time_size - 1`.
  size_t L = gsl_interp_bsearch(m_data->time.data(), t, 0, time_size - 1);

  // find the index of the last record we could read in given the size of the buffer
  size_t R = L + m_data->buffer_size - 1;

  if (R >= time_size - 1) {
    // We can read all the remaining records: no time step restriction from now on
    return {};
  }

  if (m_data->interp_type == PIECEWISE_CONSTANT) {
    // in the piece-wise constant case we can go all the way to the *next* record
    R = std::min(R + 1, time_size - 1);
    return m_data->time[R] - t;
  }

  return m_data->time[R] - t;
}

/*
 * @brief Initialize Forcing with the value at time `t`.
 *
 * @note This method does not check if an update() call is necessary!
 *
 * @param[in] t requested time
 *
 */
void Forcing::interp(double t) {

  init_interpolation({t});

  // There is only one point to interpolate at ("t" above). Here we get interpolation
  // indexes and the corresponding weight. Here L == R for the piecewise-constant
  // interpolation.
  int L = m_data->interp->left(0);
  int R = m_data->interp->right(0);
  double alpha = m_data->interp->alpha(0);

  array::AccessScope l{this};
  double ***a3 = array3();
  double  **a2 = array();

  for (auto p = m_impl->grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();
    auto *column = a3[j][i];
    // result (LHS) is a weighted average of two values.
    a2[j][i] = column[L] + alpha * (column[R] - column[L]);
  }
}


/**
 * Compute the average value over the time interval `[t, t + dt]`.
 *
 * @param t  start of the time interval, in seconds
 * @param dt length of the time interval, in seconds
 *
 */
void Forcing::average(double t, double dt) {

  // if only one record, nothing to do
  if (m_data->time.size() == 1 or
      m_data->n_records == 1 or
      dt == 0.0) {
    interp(t);
    return;
  }

  const double *data = &m_data->time[m_data->first];
  size_t data_size = m_data->n_records;
  auto type = m_data->interp_type;

  std::map<size_t, double> weights{};

  if (m_data->period > 0.0) {
    double a = t;
    double b = t + dt;
    double t0 = m_data->period_start;
    double P = m_data->period;

    // N_periods is the number of complete periods in the *middle* of the integration
    // interval
    //
    // Note that the total number of complete periods is equal to (N_periods + 1) *if*
    // delta == 0.0.
    double N_periods = 0.0;
    double delta = 0.0;
    double gamma = 0.0;
    {
      double N = std::floor((a - t0) / P);
      double M = std::floor((b - t0) / P);

      N_periods = M - (N + 1);
      delta = (a - t0) - P * N;
      gamma = (b - t0) - P * M;
    }

    if (N_periods >= 0.0) {
      assert(t0 + delta < t0 + P);
      auto W1 = integration_weights(data, data_size, type, t0 + delta, t0 + P);

      std::map<size_t, double> W2{};
      if (N_periods > 0) {
        // note: we know that t0 < t0 + P because P > 0
        W2 = integration_weights(data, data_size, type, t0, t0 + P);
      } else {
        W2 = {};
      }

      std::map<size_t, double> W3{};
      if (gamma > 0.0) {
        W3 = integration_weights(data, data_size, type, t0, t0 + gamma);
      } else {
        W3 = {};
      }

      // before the first complete period:
      weights = W1;
      // an integer number of complete periods:
      for (const auto &w : W2) {
        weights[w.first] += N_periods * w.second;
      }
      // after the last complete period:
      for (const auto &w : W3) {
        weights[w.first] += w.second;
      }
    } else {
      assert(t0 + delta < t0 + gamma);
      weights = integration_weights(data, data_size, type, t0 + delta, t0 + gamma);
    }
  } else {
    assert(dt > 0.0);
    weights = integration_weights(data, data_size, type, t, t + dt);
  }

  array::AccessScope l{this};
  double **a2 = array();
  double ***a3 = array3();

  for (auto p = m_impl->grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    a2[j][i] = 0.0;
    for (const auto &weight : weights) {
      size_t k = weight.first;
      double w = weight.second;
      a2[j][i] += w * a3[j][i][k];
    }
    a2[j][i] /= dt;
  }
}

/**
 * \brief Compute weights for the piecewise-constant interpolation.
 * This is used *both* for time-series and "snapshots".
 *
 * @param ts requested times, in seconds
 *
 */
void Forcing::init_interpolation(const std::vector<double> &ts) {

  assert(m_data->first >= 0);

  // Compute "periodized" times if necessary.
  std::vector<double> times_requested(ts.size());
  if (m_data->period > 0.0) {
    double P  = m_data->period;
    double T0 = m_data->period_start;

    for (unsigned int k = 0; k < ts.size(); ++k) {
      double t = ts[k] - T0;

      t -= std::floor(t / P) * P;

      times_requested[k] = T0 + t;
    }
  } else {
    times_requested = ts;
  }

  m_data->interp.reset(new Interpolation(m_data->interp_type,
                                         &m_data->time[m_data->first],
                                         m_data->n_records,
                                         times_requested.data(),
                                         times_requested.size()));
}

/**
 * \brief Compute values of the time-series using precomputed indices
 * (and piece-wise constant or piece-wise linear interpolation).
 *
 * @param i,j map-plane grid point
 * @param result pointer to an allocated array of the size matching the one passed to
 *               init_interpolation()
 *
 */
void Forcing::interp(int i, int j, std::vector<double> &result) {
  double ***a3 = array3();

  result.resize(m_data->interp->alpha().size());

  m_data->interp->interpolate(a3[j][i], result.data());
}

} // end of namespace array
} // end of namespace pism

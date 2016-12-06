// Copyright (C) 2008-2016 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

// Implementation of the atmosphere model using constant-in-time precipitation
// and a cosine yearly cycle for near-surface air temperatures.

#include <gsl/gsl_math.h>       // M_PI, GSL_NAN

#include "PAYearlyCycle.hh"
#include "base/util/PISMTime.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace atmosphere {

YearlyCycle::YearlyCycle(IceGrid::ConstPtr g)
  : AtmosphereModel(g) {

  m_snow_temp_july_day = m_config->get_double("atmosphere.fausto_air_temp.summer_peak_day");

  // Allocate internal IceModelVecs:
  m_air_temp_mean_annual.create(m_grid, "air_temp_mean_annual", WITHOUT_GHOSTS);
  m_air_temp_mean_annual.set_attrs("diagnostic",
                                   "mean annual near-surface air temperature (without sub-year time-dependence or forcing)",
                                   "K",
                                   "");  // no CF standard_name ??
  m_air_temp_mean_annual.metadata().set_string("source", m_reference);

  m_air_temp_mean_july.create(m_grid, "air_temp_mean_july", WITHOUT_GHOSTS);
  m_air_temp_mean_july.set_attrs("diagnostic",
                                 "mean July near-surface air temperature (without sub-year time-dependence or forcing)",
                                 "Kelvin",
                                 "");  // no CF standard_name ??
  m_air_temp_mean_july.metadata().set_string("source", m_reference);

  m_precipitation.create(m_grid, "precipitation", WITHOUT_GHOSTS);
  m_precipitation.set_attrs("model_state", "precipitation rate",
                                "kg m-2 second-1", "", 0);
  m_precipitation.metadata(0).set_string("glaciological_units", "kg m-2 year-1");
  m_precipitation.write_in_glaciological_units = true;
  m_precipitation.set_time_independent(true);
}

YearlyCycle::~YearlyCycle() {
  // empty
}

//! Reads in the precipitation data from the input file.
void YearlyCycle::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  InputOptions opts = process_input_options(m_grid->com);
  init_internal(opts.filename, opts.type == INIT_BOOTSTRAP, opts.record);
}

//! Read precipitation data from a given file.
void YearlyCycle::init_internal(const std::string &input_filename, bool do_regrid,
                                unsigned int start) {
  // read precipitation rate from file
  m_log->message(2,
             "    reading mean annual ice-equivalent precipitation rate 'precipitation'\n"
             "      from %s ... \n",
             input_filename.c_str());
  if (do_regrid == true) {
    m_precipitation.regrid(input_filename, CRITICAL); // fails if not found!
  } else {
    m_precipitation.read(input_filename, start); // fails if not found!
  }
}

void YearlyCycle::define_model_state_impl(const PIO &output) const {
  m_precipitation.define(output);
}

void YearlyCycle::write_model_state_impl(const PIO &output) const {
  m_precipitation.write(output);
}

//! Copies the stored precipitation field into result.
void YearlyCycle::mean_precipitation_impl(IceModelVec2S &result) const {
  result.copy_from(m_precipitation);
}

//! Copies the stored mean annual near-surface air temperature field into result.
void YearlyCycle::mean_annual_temp_impl(IceModelVec2S &result) const {
  result.copy_from(m_air_temp_mean_annual);
}

//! Copies the stored mean July near-surface air temperature field into result.
void YearlyCycle::mean_july_temp(IceModelVec2S &result) const {
  result.copy_from(m_air_temp_mean_july);
}

void YearlyCycle::init_timeseries_impl(const std::vector<double> &ts) const {
  // constants related to the standard yearly cycle
  const double
    julyday_fraction = m_grid->ctx()->time()->day_of_the_year_to_day_fraction(m_snow_temp_july_day);

  size_t N = ts.size();

  m_ts_times.resize(N);
  m_cosine_cycle.resize(N);
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    double tk = m_grid->ctx()->time()->year_fraction(ts[k]) - julyday_fraction;

    m_ts_times[k] = ts[k];
    m_cosine_cycle[k] = cos(2.0 * M_PI * tk);
  }
}

void YearlyCycle::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = m_precipitation(i,j);
  }
}

void YearlyCycle::temp_time_series_impl(int i, int j, std::vector<double> &result) const {

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] = m_air_temp_mean_annual(i,j) + (m_air_temp_mean_july(i,j) - m_air_temp_mean_annual(i,j)) * m_cosine_cycle[k];
  }
}

void YearlyCycle::begin_pointwise_access_impl() const {
  m_air_temp_mean_annual.begin_access();
  m_air_temp_mean_july.begin_access();
  m_precipitation.begin_access();
}

void YearlyCycle::end_pointwise_access_impl() const {
  m_air_temp_mean_annual.end_access();
  m_air_temp_mean_july.end_access();
  m_precipitation.end_access();
}
std::map<std::string, Diagnostic::Ptr> YearlyCycle::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = AtmosphereModel::diagnostics_impl();

  result["air_temp_mean_july"] = Diagnostic::Ptr(new PA_mean_july_temp(this));

  return result;
}

PA_mean_july_temp::PA_mean_july_temp(const YearlyCycle *m)
  : Diag<YearlyCycle>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "air_temp_mean_july")};

  set_attrs("mean July near-surface air temperature used in the cosine yearly cycle", "",
            "Kelvin", "Kelvin", 0);
}

IceModelVec::Ptr PA_mean_july_temp::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "air_temp_mean_july", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->mean_july_temp(*result);

  return result;
}

} // end of namespace atmosphere
} // end of namespace pism

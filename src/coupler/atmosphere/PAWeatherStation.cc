/* Copyright (C) 2014, 2015, 2016 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <gsl/gsl_math.h>

#include "PAWeatherStation.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/pism_const.hh"
#include "base/util/pism_options.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/PISMTime.hh"
#include "base/util/IceGrid.hh"
#include "base/util/io/PIO.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace atmosphere {

WeatherStation::WeatherStation(IceGrid::ConstPtr g)
  : AtmosphereModel(g),
    m_precipitation_timeseries(*g, "precipitation", m_config->get_string("time.dimension_name")),
    m_air_temp_timeseries(*g, "air_temp", m_config->get_string("time.dimension_name"))
{
  m_precipitation_timeseries.dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());
  m_precipitation_timeseries.metadata().set_string("units", "kg m-2 second-1");
  m_precipitation_timeseries.metadata().set_string("long_name",
                                        "ice-equivalent precipitation rate");

  m_air_temp_timeseries.dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());
  m_air_temp_timeseries.metadata().set_string("units", "Kelvin");
  m_air_temp_timeseries.metadata().set_string("long_name",
                                          "near-surface air temperature");
}

WeatherStation::~WeatherStation() {
  // empty
}

void WeatherStation::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2,
             "* Initializing the constant-in-space atmosphere model\n"
             "  for use with scalar data from one weather station\n"
             "  combined with lapse rate corrections...\n");

  std::string option = "-atmosphere_one_station_file";

  options::String filename(option,
                           "Specifies a file containing scalar time-series"
                           " 'precipitation' and 'air_temp'.");

  if (not filename.is_set()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Command-line option %s is required.", option.c_str());
  }

  m_log->message(2,
             "  - Reading air temperature and precipitation from '%s'...\n",
             filename->c_str());

  PIO nc(m_grid->com, "netcdf3", filename, PISM_READONLY);
  {
    m_precipitation_timeseries.read(nc, *m_grid->ctx()->time(), *m_grid->ctx()->log());
    m_air_temp_timeseries.read(nc, *m_grid->ctx()->time(), *m_grid->ctx()->log());
  }
}

MaxTimestep WeatherStation::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere weather_station");
}

void WeatherStation::update_impl(double t, double dt) {
  m_t = t;
  m_dt = dt;
}

void WeatherStation::mean_precipitation_impl(IceModelVec2S &result) const {
  const double one_week = 7 * 24 * 60 * 60;

  unsigned int N = (unsigned int)(ceil(m_dt / one_week)); // one point per week

  result.set(m_precipitation_timeseries.average(m_t, m_dt, N));
}

void WeatherStation::mean_annual_temp_impl(IceModelVec2S &result) const {
  const double one_week = 7 * 24 * 60 * 60;

  unsigned int N = (unsigned int)(ceil(m_dt / one_week)); // one point per week

  result.set(m_air_temp_timeseries.average(m_t, m_dt, N));
}

void WeatherStation::begin_pointwise_access_impl() const {
  // empty
}

void WeatherStation::end_pointwise_access_impl() const {
  // empty
}

void WeatherStation::init_timeseries_impl(const std::vector<double> &ts) const {
  size_t N = ts.size();

  m_precip_values.resize(N);
  m_air_temp_values.resize(N);

  for (unsigned int k = 0; k < N; ++k) {
    m_precip_values[k]   = m_precipitation_timeseries(ts[k]);
    m_air_temp_values[k] = m_air_temp_timeseries(ts[k]);
  }
}

void WeatherStation::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  (void)i;
  (void)j;

  result = m_precip_values;
}

void WeatherStation::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  (void)i;
  (void)j;

  result = m_air_temp_values;
}

} // end of namespace atmosphere
} // end of namespace pism

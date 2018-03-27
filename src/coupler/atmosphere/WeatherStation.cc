/* Copyright (C) 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "WeatherStation.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/Time.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace atmosphere {

WeatherStation::WeatherStation(IceGrid::ConstPtr grid)
  : AtmosphereModel(grid),
    m_precipitation_timeseries(*grid, "precipitation", m_config->get_string("time.dimension_name")),
    m_air_temp_timeseries(*grid, "air_temp", m_config->get_string("time.dimension_name"))
{
  m_precipitation_timeseries.dimension().set_string("units", grid->ctx()->time()->units_string());
  m_precipitation_timeseries.variable().set_string("units", "kg m-2 second-1");
  m_precipitation_timeseries.variable().set_string("long_name",
                                        "ice-equivalent precipitation rate");

  m_air_temp_timeseries.dimension().set_string("units", grid->ctx()->time()->units_string());
  m_air_temp_timeseries.variable().set_string("units", "Kelvin");
  m_air_temp_timeseries.variable().set_string("long_name",
                                          "near-surface air temperature");

  m_precipitation = allocate_precipitation(grid);
  m_temperature   = allocate_temperature(grid);
}

WeatherStation::~WeatherStation() {
  // empty
}

void WeatherStation::init_impl(const Geometry &geometry) {
  (void) geometry;

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

void WeatherStation::update_impl(const Geometry &geometry, double t, double dt) {
  (void) geometry;

  double one_week = 7 * 24 * 60 * 60;
  unsigned int N = (unsigned int)(ceil(dt / one_week)); // one point per week

  m_precipitation->set(m_precipitation_timeseries.average(t, dt, N));

  m_temperature->set(m_air_temp_timeseries.average(t, dt, N));
}

const IceModelVec2S& WeatherStation::mean_precipitation_impl() const {
  return *m_precipitation;
}

const IceModelVec2S& WeatherStation::mean_annual_temp_impl() const {
  return *m_temperature;
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

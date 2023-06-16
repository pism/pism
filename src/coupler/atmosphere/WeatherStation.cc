/* Copyright (C) 2014, 2015, 2016, 2017, 2018, 2020, 2021, 2022 PISM Authors
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
#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/ScalarForcing.hh"

namespace pism {
namespace atmosphere {

WeatherStation::WeatherStation(std::shared_ptr<const Grid> grid)
  : AtmosphereModel(grid) {

  m_log->message(2,
                 "* Using the constant-in-space atmosphere model\n"
                 "  for use with scalar data from one weather station\n"
                 "  combined with lapse rate corrections...\n");

  auto filename = m_config->get_string("atmosphere.one_station.file");

  if (filename.empty()) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "atmosphere.one_station.file cannot be empty");
  }

  m_log->message(2,
                 "  - Reading air temperature and precipitation from '%s'...\n",
                 filename.c_str());

  auto &ctx = *grid->ctx();

  bool periodic = false;

  m_precipitation_timeseries = std::make_shared<ScalarForcing>(ctx,
                                                               filename,
                                                               "precipitation",
                                                               "kg m-2 second-1",
                                                               "kg m-2 year-1",
                                                               "ice-equivalent precipitation rate",
                                                               periodic);

  m_air_temp_timeseries = std::make_shared<ScalarForcing>(ctx,
                                                          filename,
                                                          "air_temp",
                                                          "Kelvin",
                                                          "Kelvin",
                                                          "near-surface air temperature",
                                                          periodic);

  m_precipitation = allocate_precipitation(grid);
  m_temperature   = allocate_temperature(grid);
}

void WeatherStation::init_impl(const Geometry &geometry) {
  (void) geometry;
}

MaxTimestep WeatherStation::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere weather_station");
}

void WeatherStation::update_impl(const Geometry &geometry, double t, double dt) {
  (void) geometry;

  m_precipitation->set(m_precipitation_timeseries->average(t, dt));

  m_temperature->set(m_air_temp_timeseries->average(t, dt));
}

const array::Scalar& WeatherStation::precipitation_impl() const {
  return *m_precipitation;
}

const array::Scalar& WeatherStation::air_temperature_impl() const {
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
    m_precip_values[k]   = m_precipitation_timeseries->value(ts[k]);
    m_air_temp_values[k] = m_air_temp_timeseries->value(ts[k]);
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

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
    m_precipitation(*g, "precipitation", m_config->get_string("time_dimension_name")),
    m_air_temperature(*g, "air_temp", m_config->get_string("time_dimension_name")),
    m_precip_metadata(m_sys, "precipitation"),
    m_air_temp_metadata(m_sys, "air_temp")
{
  m_precipitation.dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());
  m_precipitation.metadata().set_string("units", "m second-1");
  m_precipitation.metadata().set_string("long_name",
                                            "ice-equivalent precipitation rate");

  m_air_temperature.dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());
  m_air_temperature.metadata().set_string("units", "Kelvin");
  m_air_temperature.metadata().set_string("long_name",
                                              "near-surface air temperature");

  m_air_temp_metadata.set_string("pism_intent", "diagnostic");
  m_air_temp_metadata.set_string("long_name", "near-surface air temperature");
  m_air_temp_metadata.set_string("units", "K");

  m_precip_metadata.set_string("pism_intent", "diagnostic");
  m_precip_metadata.set_string("long_name", "precipitation, units of ice-equivalent thickness per time");
  m_precip_metadata.set_string("units", "m second-1");
  m_precip_metadata.set_string("glaciological_units", "m year-1");
}

WeatherStation::~WeatherStation() {
  // empty
}

void WeatherStation::init() {

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
    throw RuntimeError::formatted("Command-line option %s is required.", option.c_str());
  }

  m_log->message(2,
             "  - Reading air temperature and precipitation from '%s'...\n",
             filename->c_str());

  PIO nc(m_grid->com, "netcdf3");
  nc.open(filename, PISM_READONLY);
  {
    m_precipitation.read(nc, *m_grid->ctx()->time(), *m_grid->ctx()->log());
    m_air_temperature.read(nc, *m_grid->ctx()->time(), *m_grid->ctx()->log());
  }
  nc.close();
}

MaxTimestep WeatherStation::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void WeatherStation::update_impl(double t, double dt) {
  m_t = t;
  m_dt = dt;
}

void WeatherStation::mean_precipitation(IceModelVec2S &result) {
  const double one_week = 7 * 24 * 60 * 60;

  unsigned int N = (unsigned int)(ceil(m_dt / one_week)); // one point per week

  result.set(m_precipitation.average(m_t, m_dt, N));
}

void WeatherStation::mean_annual_temp(IceModelVec2S &result) {
  const double one_week = 7 * 24 * 60 * 60;

  unsigned int N = (unsigned int)(ceil(m_dt / one_week)); // one point per week

  result.set(m_air_temperature.average(m_t, m_dt, N));
}

void WeatherStation::begin_pointwise_access() {
  // empty
}

void WeatherStation::end_pointwise_access() {
  // empty
}

void WeatherStation::init_timeseries(const std::vector<double> &ts) {
  size_t N = ts.size();

  m_precip_values.resize(N);
  m_air_temp_values.resize(N);

  for (unsigned int k = 0; k < N; ++k) {
    m_precip_values[k]   = m_precipitation(ts[k]);
    m_air_temp_values[k] = m_air_temperature(ts[k]);
  }
}

void WeatherStation::precip_time_series(int i, int j,
                                                    std::vector<double> &result) {
  (void)i;
  (void)j;

  result = m_precip_values;
}

void WeatherStation::temp_time_series(int i, int j,
                                                  std::vector<double> &result) {
  (void)i;
  (void)j;

  result = m_air_temp_values;
}

void WeatherStation::temp_snapshot(IceModelVec2S &result) {
  result.set(m_air_temperature(m_t + 0.5*m_dt));
}

void WeatherStation::add_vars_to_output_impl(const std::string &keyword,
                                               std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big" || keyword == "2dbig") {
    result.insert("air_temp");
    result.insert("precipitation");
  }
}

void WeatherStation::define_variables_impl(const std::set<std::string> &vars,
                                                  const PIO &nc, IO_Type nctype) {
  std::string order = m_config->get_string("output_variable_order");

  if (set_contains(vars, "air_temp")) {
    // don't write using glaciological units
    io::define_spatial_variable(m_air_temp_metadata, *m_grid, nc, nctype, order, false);
  }

  if (set_contains(vars, "precipitation")) {
    // do write using glaciological units
    io::define_spatial_variable(m_precip_metadata, *m_grid, nc, nctype, order, true);
  }
}

void WeatherStation::write_variables_impl(const std::set<std::string> &vars,
                                       const PIO &nc) {

  if (set_contains(vars, "air_temp")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "air_temp", WITHOUT_GHOSTS);
    tmp.metadata() = m_air_temp_metadata;

    mean_annual_temp(tmp);

    tmp.write(nc);
  }

  if (set_contains(vars, "precipitation")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "precipitation", WITHOUT_GHOSTS);
    tmp.metadata() = m_precip_metadata;

    mean_precipitation(tmp);

    tmp.write_in_glaciological_units = true;
    tmp.write(nc);
  }
}

} // end of namespace atmosphere
} // end of namespace pism

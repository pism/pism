/* Copyright (C) 2014 PISM Authors
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

#include "PAWeatherStation.hh"
#include "PISMConfig.hh"
#include "pism_const.hh"
#include "pism_options.hh"
#include "iceModelVec.hh"
#include "PISMTime.hh"

#include "error_handling.hh"

namespace pism {

PAWeatherStation::PAWeatherStation(IceGrid &g)
  : AtmosphereModel(g),
    m_precipitation(&g, "precipitation", g.config.get_string("time_dimension_name")),
    m_air_temperature(&g, "air_temp", g.config.get_string("time_dimension_name")),
    m_precip_metadata(g.config.get_unit_system(), "precipitation", m_grid),
    m_air_temp_metadata(g.config.get_unit_system(), "air_temp", m_grid)
{
  m_precipitation.get_dimension_metadata().set_units(m_grid.time->units_string());
  m_precipitation.get_metadata().set_units("m / second");
  m_precipitation.get_metadata().set_string("long_name",
                                            "ice-equivalent precipitation rate");

  m_air_temperature.get_dimension_metadata().set_units(m_grid.time->units_string());
  m_air_temperature.get_metadata().set_units("Kelvin");
  m_air_temperature.get_metadata().set_string("long_name",
                                              "near-surface air temperature");

  m_air_temp_metadata.set_string("pism_intent", "diagnostic");
  m_air_temp_metadata.set_string("long_name", "near-surface air temperature");
  m_air_temp_metadata.set_units("K");

  m_precip_metadata.set_string("pism_intent", "diagnostic");
  m_precip_metadata.set_string("long_name", "precipitation, units of ice-equivalent thickness per time");
  m_precip_metadata.set_units("m / s");
  m_precip_metadata.set_glaciological_units("m / year");
}

PAWeatherStation::~PAWeatherStation() {
  // empty
}

void PAWeatherStation::init(Vars &vars) {

  (void)vars;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  verbPrintf(2, m_grid.com,
             "* Initializing the constant-in-space atmosphere model\n"
             "  for use with scalar data from one weather station\n"
             "  combined with lapse rate corrections...\n");

  std::string filename,
    option = "-atmosphere_one_station_file";
  bool bc_file_set = false;

  {
    OptionsString(option,
                  "Specifies a file containing scalar time-series 'precipitation' and 'air_temp'.",
                  filename, bc_file_set);
  }

  if (bc_file_set == false) {
    throw RuntimeError::formatted("Command-line option %s is required.", option.c_str());
  }

  verbPrintf(2, m_grid.com,
             "  - Reading air temperature and precipitation from '%s'...\n",
             filename.c_str());

  PIO nc(m_grid.com, "netcdf3", m_grid.config.get_unit_system());
  nc.open(filename, PISM_READONLY);
  {
    m_precipitation.read(nc, m_grid.time);
    m_air_temperature.read(nc, m_grid.time);
  }
  nc.close();
}

void PAWeatherStation::update(double t, double dt) {
  m_t = t;
  m_dt = dt;
}

void PAWeatherStation::mean_precipitation(IceModelVec2S &result) {
  const double one_week = 7 * 24 * 60 * 60;

  unsigned int N = (unsigned int)(ceil(m_dt / one_week)); // one point per week

  result.set(m_precipitation.average(m_t, m_dt, N));
}

void PAWeatherStation::mean_annual_temp(IceModelVec2S &result) {
  const double one_week = 7 * 24 * 60 * 60;

  unsigned int N = (unsigned int)(ceil(m_dt / one_week)); // one point per week

  result.set(m_air_temperature.average(m_t, m_dt, N));
}

void PAWeatherStation::begin_pointwise_access() {
  // empty
}

void PAWeatherStation::end_pointwise_access() {
  // empty
}

void PAWeatherStation::init_timeseries(const std::vector<double> &ts) {
  size_t N = ts.size();

  m_precip_values.resize(N);
  m_air_temp_values.resize(N);

  for (unsigned int k = 0; k < N; ++k) {
    m_precip_values[k]   = m_precipitation(ts[k]);
    m_air_temp_values[k] = m_air_temperature(ts[k]);
  }
}

void PAWeatherStation::precip_time_series(int i, int j,
                                                    std::vector<double> &result) {
  (void)i;
  (void)j;

  result = m_precip_values;
}

void PAWeatherStation::temp_time_series(int i, int j,
                                                  std::vector<double> &result) {
  (void)i;
  (void)j;

  result = m_air_temp_values;
}

void PAWeatherStation::temp_snapshot(IceModelVec2S &result) {
  result.set(m_air_temperature(m_t + 0.5*m_dt));
}

void PAWeatherStation::add_vars_to_output(const std::string &keyword,
                                          std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big") {
    result.insert("air_temp");
    result.insert("precipitation");
  }
}

void PAWeatherStation::define_variables(const std::set<std::string> &vars,
                                                  const PIO &nc, IO_Type nctype) {
  if (set_contains(vars, "air_temp")) {
    // don't write using glaciological units
    m_air_temp_metadata.define(nc, nctype, false);
  }

  if (set_contains(vars, "precipitation")) {
    // do write using glaciological units
    m_precip_metadata.define(nc, nctype, true);
  }
}

void PAWeatherStation::write_variables(const std::set<std::string> &vars,
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

} // end of namespace pism

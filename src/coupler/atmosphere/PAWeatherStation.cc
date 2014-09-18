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

namespace pism {

PAWeatherStation::PAWeatherStation(IceGrid &g, const Config &conf)
  : AtmosphereModel(g, conf),
    m_precipitation(&g, "precipitation", conf.get_string("time_dimension_name")),
    m_air_temperature(&g, "air_temp", conf.get_string("time_dimension_name")),
    m_precip_metadata(g.get_unit_system()),
    m_air_temp_metadata(g.get_unit_system())
{
  PetscErrorCode ierr = allocate(); CHKERRCONTINUE(ierr);
}

PetscErrorCode PAWeatherStation::allocate() {
  PetscErrorCode ierr;

  m_precipitation.get_dimension_metadata().set_units(grid.time->units_string());
  m_precipitation.get_metadata().set_units("m / second");
  m_precipitation.get_metadata().set_string("long_name",
                                            "ice-equivalent precipitation rate");

  m_air_temperature.get_dimension_metadata().set_units(grid.time->units_string());
  m_air_temperature.get_metadata().set_units("Kelvin");
  m_air_temperature.get_metadata().set_string("long_name",
                                              "near-surface air temperature");

  m_air_temp_metadata.init_2d("air_temp", grid);
  m_air_temp_metadata.set_string("pism_intent", "diagnostic");
  m_air_temp_metadata.set_string("long_name", "near-surface air temperature");
  ierr = m_air_temp_metadata.set_units("K"); CHKERRQ(ierr);

  m_precip_metadata.init_2d("precipitation", grid);
  m_precip_metadata.set_string("pism_intent", "diagnostic");
  m_precip_metadata.set_string("long_name", "precipitation, units of ice-equivalent thickness per time");
  ierr = m_precip_metadata.set_units("m / s"); CHKERRQ(ierr);
  ierr = m_precip_metadata.set_glaciological_units("m / year"); CHKERRQ(ierr);

  return 0;
}

PAWeatherStation::~PAWeatherStation() {
  // empty
}

PetscErrorCode PAWeatherStation::init(Vars &vars) {
  PetscErrorCode ierr;

  (void)vars;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the constant-in-space atmosphere model\n"
                    "  for use with scalar data from one weather station\n"
                    "  combined with lapse rate corrections...\n"); CHKERRQ(ierr);

  std::string filename,
    option = "-atmosphere_one_station_file";
  bool bc_file_set = false;

  ierr = PetscOptionsBegin(grid.com, "", "Climate forcing options", ""); CHKERRQ(ierr);
  {
    ierr = OptionsString(option,
                         "Specifies a file containing scalar time-series 'precipitation' and 'air_temp'.",
                         filename, bc_file_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (bc_file_set == false) {
    PetscPrintf(grid.com,
                "PISM ERROR: Command-line option %s is required.\n", option.c_str());
    PISMEnd();
  }

  ierr = verbPrintf(2, grid.com,
                    "  - Reading air temperature and precipitation from '%s'...\n",
                    filename.c_str()); CHKERRQ(ierr);

  PIO nc(grid.com, "netcdf3", grid.get_unit_system());
  ierr = nc.open(filename, PISM_READONLY); CHKERRQ(ierr);
  {
    ierr = m_precipitation.read(nc, grid.time); CHKERRQ(ierr);
    ierr = m_air_temperature.read(nc, grid.time); CHKERRQ(ierr);
  }
  ierr = nc.close(); CHKERRQ(ierr);


  return 0;
}

PetscErrorCode PAWeatherStation::update(double t, double dt) {
  m_t = t;
  m_dt = dt;
  return 0;
}

PetscErrorCode PAWeatherStation::mean_precipitation(IceModelVec2S &result) {
  const double one_week = 7 * 24 * 60 * 60;

  unsigned int N = (unsigned int)(ceil(m_dt / one_week)); // one point per week

  PetscErrorCode ierr = result.set(m_precipitation.average(m_t, m_dt, N)); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAWeatherStation::mean_annual_temp(IceModelVec2S &result) {
  const double one_week = 7 * 24 * 60 * 60;

  unsigned int N = (unsigned int)(ceil(m_dt / one_week)); // one point per week

  PetscErrorCode ierr = result.set(m_air_temperature.average(m_t, m_dt, N)); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAWeatherStation::begin_pointwise_access() {
  // empty
  return 0;
}

PetscErrorCode PAWeatherStation::end_pointwise_access() {
  // empty
  return 0;
}

PetscErrorCode PAWeatherStation::init_timeseries(const std::vector<double> &ts) {
  size_t N = ts.size();

  m_precip_values.resize(N);
  m_air_temp_values.resize(N);

  for (unsigned int k = 0; k < N; ++k) {
    m_precip_values[k]   = m_precipitation(ts[k]);
    m_air_temp_values[k] = m_air_temperature(ts[k]);
  }

  return 0;
}

PetscErrorCode PAWeatherStation::precip_time_series(int i, int j,
                                                    std::vector<double> &result) {
  (void)i;
  (void)j;

  result = m_precip_values;
  return 0;
}

PetscErrorCode PAWeatherStation::temp_time_series(int i, int j,
                                                  std::vector<double> &result) {
  (void)i;
  (void)j;

  result = m_air_temp_values;
  return 0;
}

PetscErrorCode PAWeatherStation::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr = result.set(m_air_temperature(m_t + 0.5*m_dt)); CHKERRQ(ierr);
  return 0;
}

void PAWeatherStation::add_vars_to_output(const std::string &keyword,
                                          std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big") {
    result.insert("air_temp");
    result.insert("precipitation");
  }
}

PetscErrorCode PAWeatherStation::define_variables(const std::set<std::string> &vars,
                                                  const PIO &nc, IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp")) {
    // don't write using glaciological units
    ierr = m_air_temp_metadata.define(nc, nctype, false); CHKERRQ(ierr);
  }

  if (set_contains(vars, "precipitation")) {
    // do write using glaciological units
    ierr = m_precip_metadata.define(nc, nctype, true); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAWeatherStation::write_variables(const std::set<std::string> &vars,
                                                 const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "air_temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
    tmp.metadata() = m_air_temp_metadata;

    ierr = mean_annual_temp(tmp); CHKERRQ(ierr);

    ierr = tmp.write(nc); CHKERRQ(ierr);
  }

  if (set_contains(vars, "precipitation")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "precipitation", WITHOUT_GHOSTS); CHKERRQ(ierr);
    tmp.metadata() = m_precip_metadata;

    ierr = mean_precipitation(tmp); CHKERRQ(ierr);

    tmp.write_in_glaciological_units = true;
    ierr = tmp.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

} // end of namespace pism

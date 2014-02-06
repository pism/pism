// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#include "PAAnomaly.hh"
#include "PISMConfig.hh"
#include "IceGrid.hh"
#include <assert.h>

PAAnomaly::PAAnomaly(IceGrid &g, const PISMConfig &conf, PISMAtmosphereModel* in)
  : PGivenClimate<PAModifier,PISMAtmosphereModel>(g, conf, in),
    air_temp(g.get_unit_system()),
    precipitation(g.get_unit_system())
{
  PetscErrorCode ierr = allocate_PAAnomaly(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

PetscErrorCode PAAnomaly::allocate_PAAnomaly() {
  PetscErrorCode ierr;

  option_prefix	 = "-atmosphere_anomaly";

  // will be de-allocated by the parent's destructor
  air_temp_anomaly      = new IceModelVec2T;
  precipitation_anomaly = new IceModelVec2T;

  m_fields["air_temp_anomaly"]      = air_temp_anomaly;
  m_fields["precipitation_anomaly"] = precipitation_anomaly;

  ierr = process_options(); CHKERRQ(ierr);

  std::map<std::string, std::string> standard_names;
  ierr = set_vec_parameters(standard_names); CHKERRQ(ierr);

  ierr = air_temp_anomaly->create(grid, "air_temp_anomaly", false); CHKERRQ(ierr);
  ierr = air_temp_anomaly->set_attrs("climate_forcing",
                                     "anomaly of the near-surface air temperature",
                                     "Kelvin", ""); CHKERRQ(ierr);

  ierr = precipitation_anomaly->create(grid, "precipitation_anomaly", false); CHKERRQ(ierr);
  ierr = precipitation_anomaly->set_attrs("climate_forcing",
                                          "anomaly of the ice-equivalent precipitation rate",
                                          "m s-1", ""); CHKERRQ(ierr);
  ierr = precipitation_anomaly->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  precipitation_anomaly->write_in_glaciological_units = true;

  air_temp.init_2d("air_temp", grid);
  air_temp.set_string("pism_intent", "diagnostic");
  air_temp.set_string("long_name", "near-surface air temperature");
  ierr = air_temp.set_units("K"); CHKERRQ(ierr);

  precipitation.init_2d("precipitation", grid);
  precipitation.set_string("pism_intent", "diagnostic");
  precipitation.set_string("long_name", "precipitation, units of ice-equivalent thickness per time");
  ierr = precipitation.set_units("m / s"); CHKERRQ(ierr);
  ierr = precipitation.set_glaciological_units("m / year"); CHKERRQ(ierr);

  return 0;
}

PAAnomaly::~PAAnomaly()
{
  // empty
}

PetscErrorCode PAAnomaly::init(PISMVars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  assert(input_model != NULL);
  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the -atmosphere ...,anomaly code...\n"); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "    reading anomalies from %s ...\n",
                    filename.c_str()); CHKERRQ(ierr);

  ierr = air_temp_anomaly->init(filename, bc_period, bc_reference_time); CHKERRQ(ierr);
  ierr = precipitation_anomaly->init(filename, bc_period, bc_reference_time); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAAnomaly::update(double my_t, double my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = precipitation_anomaly->average(m_t, m_dt); CHKERRQ(ierr);
  ierr = air_temp_anomaly->average(m_t, m_dt); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PAAnomaly::mean_precipitation(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->mean_precipitation(result); CHKERRQ(ierr);

  return result.add(1.0, *precipitation_anomaly);
}

PetscErrorCode PAAnomaly::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->mean_annual_temp(result); CHKERRQ(ierr);

  return result.add(1.0, *air_temp_anomaly);
}

PetscErrorCode PAAnomaly::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->temp_snapshot(result); CHKERRQ(ierr);

  return result.add(1.0, *air_temp_anomaly);
}


PetscErrorCode PAAnomaly::begin_pointwise_access() {
  PetscErrorCode ierr = input_model->begin_pointwise_access(); CHKERRQ(ierr);
  ierr = air_temp_anomaly->begin_access(); CHKERRQ(ierr);
  ierr = precipitation_anomaly->begin_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAAnomaly::end_pointwise_access() {
  PetscErrorCode ierr = input_model->end_pointwise_access(); CHKERRQ(ierr);
  ierr = precipitation_anomaly->end_access(); CHKERRQ(ierr);
  ierr = air_temp_anomaly->end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAAnomaly::init_timeseries(double *ts, unsigned int N) {
  PetscErrorCode ierr;
  ierr = input_model->init_timeseries(ts, N); CHKERRQ(ierr);

  ierr = air_temp_anomaly->init_interpolation(ts, N); CHKERRQ(ierr);

  ierr = precipitation_anomaly->init_interpolation(ts, N); CHKERRQ(ierr);

  m_ts_times.resize(N);
  for (unsigned int k = 0; k < m_ts_times.size(); ++k)
    m_ts_times[k] = ts[k];
  
  return 0;
}

PetscErrorCode PAAnomaly::temp_time_series(int i, int j, double *result) {
  PetscErrorCode ierr;

  ierr = input_model->temp_time_series(i, j, result); CHKERRQ(ierr);

  m_temp_anomaly.reserve(m_ts_times.size());
  ierr = air_temp_anomaly->interp(i, j, &m_temp_anomaly[0]); CHKERRQ(ierr);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k)
    result[k] += m_temp_anomaly[k];

  return 0;
}

PetscErrorCode PAAnomaly::precip_time_series(int i, int j, double *result) {
  PetscErrorCode ierr;

  ierr = input_model->precip_time_series(i, j, result); CHKERRQ(ierr);

  m_mass_flux_anomaly.reserve(m_ts_times.size());
  ierr = precipitation_anomaly->interp(i, j, &m_mass_flux_anomaly[0]); CHKERRQ(ierr);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k)
    result[k] += m_mass_flux_anomaly[k];

  return 0;
}

void PAAnomaly::add_vars_to_output(std::string keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("air_temp");
    result.insert("precipitation");
  }
}


PetscErrorCode PAAnomaly::define_variables(std::set<std::string> vars, const PIO &nc,
                                           PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp")) {
    ierr = air_temp.define(nc, nctype, false); CHKERRQ(ierr);
    vars.erase("air_temp");
  }

  if (set_contains(vars, "precipitation")) {
    ierr = precipitation.define(nc, nctype, true); CHKERRQ(ierr);
    vars.erase("precipitation");
  }

  ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PAAnomaly::write_variables(std::set<std::string> vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "air_temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
    tmp.metadata() = air_temp;

    ierr = mean_annual_temp(tmp); CHKERRQ(ierr);

    ierr = tmp.write(nc); CHKERRQ(ierr);

    vars.erase("air_temp");
  }

  if (set_contains(vars, "precipitation")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "precipitation", WITHOUT_GHOSTS); CHKERRQ(ierr);
    tmp.metadata() = precipitation;

    ierr = mean_precipitation(tmp); CHKERRQ(ierr);

    tmp.write_in_glaciological_units = true;
    ierr = tmp.write(nc); CHKERRQ(ierr);

    vars.erase("precipitation");
  }

  ierr = input_model->write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}


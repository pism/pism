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

namespace pism {

PAAnomaly::PAAnomaly(IceGrid &g, AtmosphereModel* in)
  : PGivenClimate<PAModifier,AtmosphereModel>(g, in),
    air_temp(g.config.get_unit_system(), "air_temp", m_grid),
    precipitation(g.config.get_unit_system(), "precipitation", m_grid)
{
  option_prefix  = "-atmosphere_anomaly";

  // will be de-allocated by the parent's destructor
  air_temp_anomaly      = new IceModelVec2T;
  precipitation_anomaly = new IceModelVec2T;

  m_fields["air_temp_anomaly"]      = air_temp_anomaly;
  m_fields["precipitation_anomaly"] = precipitation_anomaly;

  process_options();

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  air_temp_anomaly->create(m_grid, "air_temp_anomaly", false);
  air_temp_anomaly->set_attrs("climate_forcing",
                              "anomaly of the near-surface air temperature",
                              "Kelvin", "");

  precipitation_anomaly->create(m_grid, "precipitation_anomaly", false);
  precipitation_anomaly->set_attrs("climate_forcing",
                                   "anomaly of the ice-equivalent precipitation rate",
                                   "m s-1", "");
  precipitation_anomaly->set_glaciological_units("m year-1");
  precipitation_anomaly->write_in_glaciological_units = true;

  air_temp.set_string("pism_intent", "diagnostic");
  air_temp.set_string("long_name", "near-surface air temperature");
  air_temp.set_units("K");

  precipitation.set_string("pism_intent", "diagnostic");
  precipitation.set_string("long_name", "precipitation, units of ice-equivalent thickness per time");
  precipitation.set_units("m / s");
  precipitation.set_glaciological_units("m / year");
}

PAAnomaly::~PAAnomaly()
{
  // empty
}

void PAAnomaly::init() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  assert(input_model != NULL);
  input_model->init();

  verbPrintf(2, m_grid.com,
             "* Initializing the -atmosphere ...,anomaly code...\n");

  verbPrintf(2, m_grid.com,
             "    reading anomalies from %s ...\n",
             filename.c_str());

  air_temp_anomaly->init(filename, bc_period, bc_reference_time);
  precipitation_anomaly->init(filename, bc_period, bc_reference_time);
}

void PAAnomaly::update(double my_t, double my_dt) {
  update_internal(my_t, my_dt);

  precipitation_anomaly->average(m_t, m_dt);
  air_temp_anomaly->average(m_t, m_dt);
}


void PAAnomaly::mean_precipitation(IceModelVec2S &result) {
  input_model->mean_precipitation(result);

  result.add(1.0, *precipitation_anomaly);
}

void PAAnomaly::mean_annual_temp(IceModelVec2S &result) {
  input_model->mean_annual_temp(result);

  result.add(1.0, *air_temp_anomaly);
}

void PAAnomaly::temp_snapshot(IceModelVec2S &result) {
  input_model->temp_snapshot(result);

  result.add(1.0, *air_temp_anomaly);
}


void PAAnomaly::begin_pointwise_access() {
  input_model->begin_pointwise_access();
  air_temp_anomaly->begin_access();
  precipitation_anomaly->begin_access();
}

void PAAnomaly::end_pointwise_access() {
  input_model->end_pointwise_access();
  precipitation_anomaly->end_access();
  air_temp_anomaly->end_access();
}

void PAAnomaly::init_timeseries(const std::vector<double> &ts) {
  input_model->init_timeseries(ts);

  air_temp_anomaly->init_interpolation(ts);

  precipitation_anomaly->init_interpolation(ts);

  m_ts_times = ts;
}

void PAAnomaly::temp_time_series(int i, int j, std::vector<double> &result) {
  input_model->temp_time_series(i, j, result);

  m_temp_anomaly.reserve(m_ts_times.size());
  air_temp_anomaly->interp(i, j, m_temp_anomaly);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] += m_temp_anomaly[k];
  }
}

void PAAnomaly::precip_time_series(int i, int j, std::vector<double> &result) {
  input_model->precip_time_series(i, j, result);

  m_mass_flux_anomaly.reserve(m_ts_times.size());
  precipitation_anomaly->interp(i, j, m_mass_flux_anomaly);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] += m_mass_flux_anomaly[k];
  }
}

void PAAnomaly::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("air_temp");
    result.insert("precipitation");
  }
}

void PAAnomaly::define_variables(const std::set<std::string> &vars_input, const PIO &nc,
                                           IO_Type nctype) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "air_temp")) {
    air_temp.define(nc, nctype, false);
    vars.erase("air_temp");
  }

  if (set_contains(vars, "precipitation")) {
    precipitation.define(nc, nctype, true);
    vars.erase("precipitation");
  }

  input_model->define_variables(vars, nc, nctype);
}


void PAAnomaly::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "air_temp")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "air_temp", WITHOUT_GHOSTS);
    tmp.metadata() = air_temp;

    mean_annual_temp(tmp);

    tmp.write(nc);

    vars.erase("air_temp");
  }

  if (set_contains(vars, "precipitation")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "precipitation", WITHOUT_GHOSTS);
    tmp.metadata() = precipitation;

    mean_precipitation(tmp);

    tmp.write_in_glaciological_units = true;
    tmp.write(nc);

    vars.erase("precipitation");
  }

  input_model->write_variables(vars, nc);
}


} // end of namespace pism

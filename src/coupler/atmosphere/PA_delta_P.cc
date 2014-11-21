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

#include "PA_delta_P.hh"
#include "PISMConfig.hh"

namespace pism {

PA_delta_P::PA_delta_P(IceGrid &g, const Config &conf, AtmosphereModel* in)
  : PScalarForcing<AtmosphereModel,PAModifier>(g, conf, in),
    air_temp(g.get_unit_system(), "air_temp", grid),
    precipitation(g.get_unit_system(), "precipitation", grid)
{
  offset = NULL;
  PetscErrorCode ierr = allocate_PA_delta_P(); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    throw std::runtime_error("PA_delta_P allocation failed");
  }
}

PetscErrorCode PA_delta_P::allocate_PA_delta_P() {
  PetscErrorCode ierr;

  option_prefix = "-atmosphere_delta_P";
  offset_name = "delta_P";
  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));
  offset->get_metadata().set_units("m / second");
  offset->get_metadata().set_glaciological_units("m / year");
  offset->get_metadata().set_string("long_name",
                                    "precipitation offsets, units of ice-equivalent thickness");
  offset->get_dimension_metadata().set_units(grid.time->units_string());


  air_temp.set_string("pism_intent", "diagnostic");
  air_temp.set_string("long_name", "near-surface air temperature");
  ierr = air_temp.set_units("K"); CHKERRQ(ierr);

  precipitation.set_string("pism_intent", "diagnostic");
  precipitation.set_string("long_name", "precipitation, units of ice-equivalent thickness per time");
  ierr = precipitation.set_units("m / s"); CHKERRQ(ierr);
  ierr = precipitation.set_glaciological_units("m / year"); CHKERRQ(ierr);

  return 0;
}

PA_delta_P::~PA_delta_P()
{
  // empty
}

void PA_delta_P::init(Vars &vars) {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  input_model->init(vars);

  verbPrintf(2, grid.com,
             "* Initializing precipitation forcing using scalar offsets...\n");

  init_internal();
}

void PA_delta_P::init_timeseries(const std::vector<double> &ts) {
  PAModifier::init_timeseries(ts);

  m_offset_values.resize(m_ts_times.size());
  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    m_offset_values[k] = (*offset)(m_ts_times[k]);
  }
}



void PA_delta_P::mean_precipitation(IceModelVec2S &result) {
  input_model->mean_precipitation(result);
  offset_data(result);
}

void PA_delta_P::precip_time_series(int i, int j, std::vector<double> &result) {
  input_model->precip_time_series(i, j, result);
  
  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] += m_offset_values[k];
  }
}

void PA_delta_P::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("air_temp");
    result.insert("precipitation");
  }
}


void PA_delta_P::define_variables(const std::set<std::string> &vars_input, const PIO &nc,
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


void PA_delta_P::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "air_temp")) {
    IceModelVec2S tmp;
    tmp.create(grid, "air_temp", WITHOUT_GHOSTS);
    tmp.metadata() = air_temp;

    mean_annual_temp(tmp);

    tmp.write(nc);

    vars.erase("air_temp");
  }

  if (set_contains(vars, "precipitation")) {
    IceModelVec2S tmp;
    tmp.create(grid, "precipitation", WITHOUT_GHOSTS);
    tmp.metadata() = precipitation;

    mean_precipitation(tmp);

    tmp.write_in_glaciological_units = true;
    tmp.write(nc);

    vars.erase("precipitation");
  }

  input_model->write_variables(vars, nc);
}

} // end of namespace pism

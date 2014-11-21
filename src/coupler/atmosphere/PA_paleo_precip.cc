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

#include "PA_paleo_precip.hh"
#include "PISMConfig.hh"

namespace pism {

PA_paleo_precip::PA_paleo_precip(IceGrid &g, const Config &conf, AtmosphereModel* in)
  : PScalarForcing<AtmosphereModel,PAModifier>(g, conf, in),
    air_temp(g.get_unit_system(), "air_temp", grid),
    precipitation(g.get_unit_system(), "precipitation", grid)
{
  offset = NULL;
  PetscErrorCode ierr = allocate_PA_paleo_precip(); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    throw std::runtime_error("PA_paleo_precip allocation failed");
  }
}

PetscErrorCode PA_paleo_precip::allocate_PA_paleo_precip() {
  PetscErrorCode ierr;

  option_prefix = "-atmosphere_paleo_precip";
  offset_name = "delta_T";
  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));
  offset->get_metadata().set_units("Kelvin");
  offset->get_metadata().set_string("long_name", "air temperature offsets");
  offset->get_dimension_metadata().set_units(grid.time->units_string());

  air_temp.set_string("pism_intent", "diagnostic");
  air_temp.set_string("long_name", "near-surface air temperature");
  ierr = air_temp.set_units("K"); CHKERRQ(ierr);

  precipitation.set_string("pism_intent", "diagnostic");
  precipitation.set_string("long_name", "precipitation, units of ice-equivalent thickness per time");
  ierr = precipitation.set_units("m / s"); CHKERRQ(ierr);
  ierr = precipitation.set_glaciological_units("m / year"); CHKERRQ(ierr);

  m_precipexpfactor = config.get("precip_exponential_factor_for_temperature");

  return 0;
}

PA_paleo_precip::~PA_paleo_precip()
{
  // empty
}

void PA_paleo_precip::init(Vars &vars) {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  input_model->init(vars);

  verbPrintf(2, grid.com,
             "* Initializing paleo-precipitation correction using temperature offsets...\n");

  init_internal();
}

void PA_paleo_precip::init_timeseries(const std::vector<double> &ts) {

  PAModifier::init_timeseries(ts);

  size_t N = ts.size();

  m_scaling_values.resize(N);
  for (unsigned int k = 0; k < N; ++k) {
    m_scaling_values[k] = exp(m_precipexpfactor * (*offset)(m_ts_times[k]));
  }
}

void PA_paleo_precip::mean_precipitation(IceModelVec2S &result) {
  input_model->mean_precipitation(result);
  result.scale(exp(m_precipexpfactor * (*offset)(m_t + 0.5 * m_dt)));
}

void PA_paleo_precip::precip_time_series(int i, int j, std::vector<double> &result) {
  input_model->precip_time_series(i, j, result);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] *= m_scaling_values[k];
  }
}

void PA_paleo_precip::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("air_temp");
    result.insert("precipitation");
  }
}


void PA_paleo_precip::define_variables(const std::set<std::string> &vars_input, const PIO &nc,
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


void PA_paleo_precip::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
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

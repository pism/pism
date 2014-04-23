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

PA_paleo_precip::PA_paleo_precip(IceGrid &g, const PISMConfig &conf, PISMAtmosphereModel* in)
  : PScalarForcing<PISMAtmosphereModel,PAModifier>(g, conf, in),
    air_temp(g.get_unit_system()),
    precipitation(g.get_unit_system())
{
  offset = NULL;
  PetscErrorCode ierr = allocate_PA_paleo_precip(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

PetscErrorCode PA_paleo_precip::allocate_PA_paleo_precip() {
  PetscErrorCode ierr;

  option_prefix = "-atmosphere_paleo_precip";
  offset_name = "delta_T";
  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));
  offset->set_units("Kelvin", "Kelvin");
  offset->set_dimension_units(grid.time->units_string(), "");
  offset->set_attr("long_name", "air temperature offsets");

  air_temp.init_2d("air_temp", grid);
  air_temp.set_string("pism_intent", "diagnostic");
  air_temp.set_string("long_name", "near-surface air temperature");
  ierr = air_temp.set_units("K"); CHKERRQ(ierr);

  precipitation.init_2d("precipitation", grid);
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

PetscErrorCode PA_paleo_precip::init(PISMVars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing paleo-precipitation correction using temperature offsets...\n"); CHKERRQ(ierr);

  ierr = init_internal(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PA_paleo_precip::init_timeseries(double *ts, unsigned int N) {
  PetscErrorCode ierr;

  ierr = PAModifier::init_timeseries(ts, N); CHKERRQ(ierr);

  m_scaling_values.resize(N);
  for (unsigned int k = 0; k < N; ++k)
    m_scaling_values[k] = exp( m_precipexpfactor * (*offset)(m_ts_times[k]));

  return 0;
}

PetscErrorCode PA_paleo_precip::mean_precipitation(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->mean_precipitation(result);
  ierr = result.scale(exp( m_precipexpfactor * (*offset)(m_t + 0.5 * m_dt) )); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PA_paleo_precip::precip_time_series(int i, int j, double *result) {
  PetscErrorCode ierr = input_model->precip_time_series(i, j, result); CHKERRQ(ierr);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k)
    result[k] *= m_scaling_values[k];

  return 0;
}

void PA_paleo_precip::add_vars_to_output(std::string keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("air_temp");
    result.insert("precipitation");
  }
}


PetscErrorCode PA_paleo_precip::define_variables(std::set<std::string> vars, const PIO &nc,
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


PetscErrorCode PA_paleo_precip::write_variables(std::set<std::string> vars, const PIO &nc) {
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

} // end of namespace pism

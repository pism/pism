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

#include "PA_frac_P.hh"
#include "PISMConfig.hh"

namespace pism {

PA_frac_P::PA_frac_P(IceGrid &g, const Config &conf, AtmosphereModel* in)
  : PScalarForcing<AtmosphereModel,PAModifier>(g, conf, in),
    air_temp(g.get_unit_system()),
    precipitation(g.get_unit_system())
{
  offset = NULL;
  PetscErrorCode ierr = allocate_PA_frac_P(); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    throw std::runtime_error("PA_frac_P allocation failed");
  }
}

PetscErrorCode PA_frac_P::allocate_PA_frac_P() {
  PetscErrorCode ierr;

  option_prefix = "-atmosphere_frac_P";
  offset_name = "frac_P";
  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));
  offset->get_metadata().set_units("1");
  offset->get_metadata().set_string("long_name", "precipitation multiplier, pure fraction");
  offset->get_dimension_metadata().set_units(grid.time->units_string());

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

PA_frac_P::~PA_frac_P()
{
  // empty; "offset" is deleted by ~PScalarForcing().
}

PetscErrorCode PA_frac_P::init(Vars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing precipitation forcing using scalar multipliers...\n"); CHKERRQ(ierr);

  ierr = init_internal(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PA_frac_P::init_timeseries(const std::vector<double> &ts) {
  PetscErrorCode ierr;

  ierr = PAModifier::init_timeseries(ts); CHKERRQ(ierr);

  m_offset_values.resize(m_ts_times.size());
  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    m_offset_values[k] = (*offset)(m_ts_times[k]);
  }

  return 0;
}

PetscErrorCode PA_frac_P::mean_precipitation(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->mean_precipitation(result);
  ierr = scale_data(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PA_frac_P::precip_time_series(int i, int j, std::vector<double> &result) {
  PetscErrorCode ierr = input_model->precip_time_series(i, j, result); CHKERRQ(ierr);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] *= m_offset_values[k];
  }
  return 0;
}

void PA_frac_P::add_vars_to_output(const std::string &keyword,
                                   std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("air_temp");
    result.insert("precipitation");
  }
}

PetscErrorCode PA_frac_P::define_variables(const std::set<std::string> &vars_input,
                                           const PIO &nc, IO_Type nctype) {
  std::set<std::string> vars = vars_input;
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

PetscErrorCode PA_frac_P::write_variables(const std::set<std::string> &vars_input,
                                          const PIO &nc) {
  std::set<std::string> vars = vars_input;
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

// Copyright (C) 2011, 2012 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include "PA_delta_T.hh"

/// delta_T forcing of near-surface air temperatures

PA_delta_T::PA_delta_T(IceGrid &g, const NCConfigVariable &conf, PISMAtmosphereModel* in)
  : PScalarForcing<PISMAtmosphereModel,PAModifier>(g, conf, in)
{
  option_prefix = "-atmosphere_delta_T";
  offset_name = "delta_T";
  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));
  offset->set_units("Kelvin", "");
  offset->set_dimension_units(grid.time->units(), "");
  offset->set_attr("long_name", "near-surface air temperature offsets");
}

PetscErrorCode PA_delta_T::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing near-surface air temperature forcing using scalar offsets...\n"); CHKERRQ(ierr);

  ierr = init_internal(); CHKERRQ(ierr);

  air_temp.init_2d("air_temp", grid);
  air_temp.set_string("pism_intent", "diagnostic");
  air_temp.set_string("long_name", "near-surface air temperature");
  ierr = air_temp.set_units("K"); CHKERRQ(ierr);

  precipitation.init_2d("precipitation", grid);
  precipitation.set_string("pism_intent", "diagnostic");
  precipitation.set_string("long_name", "near-surface air temperature");
  ierr = precipitation.set_units("m / s"); CHKERRQ(ierr);
  ierr = precipitation.set_glaciological_units("m / year"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PA_delta_T::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->mean_annual_temp(result); CHKERRQ(ierr);
  ierr = offset_data(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PA_delta_T::temp_time_series(int i, int j, int N,
                                             PetscReal *ts, PetscReal *values) {
  PetscErrorCode ierr = input_model->temp_time_series(i, j, N, ts, values); CHKERRQ(ierr);
  
  if (offset) {
    for (int k = 0; k < N; ++k)
      values[k] += (*offset)(ts[k]);
  }

  return 0;
}

PetscErrorCode PA_delta_T::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->temp_snapshot(result); CHKERRQ(ierr);
  ierr = offset_data(result); CHKERRQ(ierr);
  return 0;
}

void PA_delta_T::add_vars_to_output(string keyword,
                                    map<string,NCSpatialVariable> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result["air_temp"] = air_temp;
    result["precipitation"] = precipitation;
  }
}


PetscErrorCode PA_delta_T::define_variables(set<string> vars, const PIO &nc,
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


PetscErrorCode PA_delta_T::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "air_temp", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(air_temp, 0); CHKERRQ(ierr);

    ierr = mean_annual_temp(tmp); CHKERRQ(ierr);

    ierr = tmp.write(filename); CHKERRQ(ierr);

    vars.erase("air_temp");
  }

  if (set_contains(vars, "precipitation")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "precipitation", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(air_temp, 0); CHKERRQ(ierr);

    ierr = mean_precipitation(tmp); CHKERRQ(ierr);

    tmp.write_in_glaciological_units = true;
    ierr = tmp.write(filename); CHKERRQ(ierr);

    vars.erase("precipitation");
  }

  ierr = input_model->write_variables(vars, filename); CHKERRQ(ierr);

  return 0;
}

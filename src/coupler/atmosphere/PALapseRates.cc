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

#include "PALapseRates.hh"

PetscErrorCode PALapseRates::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool precip_lapse_rate_set;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "  [using air temperature and precipitation lapse corrections]\n"); CHKERRQ(ierr);

  ierr = init_internal(vars); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Lapse rate options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsReal("-precip_lapse_rate",
                           "Elevation lapse rate for the surface mass balance, in m/year per km",
                           precip_lapse_rate, precip_lapse_rate_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "   air temperature lapse rate: %3.3f K per km\n"
                    "   precipitation lapse rate:   %3.3f m/year per km\n",
                    temp_lapse_rate, precip_lapse_rate); CHKERRQ(ierr);

  temp_lapse_rate = convert(temp_lapse_rate, "K/km", "K/m");

  precip_lapse_rate = convert(precip_lapse_rate, "m/year / km", "m/s / m");

  precipitation.init_2d("precipitation", grid);
  precipitation.set_string("pism_intent", "diagnostic");
  precipitation.set_string("long_name",
                    "ice-equivalent precipitation rate with a lapse-rate correction");
  ierr = precipitation.set_units("m s-1"); CHKERRQ(ierr);
  ierr = precipitation.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  air_temp.init_2d("air_temp", grid);
  air_temp.set_string("pism_intent", "diagnostic");
  air_temp.set_string("long_name",
                      "near-surface air temperature with a lapse-rate correction");
  ierr = air_temp.set_units("K"); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PALapseRates::mean_precipitation(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = input_model->mean_precipitation(result); CHKERRQ(ierr);
  ierr = lapse_rate_correction(result, precip_lapse_rate); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PALapseRates::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = input_model->mean_annual_temp(result); CHKERRQ(ierr);
  ierr = lapse_rate_correction(result, temp_lapse_rate); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PALapseRates::begin_pointwise_access() {
  PetscErrorCode ierr;
  ierr = input_model->begin_pointwise_access(); CHKERRQ(ierr);
  ierr = reference_surface.begin_access(); CHKERRQ(ierr);
  ierr = surface->begin_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PALapseRates::end_pointwise_access() {
  PetscErrorCode ierr;
  ierr = input_model->end_pointwise_access(); CHKERRQ(ierr);
  ierr = reference_surface.end_access(); CHKERRQ(ierr);
  ierr = surface->end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PALapseRates::temp_time_series(int i, int j, int N,
                                              PetscReal *ts, PetscReal *values) {
  PetscErrorCode ierr;
  vector<PetscScalar> usurf(N);

  ierr = input_model->temp_time_series(i, j, N, ts, values); CHKERRQ(ierr);

  ierr = reference_surface.interp(i, j, N, ts, &usurf[0]); CHKERRQ(ierr);

  for (int m = 0; m < N; ++m) {
    values[m] -= temp_lapse_rate * ((*surface)(i, j) - usurf[m]);
  }

  return 0;
}

PetscErrorCode PALapseRates::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = input_model->temp_snapshot(result); CHKERRQ(ierr);
  ierr = lapse_rate_correction(result, temp_lapse_rate); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PALapseRates::define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp")) {
    ierr = air_temp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "precipitation")) {
    ierr = precipitation.define(nc, nctype, true); CHKERRQ(ierr);
  }

  ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PALapseRates::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "air_temp", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(air_temp, 0); CHKERRQ(ierr);

    ierr = temp_snapshot(tmp); CHKERRQ(ierr);

    ierr = tmp.write(nc); CHKERRQ(ierr);

    vars.erase("air_temp");
  }

  if (set_contains(vars, "precipitation")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "precipitation", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(precipitation, 0); CHKERRQ(ierr);

    ierr = mean_precipitation(tmp); CHKERRQ(ierr);
    tmp.write_in_glaciological_units = true;
    ierr = tmp.write(nc); CHKERRQ(ierr);

    vars.erase("precipitation");
  }

  ierr = input_model->write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}

void PALapseRates::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result["air_temp"] = air_temp;
    result["precipitation"]   = precipitation;
  }
}

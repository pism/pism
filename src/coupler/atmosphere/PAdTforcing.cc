// Copyright (C) 2011 PISM Authors
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

#include "PAdTforcing.hh"

/// -dTforcing of near-surface air temperatures

PAdTforcing::PAdTforcing(IceGrid &g, const NCConfigVariable &conf, PISMAtmosphereModel* in)
  : PScalarForcing<PISMAtmosphereModel,PAModifier>(g, conf, in)
{
  option = "-dTforcing";
  offset_name = "delta_T";
  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));
  offset->set_units("Celsius", "");
  offset->set_dimension_units("seconds", "");
  offset->set_attr("long_name", "near-surface air temperature offsets");
}

PetscErrorCode PAdTforcing::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing near-surface air temperature forcing using scalar offsets...\n"); CHKERRQ(ierr);

  ierr = init_internal(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAdTforcing::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->mean_annual_temp(result); CHKERRQ(ierr);
  ierr = offset_data(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAdTforcing::temp_time_series(int i, int j, int N,
                                             PetscReal *ts, PetscReal *values) {
  PetscErrorCode ierr = input_model->temp_time_series(i, j, N, ts, values); CHKERRQ(ierr);
  
  if (offset) {
    for (int k = 0; k < N; ++k)
      values[k] += (*offset)(ts[k]);
  }

  return 0;
}

PetscErrorCode PAdTforcing::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->temp_snapshot(result); CHKERRQ(ierr);
  ierr = offset_data(result); CHKERRQ(ierr);
  return 0;
}

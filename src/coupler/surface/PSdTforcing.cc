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

#include "PSdTforcing.hh"

/// -dTforcing of ice surface temperatures

PSdTforcing::PSdTforcing(IceGrid &g, const NCConfigVariable &conf, PISMSurfaceModel* in)
  : PScalarForcing<PISMSurfaceModel,PSModifier>(g, conf, in)
{
  option = "-dTforcing";
  offset_name = "delta_T";
  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));
  offset->set_units("Celsius", "");
  offset->set_dimension_units(grid.time->units(), "");
  offset->set_attr("long_name", "ice-surface temperature offsets");
}

PetscErrorCode PSdTforcing::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing ice-surface temperature forcing using scalar offsets...\n"); CHKERRQ(ierr);

  ierr = init_internal(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSdTforcing::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = input_model->ice_surface_temperature(result); CHKERRQ(ierr);
  ierr = offset_data(result); CHKERRQ(ierr);
  return 0;
}

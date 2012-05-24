// Copyright (C) 2012 PISM Authors
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

#include "PACosineYearlyCycle.hh"
#include "pism_options.hh"

PetscErrorCode PACosineYearlyCycle::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool flag;
  string input_file;

  variables = &vars;
  snow_temp_july_day = config.get("snow_temp_july_day");

  ierr = verbPrintf(2, grid.com,
		    "* Initializing the 'cosine yearly cycle' atmosphere model (-atmosphere yearly_cycle)...\n");
  CHKERRQ(ierr);

  {
    // Allocate storage:
    ierr = air_temp_mean_annual.create(grid, "air_temp_mean_annual", false); CHKERRQ(ierr);
    ierr = air_temp_mean_annual.set_attrs("climate_state",
                                          "mean annual near-surface air temperature",
                                          "Kelvin",
                                          ""); CHKERRQ(ierr);
    air_temp_mean_annual.time_independent = true;

    ierr = air_temp_mean_july.create(grid, "air_temp_mean_july", false); CHKERRQ(ierr);
    ierr = air_temp_mean_july.set_attrs("climate_state",
                                        "mean July near-surface air temperature",
                                        "Kelvin",
                                        ""); CHKERRQ(ierr);
    air_temp_mean_july.time_independent = true;

    ierr = precipitation.create(grid, "precipitation", false); CHKERRQ(ierr);
    ierr = precipitation.set_attrs("climate_state",
                                   "mean annual ice-equivalent precipitation rate",
                                   "m s-1",
                                   ""); CHKERRQ(ierr);
    ierr = precipitation.set_glaciological_units("m year-1");
    precipitation.write_in_glaciological_units = true;
    precipitation.time_independent = true;

  }

  ierr = PetscOptionsBegin(grid.com, "", "Options controlling '-atmosphere yearly_cycle'",
                           ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-atmosphere_yearly_cycle_file", "PACosineYearlyCycle input file name",
                             input_file, flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (flag == false) {
    PetscPrintf(grid.com,
                "PISM ERROR: Please specify an '-atmosphere yearly_cycle' input file\n"
                "            using the -atmosphere_yearly_cycle_file option.\n");
    PISMEnd();
  }

  ierr = verbPrintf(2, grid.com,
                    "  Reading mean annual air temperature, mean July air temperature, and\n"
                    "  precipitation fields from '%s'...\n", input_file.c_str()); CHKERRQ(ierr);

  ierr = air_temp_mean_annual.regrid(input_file, true); CHKERRQ(ierr);
  ierr = air_temp_mean_july.regrid(input_file, true); CHKERRQ(ierr);
  ierr = precipitation.regrid(input_file, true); CHKERRQ(ierr);

  air_temp_snapshot.init_2d("air_temp_snapshot", grid);
  air_temp_snapshot.set_string("pism_intent", "diagnostic");
  air_temp_snapshot.set_string("long_name",
                               "snapshot of the near-surface air temperature");
  ierr = air_temp_snapshot.set_units("Kelvin"); CHKERRQ(ierr);


  return 0;
}


PetscErrorCode PACosineYearlyCycle::update(PetscReal my_t, PetscReal my_dt) {
  t = my_t;
  dt = my_dt;
  return 0;
}


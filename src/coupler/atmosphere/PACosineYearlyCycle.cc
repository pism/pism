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
#include "Timeseries.hh"
#include "PISMTime.hh"
#include "pism_options.hh"

PACosineYearlyCycle::~PACosineYearlyCycle() {
  if (A != NULL)
    delete A;
}

PetscErrorCode PACosineYearlyCycle::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool input_file_flag, scaling_flag;
  string input_file, scaling_file;

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
    ierr = PISMOptionsString("-atmosphere_yearly_cycle_file",
                             "PACosineYearlyCycle input file name",
                             input_file, input_file_flag); CHKERRQ(ierr);
    ierr = PISMOptionsString("-atmosphere_yearly_cycle_scaling_file",
                             "PACosineYearlyCycle amplitude scaling input file name",
                             scaling_file, scaling_flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (input_file_flag == false) {
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

  if (scaling_flag) {
    ierr = verbPrintf(2, grid.com,
                      "  Reading cosine yearly cycle amplitude scaling from '%s'...\n",
                      scaling_file.c_str()); CHKERRQ(ierr);

    A = new Timeseries(&grid, "amplitude_scaling",
                       config.get_string("time_dimension_name"));
    A->set_units("1", "1");
    A->set_dimension_units(grid.time->units(), "");
    A->set_attr("long_name", "cosine yearly cycle amplitude scaling");

    ierr = A->read(scaling_file, grid.time->use_reference_date()); CHKERRQ(ierr);
  } else {
    if (A != NULL)
      delete A;
    A = NULL;
  }

  return 0;
}


PetscErrorCode PACosineYearlyCycle::update(PetscReal my_t, PetscReal my_dt) {
  t = my_t;
  dt = my_dt;
  return 0;
}

PetscErrorCode PACosineYearlyCycle::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr;
  const PetscReal
    sperd = 8.64e4, // exact number of seconds per day
    julyday_fraction = (sperd / secpera) * snow_temp_july_day;

  double T = grid.time->year_fraction(t + 0.5 * dt) - julyday_fraction,
    scaling = 1;

  if (A != NULL) {
    scaling = (*A)(t + 0.5 * dt);
  }

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.begin_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_july.begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i,j) = air_temp_mean_annual(i,j) +
        scaling * (air_temp_mean_july(i,j) - air_temp_mean_annual(i,j)) * cos(2.0 * pi * T);
    }
  }

  ierr = air_temp_mean_july.end_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PACosineYearlyCycle::temp_time_series(int i, int j, int N,
                                                     PetscReal *ts, PetscReal *values) {
  // constants related to the standard yearly cycle
  const PetscReal
    sperd = 8.64e4, // exact number of seconds per day
    julyday_fraction = (sperd / secpera) * snow_temp_july_day;

  for (PetscInt k = 0; k < N; ++k) {
    double tk = grid.time->year_fraction(ts[k]) - julyday_fraction,
      scaling = 1;

    if (A != NULL) {
      scaling = (*A)(ts[k]);
    }

    values[k] = air_temp_mean_annual(i,j) +
      scaling * (air_temp_mean_july(i,j) - air_temp_mean_annual(i,j)) * cos(2.0 * pi * tk);
  }

  return 0;
}

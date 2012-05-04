// Copyright (C) 2007-2012 Ed Bueler and Nathan Shemonski and Constantine Khroulev
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

#include "PAEismintGreenland.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"

PA_EISMINT_Greenland::PA_EISMINT_Greenland(IceGrid &g, const NCConfigVariable &conf)
  : PAYearlyCycle(g, conf) {
  do_greenhouse_warming = false;
  greenhouse_warming_start_year = 0.0;
}

PetscErrorCode PA_EISMINT_Greenland::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  if ((fabs(my_t - t) < 1e-12) &&
      (fabs(my_dt - dt) < 1e-12))
    return 0;

  t  = my_t;
  dt = my_dt;

  ierr = surfelev->begin_access();   CHKERRQ(ierr);
  ierr = lat->begin_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.begin_access();  CHKERRQ(ierr);
  ierr = air_temp_mean_july.begin_access();  CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscReal Z = PetscMax((*surfelev)(i,j), 20 * ((*lat)(i,j) - 65));

      air_temp_mean_annual(i,j) = 49.13 - 0.007992 * Z - 0.7576 * (*lat)(i,j) + 273.15;
      air_temp_mean_july(i,j) = 273.15 + 30.38 - 0.006277 * (*surfelev)(i,j) - 0.3262 * (*lat)(i,j);
    }
  }
  ierr = surfelev->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.end_access();  CHKERRQ(ierr);
  ierr = air_temp_mean_july.end_access();  CHKERRQ(ierr);

  if (do_greenhouse_warming == PETSC_TRUE) {
    const PetscScalar shift = greenhouse_shift(my_t, my_dt);
    ierr = air_temp_mean_annual.shift(shift); CHKERRQ(ierr);
    ierr = air_temp_mean_july.shift(shift); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PA_EISMINT_Greenland::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
		    "* Initializing Greenland atmosphere model based on the EISMINT Greenland (C. Ritz, 1997)\n"
		    "  air temperature parameterization and using stored time-independent precipitation...\n"); CHKERRQ(ierr);

  reference = "Ritz, C. (1997). EISMINT Intercomparison Experiment: Comparison of existing Greenland models."
    " URL: http://homepages.vub.ac.be/~phuybrec/eismint/greenland.html";

  ierr = PAYearlyCycle::init(vars); CHKERRQ(ierr);

  // initialize pointers to fields the parameterization depends on:
  surfelev = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (!surfelev) SETERRQ(grid.com, 1, "ERROR: surface_altitude is not available");

  lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
  if (!lat) SETERRQ(grid.com, 1, "ERROR: latitude is not available");

  PetscBool gwl3_start_set;
  PetscReal  gwl3_start_year;
  ierr = PetscOptionsGetReal(PETSC_NULL, "-gwl3_start_year", &gwl3_start_year, &gwl3_start_set);
  CHKERRQ(ierr);

  if (gwl3_start_set) {  // do GWL3, and set start year
    do_greenhouse_warming = true;
    greenhouse_warming_start_year = gwl3_start_year;

    ierr = verbPrintf(2, grid.com,
		      "    turning on GWL3 warming scenario at year %3.3f...\n",
		      gwl3_start_year);
  }

  return 0;
}

PetscReal PA_EISMINT_Greenland::greenhouse_shift(PetscReal my_t, PetscReal my_dt) {
  // compute age back to start of GWL3; use midpoint of interval as the time
  PetscScalar age_years = grid.time->seconds_to_years(my_t + 0.5 * my_dt) - greenhouse_warming_start_year;

  if (age_years <= 0.0) {
    return 0.0; // before time 0.0, no warming
  } else {
    if (age_years <= 80.0) {
      return age_years * 0.035;
    } else if (age_years <= 500.0) {
      return 2.8 + (age_years - 80.0) * 0.0017;
    } else { // after 500 years, constant amount
      return 3.514;
    }
  }
}

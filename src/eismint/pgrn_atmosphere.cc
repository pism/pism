// Copyright (C) 2007-2011 Ed Bueler and Nathan Shemonski and Constantine Khroulev
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

# include "pgrn_atmosphere.hh"

PA_EISMINT_Greenland::PA_EISMINT_Greenland(IceGrid &g, const NCConfigVariable &conf)
  : PAYearlyCycle(g, conf) {
  do_greenhouse_warming = false;
  greenhouse_warming_start_year = 0.0;
}

PetscErrorCode PA_EISMINT_Greenland::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  if ((fabs(t_years - t) < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  t  = t_years;
  dt = dt_years;

  ierr = surfelev->begin_access();   CHKERRQ(ierr);
  ierr = lat->begin_access(); CHKERRQ(ierr);
  ierr = temp_ma.begin_access();  CHKERRQ(ierr);
  ierr = temp_mj.begin_access();  CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscReal Z = PetscMax((*surfelev)(i,j), 20 * ((*lat)(i,j) - 65));

      temp_ma(i,j) = 49.13 - 0.007992 * Z - 0.7576 * (*lat)(i,j) + 273.15;
      temp_mj(i,j) = 273.15 + 30.38 - 0.006277 * (*surfelev)(i,j) - 0.3262 * (*lat)(i,j);
    }
  }  
  ierr = surfelev->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = temp_ma.end_access();  CHKERRQ(ierr);
  ierr = temp_mj.end_access();  CHKERRQ(ierr);

  if (do_greenhouse_warming == PETSC_TRUE) {
    const PetscScalar shift = greenhouse_shift(t_years, dt_years);
    ierr = temp_ma.shift(shift); CHKERRQ(ierr);
    ierr = temp_mj.shift(shift); CHKERRQ(ierr);
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
  if (!surfelev) SETERRQ(1, "ERROR: surface_altitude is not available");

  lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
  if (!lat) SETERRQ(1, "ERROR: latitude is not available");

  PetscTruth gwl3_start_set;
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

PetscReal PA_EISMINT_Greenland::greenhouse_shift(PetscReal t_years, PetscReal dt_years) {
  // compute age back to start of GWL3; use midpoint of interval as the time
  PetscScalar age = (t_years + 0.5 * dt_years) - greenhouse_warming_start_year;
  if (age <= 0.0) {
    return 0.0; // before time 0.0, no warming
  } else {
    if (age <= 80.0) {
      return age * 0.035;
    } else if (age <= 500.0) {
      return 2.8 + (age - 80.0) * 0.0017;
    } else { // after 500 years, constant amount
      return 3.514;
    }
  }
}

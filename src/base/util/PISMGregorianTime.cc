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

#include "PISMGregorianTime.hh"

PISMGregorianTime::PISMGregorianTime(MPI_Comm c, const NCConfigVariable &conf)
  : PISMTime(c, conf) {

  calendar_string = "gregorian";  // only "gregorian" is supported by this class

  run_start = config.get("start_year", "years", "seconds");
  run_end   = run_start + config.get("run_length_years", "years", "seconds");

  time_in_seconds = run_start;
}

PetscErrorCode PISMGregorianTime::init() {
  PetscErrorCode ierr;

  ierr = PetscOptionsBegin(grid.com, "", "PISM time options", ""); CHKERRQ(ierr);
  {
    // -start_date
    // -end_date
    // -run_length
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMGregorianTime::interval_to_seconds(string interval, double &result) {
  PetscErrorCode ierr;
  utUnit ut_unit, internal_unit;
  double slope, intercept;
  int errcode;

  errcode = utScan("seconds", &internal_unit);
  if (errcode != 0)
    SETERRQ(PETSC_COMM_SELF, 1, "utScan(\"seconds\", ...) failed");

  errcode = utScan(interval.c_str(), &ut_unit);
  if (errcode != 0) {
    PetscPrintf(com, "PISM ERROR: can't parse '%s'\n", interval.c_str());
    PISMEnd();
  }

  if (utIsTime(&ut_unit) == 1) {
    errcode = utConvert(&ut_unit, &internal_unit, &slope, &intercept);
    if (errcode != 0) {
      PetscPrintf(com, "PISM ERROR: can't convert '%s' to seconds\n", interval.c_str());
      PISMEnd();
    }

    result = slope;
  } else {
    errcode = utScan("1", &internal_unit);
    if(errcode != 0)
      SETERRQ(PETSC_COMM_SELF, 1, "utScan(\"1\", ...) failed");

    errcode = utConvert(&ut_unit, &internal_unit, &slope, &intercept);
    if (errcode != 0) {
      PetscPrintf(com, "PISM ERROR: can't parse '%s'\n", interval.c_str());
      PISMEnd();
    }

    result = slope;
  }

  return 0;
}


// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

#include <cstring>
#include <cstdio>
#include <petscda.h>
#include "iceModel.hh"


bool IceModel::hasSuffix(const char* fname, const char *suffix) const {
  int flen = strlen(fname);
  int slen = strlen(suffix);
  if (strcmp(fname + flen - slen, suffix) == 0) {
    return true;
  } else {
    return false;
  }
}


//! Initialize from a saved PISM model state (in NetCDF format).
/*! 
Calls initFromFile_netCDF() to do the actual work.
 */
PetscErrorCode IceModel::initFromFile(const char *fname) {
  PetscErrorCode  ierr;

  if (hasSuffix(fname, ".pb") == true) {
    SETERRQ1(1,"ERROR: .pb format no longer supported; cannot initialize from file %s", fname);
  }
  
  if (hasSuffix(fname, ".nc") == false) {
    ierr = verbPrintf(1,grid.com,
       "WARNING:  Unknown file format for %s.  Trying to read as NetCDF.\n",fname); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2,grid.com,"initializing from NetCDF format file  %s  ...\n",
                     fname); CHKERRQ(ierr);
  ierr = initFromFile_netCDF(fname); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModel::setStartRunEndYearsFromOptions(const PetscTruth grid_p_year_VALID) {
  PetscErrorCode ierr;

  // read options about year of start, year of end, number of run years;
  // note grid.p->year has already been set from input file
  PetscScalar usrStartYear, usrEndYear, usrRunYears;
  PetscTruth ysSet, yeSet, ySet;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ys", &usrStartYear, &ysSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ye", &usrEndYear, &yeSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-y", &usrRunYears, &ySet); CHKERRQ(ierr);
  if (ysSet == PETSC_TRUE) {
    // user option overwrites data
    ierr = setStartYear(usrStartYear); CHKERRQ(ierr);
    grid.p->year = usrStartYear;
  } else if (grid_p_year_VALID == PETSC_TRUE) {
    ierr = setStartYear(grid.p->year); CHKERRQ(ierr);
  } // else do nothing; defaults are set
  if (yeSet == PETSC_TRUE) {
    if (usrEndYear < startYear) {
      SETERRQ(1,
        "ERROR: -ye value less than -ys value (or input file year or default).\n"
        "PISM cannot run backward in time");
    }
    if (ySet == PETSC_TRUE) {
      ierr = verbPrintf(1,grid.com,"WARNING: -y option ignored.  -ye used instead.\n"); CHKERRQ(ierr);
    }
    ierr = setEndYear(usrEndYear); CHKERRQ(ierr);
  } else if (ySet == PETSC_TRUE) {
    ierr = setEndYear(usrRunYears + startYear); CHKERRQ(ierr);
  } else {
    ierr = setEndYear(DEFAULT_RUN_YEARS + startYear); CHKERRQ(ierr);
  }
  
  yearsStartRunEndDetermined = PETSC_TRUE;
  return 0;
}

  
PetscErrorCode  IceModel::writeFiles(const char* defaultbasename) {
  PetscErrorCode ierr = writeFiles(defaultbasename,PETSC_FALSE); CHKERRQ(ierr);
  return 0;
}

  
//! Save model state in NetCDF format (and save variables in Matlab format if desired).
/*! 
Optionally allows saving of full velocity field.

Calls dumpToFile_netCDF(), dumpToFile_diagnostic_netCDF(), and writeMatlabVars() to do 
the actual work.
 */
PetscErrorCode  IceModel::writeFiles(const char* defaultbasename,
                                     const PetscTruth forceFullDiagnostics) {
  PetscErrorCode ierr;

  char b[PETSC_MAX_PATH_LEN];
  char fmt[PETSC_MAX_PATH_LEN] = "n";
  char ncf[PETSC_MAX_PATH_LEN]; // netCDF format

  if (doPDD == PETSC_TRUE) { // want to save snow accumulation map, not net accumulation
    ierr = putBackSnowAccumPDD(); CHKERRQ(ierr);
  }
  
  ierr = stampHistoryEnd(); CHKERRQ(ierr);

  // Use the defaults passed from the driver if not specified on command line.
  // We should leave space for a suffix and null byte
  strncpy(b, defaultbasename, PETSC_MAX_PATH_LEN-4);
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", b, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-of", fmt, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);

  if (strchr(fmt, 'p') != NULL) {
    strcat(b,"_pb");  // will write basename_pb.nc
    strcat(fmt,"n");
    ierr = verbPrintf(1, grid.com, 
       "\nWARNING: .pb format no longer supported; writing to NetCDF file %s.nc\n",b);
       CHKERRQ(ierr);
  }

  if (strchr(fmt, 'n') != NULL) {
    strcpy(ncf, b);
    strcat(ncf, ".nc");
    PetscTruth userWantsFull;
    ierr = PetscOptionsHasName(PETSC_NULL, "-full3Dout", &userWantsFull); CHKERRQ(ierr);
    if ((forceFullDiagnostics == PETSC_TRUE) || (userWantsFull == PETSC_TRUE)) {
      ierr = verbPrintf(1, grid.com, 
            "Writing model state, with full 3D velocities, to file `%s'", ncf); CHKERRQ(ierr);
      ierr = dumpToFile_diagnostic_netCDF(ncf); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(1, grid.com, "Writing model state to file `%s'", ncf); CHKERRQ(ierr);
      ierr = dumpToFile_netCDF(ncf); CHKERRQ(ierr);
    }
  }

  if (strchr(fmt, 'm') != NULL) {
    ierr = verbPrintf(1, grid.com, 
       "\nWARNING: .m format no longer supported with '-o'; use '-mato' and '-matv'\n");
       CHKERRQ(ierr);
  }

  // write out individual variables out to Matlab file
  char       matf[PETSC_MAX_PATH_LEN];
  PetscTruth matoSet, matvSet;
  ierr = PetscOptionsGetString(PETSC_NULL, "-mato", matf, PETSC_MAX_PATH_LEN, &matoSet); 
           CHKERRQ(ierr);
  if (matoSet == PETSC_FALSE) {// put default name in matf; perhaps user set "-matv" only
    strcpy(matf, "pism_views");
  }
  strcpy(matlabOutVars, "\0");
  ierr = PetscOptionsGetString(PETSC_NULL, "-matv", matlabOutVars, PETSC_MAX_PATH_LEN, &matvSet); 
            CHKERRQ(ierr);
  if (matvSet == PETSC_TRUE) {
    strcat(matf, ".m");
    ierr = verbPrintf(1, grid.com, 
       " ... writing variables %s to Matlab file `%s'", matlabOutVars, matf); CHKERRQ(ierr);
    ierr = writeMatlabVars(matf); CHKERRQ(ierr); // see iMmatlab.cc
  }

  return 0;
}


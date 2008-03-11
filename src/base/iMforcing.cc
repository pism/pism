// Copyright (C) 2008 Nathan Shemonski and Ed Bueler
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
#include "forcing.hh"
#include "iceModel.hh"


//! Initialize the forcing data if either -dTforcing or -dSLforcing is set.
PetscErrorCode IceModel::initForcingFromOptions() {
  PetscErrorCode ierr;
  char dTFile[PETSC_MAX_PATH_LEN], dSLFile[PETSC_MAX_PATH_LEN];
  PetscTruth dTforceSet, dSLforceSet;

  if (dTforcing != PETSC_NULL) {
    SETERRQ(1, "dTforcing should be PETSC_NULL at start of initForcingFromOptions()\n");
  }
  if (dSLforcing != PETSC_NULL) {
    SETERRQ(2, "dSLforcing should be PETSC_NULL at start of initForcingFromOptions()\n");
  }
  ierr = PetscOptionsGetString(PETSC_NULL, "-dTforcing", dTFile,
                               PETSC_MAX_PATH_LEN, &dTforceSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-dSLforcing", dSLFile,
                               PETSC_MAX_PATH_LEN, &dSLforceSet); CHKERRQ(ierr);

  if (dTforceSet == PETSC_TRUE) {
    dTforcing = new IceSheetForcing;
    int stat, ncid = 0;
    ierr = verbPrintf(2, grid.com, 
         "reading delta T data from forcing file %s ...\n", dTFile); 
         CHKERRQ(ierr);
    if (grid.rank == 0) {
      stat = nc_open(dTFile, 0, &ncid); CHKERRQ(nc_check(stat));
    }
    ierr = dTforcing->readCoreClimateData(grid.com, grid.rank, ncid, 
         grid.year,ISF_DELTA_T); CHKERRQ(ierr);
    if (grid.rank == 0) {
      stat = nc_close(ncid); CHKERRQ(nc_check(stat));
    }
  }

  if (dSLforceSet == PETSC_TRUE) {
    bedDiffSeaLevel = 0.0;
    dSLforcing = new IceSheetForcing;
    int stat, ncid = 0;
    ierr = verbPrintf(2, grid.com, 
         "reading delta sea level data from forcing file %s ...\n", 
         dSLFile); CHKERRQ(ierr);
    if (grid.rank == 0) {
      stat = nc_open(dSLFile, 0, &ncid); CHKERRQ(nc_check(stat));
    }
    ierr = dSLforcing->readCoreClimateData(grid.com, grid.rank, 
             ncid, grid.year,ISF_DELTA_SEA_LEVEL); CHKERRQ(ierr);
    if (grid.rank == 0) {
      stat = nc_close(ncid); CHKERRQ(nc_check(stat));
    }
  }

  return 0;
}


//! Read in new delta T or delta sea level from forcing data.
PetscErrorCode IceModel::updateForcing() {
  PetscErrorCode ierr;

  if (dTforcing != PETSC_NULL) {
    PetscScalar TsChange;
    ierr = dTforcing->updateFromCoreClimateData(grid.year,&TsChange); CHKERRQ(ierr);
    ierr = VecShift(vTs,TsChange); CHKERRQ(ierr);
  }

  if (dSLforcing != PETSC_NULL) {
    PetscScalar seaLevelChange;
    ierr = dSLforcing->updateFromCoreClimateData(grid.year,&seaLevelChange); CHKERRQ(ierr);
    ierr = VecShift(vbed,bedDiffSeaLevel - seaLevelChange); CHKERRQ(ierr);
    bedDiffSeaLevel = seaLevelChange;
    // communicate ghosted bed (bed slope matters); needed with VecShift?
    ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
    if (seaLevelChange != 0) {
      updateSurfaceElevationAndMask();
    }
  }
  return 0;
}


//! De-allocate the forcing stuff.
PetscErrorCode IceModel::forcingCleanup() {
  PetscErrorCode ierr;
  if (dTforcing != PETSC_NULL) {
    delete dTforcing;
    dTforcing = PETSC_NULL;
  }
  if (dSLforcing != PETSC_NULL) {
    ierr = VecShift(vbed, bedDiffSeaLevel); CHKERRQ(ierr);
    delete dSLforcing;
    dSLforcing = PETSC_NULL;
  }
  return 0;
}


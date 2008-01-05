// Copyright (C) 2008 Ed Bueler
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
#include <cstdlib>
#include <petscda.h>
#include <netcdf.h>
#include "nc_util.hh"


IceModelVec::IceModelVec() {

  allocated = PETSC_FALSE;
  v = PETSC_NULL;
  da = PETSC_NULL;
  grid = PETSC_NULL;

  strcpy(long_name,"UNKNOWN long_name");
  strcpy(units,"UNKNOWN units");
  strcpy(pism_intent,"UNKNOWN pism_intent");

  has_standard_name = PETSC_FALSE;
  strcpy(standard_name,"UNKNOWN standard_name");
  
  ncid = -9999;
};


IceModelVec::~IceModelVec() {

  if (allocated == PETSC_TRUE) {
    destroy();
    allocated = PETSC_FALSE;
  }
};


PetscErrorCode  IceModelVec::destroy() {
  PetscErrorCode ierr = VecDestroy(v); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec::checkAllocated() {
  if (allocated == PETSC_FALSE) {
    SETERRQ1(1,"IceModelVec ERROR: IceModelVec NOT allocated (long_name = %s)\n",
             long_name);
  }
  return 0;
}


/*
PetscErrorCode  IceModelVec::checkSelfOwnsIt(const PetscInt i, const PetscInt j) {
  if (allocated == PETSC_FALSE) {
    SETERRQ3(1,"IceModelVec ERROR: (i,j)=(%d,%d) not in ownership range of processor %d\n",
             i,j,grid->rank);
  }
  return 0;
}
*/

PetscErrorCode  IceModelVec::beginGhostComm() {
  CHKERRQ(checkAllocated());
  PetscErrorCode ierr;
  ierr = DALocalToLocalBegin(da, v, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec::endGhostComm() {
  CHKERRQ(checkAllocated());
  PetscErrorCode ierr;
  ierr = DALocalToLocalEnd(da, v, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}



/****************** 3d Vecs ********************/

IceModelVec3::IceModelVec3() : IceModelVec() {
  Mz = -1;
  dz = -1.0;
  levels = PETSC_NULL;
  
  array = PETSC_NULL;
};


PetscErrorCode  IceModelVec3::destroy() {
  if (Mz > 0) {
    delete [] levels;
  }
  PetscErrorCode ierr = IceModelVec::destroy(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::create(IceGrid* mygrid);
  if (allocated == PETSC_TRUE) {
    SETERRQ1(1,"IceModelVec3 ERROR: IceModelVec3 already allocated; has long_name = %s\n",
             long_name);
  }
  
  grid = mygrid;
  da = grid->da3;
  Mz = (grid->p)->Mz;
  dz = (grid->p)->dz;

  levels = new PetscScalar[Mz];
  for (PetscInt k=0; k < Mz; k++) {
    levels[k] = dz * ((PetscScalar) k);
  }
  
  PetscErrorCode ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  allocated = PETSC_TRUE;
  return 0;
}


PetscErrorCode  IceModelVec3::needAccessToVals() {
  CHKERRQ(checkAllocated());
  PetscErrorCode ierr = DAVecGetArray(da, v, &array); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::doneAccessToVals() {
  CHKERRQ(checkAllocated());
  PetscErrorCode ierr = DAVecRestoreArray(da, v, &array); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::checkHaveArray() {
  CHKERRQ(checkAllocated());
  if (array == PETSC_NULL) {
    SETERRQ(1,"IceModelVec3 ERROR: array not available\n"
              "  (REMEMBER TO RUN needAccessToVals() before access and doneAccessToVals after access)\n");
  }
  return 0;
}


PetscErrorCode  IceModelVec3::isLegalLevel(const PetscScalar z) {
  if (z < 0.0) {
    SETERRQ(1,"IceModelVec3 ERROR: level z is below base of ice (z must be nonnegative); ENDING!\n");
  }
  if (z > levels[Mz - 1]) {
    SETERRQ(2,"IceModelVec3 ERROR: level z is above top of computational grid for ice; ENDING!\n");
  }
  return 0;
}


//! Compute value of scalar quantity at level z (m) above base of ice.
PetscScalar     IceModelVec3::getValZ(const PetscInt i, const PetscInt j, const PetscScalar z) {
  // use linear interpolation
  checkHaveArray();
  isLegalLevel(z);
  if (z == levels[Mz - 1])
    return array[i][j][Mz - 1];
  const PetscInt     kbz = static_cast<PetscInt>(floor(z/dz));  // k value just below z
  const PetscScalar  val_kbz = array[i][j][kbz];
  return val_kbz + ( (z - levels[kbz]) / dz ) * (array[i][j][kbz + 1] - val_kbz);
}


PetscErrorCode  IceModelVec3::getHorSlice(Vec &gslice, const PetscScalar z) {
  PetscErrorCode ierr;
  PetscScalar    **slice_val;
  ierr = DAVecGetArray(grid->da2, gslice, &slice_val); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      slice_val[i][j] = getValZ(i,j,z);
    }
  }
  ierr = DAVecRestoreArray(grid->da2, gslice, &slice_val); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::getSurfaceValues(Vec &gsurf, PetscScalar **H) {
  PetscErrorCode ierr;
  PetscScalar    **surf_val;
  ierr = DAVecGetArray(grid->da2, gsurf, &surf_val); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      surf_val[i][j] = getValZ(i,j,H[i][j]);
    }
  }
  ierr = DAVecRestoreArray(grid->da2, gsurf, &surf_val); CHKERRQ(ierr);
  return 0;
}


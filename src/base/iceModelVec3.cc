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

#include "iceModelVec.hh"

// this file contains method for derived class IceModelVec3

// methods for base class IceModelVec and derived class IceModelVec2 and IceModelVec2Box 
// are in "iceModelVec.cc"

IceModelVec3::IceModelVec3() : IceModelVec() {
};


//! Allocate a DA and a Vec from information in IceGrid.
PetscErrorCode  IceModelVec3::create(IceGrid &mygrid, const char my_varname[], bool local) {
  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec3 with varname='%s' already allocated\n",varname);
  }
  
  grid = &mygrid;
  
  PetscInt       M, N, m, n;
  PetscErrorCode ierr;
  ierr = DAGetInfo(mygrid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate3d(mygrid.com, DA_YZPERIODIC, DA_STENCIL_STAR, mygrid.p->Mz, N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da); CHKERRQ(ierr);
  IOwnDA = true;

  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  strcpy(varname,my_varname);
  return 0;
}


//! Allocate a DA and a Vec from information in IceGrid; use an existing DA from an existing IceModelVec3.
PetscErrorCode  IceModelVec3::createSameDA(IceModelVec3 imv3_source,
                                           IceGrid &mygrid, const char my_varname[], bool local) {
  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec3 with varname='%s' already allocated\n",varname);
  }
  
  grid = &mygrid;
  
  da = imv3_source.da;
  IOwnDA = false;

  PetscErrorCode ierr;
  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  strcpy(varname,my_varname);
  return 0;
}


PetscErrorCode  IceModelVec3::beginGhostCommTransfer(IceModelVec3 imv3_source) {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec3! (has varname='%s')\n",
               varname);
  }
//  if (!imv3_source.localp) {
//    SETERRQ1(2,"makes no sense to communicate ghosts for GLOBAL IceModelVec! (has varname='%s')\n",
//               imv3_source.varname);
//  }
  if (imv3_source.localp) {
    SETERRQ1(2,"source IceModelVec3 must be GLOBAL! (has varname='%s')\n",
               imv3_source.varname);
  }
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = imv3_source.checkAllocated(); CHKERRQ(ierr);
//  ierr = DALocalToLocalBegin(da, imv3_source.v, INSERT_VALUES, v);  CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, imv3_source.v, INSERT_VALUES, v);  CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::endGhostCommTransfer(IceModelVec3 imv3_source) {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec3! (has varname='%s')\n",
               varname);
  }
//  if (!imv3_source.localp) {
//    SETERRQ1(2,"makes no sense to communicate ghosts for GLOBAL IceModelVec! (has varname='%s')\n",
//               imv3_source.varname);
//  }
  if (imv3_source.localp) {
    SETERRQ1(2,"source IceModelVec3 must be GLOBAL! (has varname='%s')\n",
               imv3_source.varname);
  }
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = imv3_source.checkAllocated(); CHKERRQ(ierr);
//  ierr = DALocalToLocalEnd(da, imv3_source.v, INSERT_VALUES, v); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, imv3_source.v, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::isLegalLevel(const PetscScalar z) {
  if (z < 0.0 - 1.0e-6) {
    SETERRQ2(1,"level z = %5.4f is below base of ice (z must be nonnegative);\n"
               "  IceModelVec3 has varname='%s'; ENDING!\n",
              z,varname);
  }
  if (z > (grid->p)->Lz + 1.0e-6) {
    SETERRQ3(2,"level z = %10.8f is above top of computational grid Lz = %10.8f;\n"
               "  IceModelVec3 has varname='%s'; ENDING!\n",
              z, (grid->p)->Lz,varname);
  }
  return 0;
}


//! Set values of an ice scalar quantity by linear <i>interpolation</i> from given values in a given column.
/*!
Input arrays \c levelsIN and \c valsIN must be allocated arrays of \c nlevels scalars
(\c PetscScalar).  Upon completion, internal storage will hold values derived from 
linearly interpolating the input values.

\c levelsIN must be strictly increasing.

Piecewise linear interpolation is used and the input values must span a sufficient range
of \f$z\f$ values so that all stored values, at heights in \c zlevels, can be determined 
by interpolation; extrapolation is not allowed.  Therefore <tt>(levelsIN[0] <= 0.0)</tt> 
and <tt>(levelsIN[nlevels-1] >= Lz)</tt> must both be true.
 */
PetscErrorCode  IceModelVec3::setValColumn(
                     const PetscInt i, const PetscInt j, 
                     const PetscInt nlevels, PetscScalar *levelsIN,
                     PetscScalar *valsIN) {

  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // check if in ownership ?

  if (levelsIN[0] > 0.0 + 1.0e-3) {
    SETERRQ2(1,"levelsIN[0]=%10.9f is above base of ice at z=0 so *interpolation*\n"
              "   is impossible; IceModelVec3 has varname='%s';  ENDING!\n",
              levelsIN[0],varname);
  }
  if (levelsIN[nlevels - 1] < (grid->p)->Lz - 1.0e-3) {
    SETERRQ3(2,"levelsIN[nlevels-1] = %10.9f is below top of computational domain\n"
               "   at z=Lz=%10.9f, so *interpolation* is impossible;\n"
               "   IceModelVec3 has varname='%s';  ENDING!\n",
               levelsIN[nlevels-1],(grid->p)->Lz,varname);
  }
  for (PetscInt k=0; k < nlevels - 1; k++) {
    if (levelsIN[k] >= levelsIN[k+1]) {
      SETERRQ2(3,"levelsIN not *strictly increasing* at index %d;\n"
                 "    IceModelVec3 has varname='%s';  ENDING!\n",
                 k,varname);
    }
  }

  PetscScalar *levels;
  levels = grid->zlevels;
  PetscScalar ***arr = (PetscScalar***) array;
  
  PetscInt mcurr = 0;
  for (PetscInt k=0; k < (grid->p)->Mz; k++) {
    while (levelsIN[mcurr+1] < levels[k]) {
      mcurr++;
    }
    const PetscScalar increment = (levels[k] - levelsIN[mcurr]) / (levelsIN[mcurr+1] - levelsIN[mcurr]);
    arr[i][j][k] = valsIN[mcurr] +  increment * (valsIN[mcurr+1] - valsIN[mcurr]);
  }

  return 0;
}


//! Set all values of scalar quantity to given a single value in a particular column.
PetscErrorCode  IceModelVec3::setToConstantColumn(const PetscInt i, const PetscInt j, 
                                                  const PetscScalar c) {

  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  // check if in ownership ?
  PetscScalar ***arr = (PetscScalar***) array;
  for (PetscInt k=0; k < (grid->p)->Mz; k++) {
    arr[i][j][k] = c;
  }
  return 0;
}


//! Return value of scalar quantity at level z (m) above base of ice (by interpolation).
PetscScalar     IceModelVec3::getValZ(const PetscInt i, const PetscInt j, const PetscScalar z) {
  // use linear interpolation
  if (checkHaveArray() != 0) {
    PetscPrintf(PETSC_COMM_SELF, 
       "IceModelVec3 getValZ(): array was not allocated (so says IceModelVec::checkHaveArray())\n");
    PetscEnd();
  }
  if (isLegalLevel(z) != 0) {
    PetscPrintf(PETSC_COMM_SELF, 
       "IceModelVec3 getValZ(): level was not equal (so says isLegalLevel())\n");
    PetscEnd();
  }
  PetscScalar ***arr = (PetscScalar***) array;
  if (z >= (grid->p)->Lz)
    return arr[i][j][(grid->p)->Mz - 1];
  else if (z <= 0.0)
    return arr[i][j][0];

  PetscScalar* levels = grid->zlevels;
  PetscInt mcurr = 0;
  while (levels[mcurr+1] < z) {
    mcurr++;
  }
  const PetscScalar incr = (z - levels[mcurr]) / (levels[mcurr+1] - levels[mcurr]);
  const PetscScalar valm = arr[i][j][mcurr];
  return valm + incr * (arr[i][j][mcurr+1] - valm);
}


//! Return values on planar star stencil of scalar quantity at level z.
PetscErrorCode   IceModelVec3::getPlaneStarZ(const PetscInt i, const PetscInt j, const PetscScalar z,
                                             planeStar *star) {
  PetscErrorCode ierr;
  // use linear interpolation
  ierr = checkHaveArray();  CHKERRQ(ierr);
  ierr = isLegalLevel(z);  CHKERRQ(ierr);
  // check ownership here?
  if (!localp) {
    SETERRQ1(1,"IceModelVec3 ERROR: IceModelVec3 with varname='%s' is GLOBAL and cannot do getPlaneStarZ()\n",
             varname);
  }

  PetscInt     kbz;
  PetscScalar  incr;
  if (z >= (grid->p)->Lz) {
    kbz = (grid->p)->Mz - 1;
    incr = 0.0;
  } else if (z <= 0.0) {
    kbz = 0;
    incr = 0.0;
  } else {
    PetscScalar* levels = grid->zlevels;
    kbz = 0;
    while (levels[kbz+1] < z) {
      kbz++;
    }
    incr = (z - levels[kbz]) / (levels[kbz+1] - levels[kbz]);
  }

  PetscScalar ***arr = (PetscScalar***) array;

  star->ij = arr[i][j][kbz] + incr * (arr[i][j][kbz + 1] - arr[i][j][kbz]);
  star->ip1 = arr[i+1][j][kbz] + incr * (arr[i+1][j][kbz + 1] - arr[i+1][j][kbz]);
  star->im1 = arr[i-1][j][kbz] + incr * (arr[i-1][j][kbz + 1] - arr[i-1][j][kbz]);
  star->jp1 = arr[i][j+1][kbz] + incr * (arr[i][j+1][kbz + 1] - arr[i][j+1][kbz]);
  star->jm1 = arr[i][j-1][kbz] + incr * (arr[i][j-1][kbz + 1] - arr[i][j-1][kbz]);
  return 0;
}


//! Return values of ice scalar quantity at given levels (m) above base of ice, using piecewise linear interpolation.
/*!
Input array \c levelsIN must be an allocated array of \c nlevels scalars (\c PetscScalar).

\c levelsIN must be strictly increasing and in the range \f$0 <= z <= \mathtt{grid.p->Lz}\f$.

Return array \c valsOUT must be an allocated array of \c nlevels scalars (\c PetscScalar).
Upon return, \c valsOUT will be filled with values of scalar quantity at the \f$z\f$ values in \c levelsIN.
 */
PetscErrorCode  IceModelVec3::getValColumn(const PetscInt i, const PetscInt j, 
                     const PetscInt nlevelsIN, PetscScalar *levelsIN, PetscScalar *valsOUT) {
  
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // check if in ownership ?

  ierr = isLegalLevel(levelsIN[0]); CHKERRQ(ierr);
  ierr = isLegalLevel(levelsIN[nlevelsIN - 1]); CHKERRQ(ierr);
  for (PetscInt k=0; k < nlevelsIN - 1; k++) {
    if (levelsIN[k] >= levelsIN[k+1]) {
      SETERRQ2(1,"levelsIN not *strictly increasing* at index %d\n"
                 "    (IceModelVec3 with varname='%s')  ENDING!\n",k,varname);
    }
  }

  PetscScalar* levels = grid->zlevels;
  PetscScalar ***arr = (PetscScalar***) array;
  
  PetscInt mcurr = 0;
  for (PetscInt k = 0; k < nlevelsIN; k++) {
    while (levels[mcurr+1] < levelsIN[k]) {
      mcurr++;
    }
    const PetscScalar incr = (levelsIN[k] - levels[mcurr]) / (levels[mcurr+1] - levels[mcurr]);
    const PetscScalar valm = arr[i][j][mcurr];
    valsOUT[k] = valm + incr * (arr[i][j][mcurr+1] - valm);
  }

  return 0;
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


PetscErrorCode  IceModelVec3::getSurfaceValuesVec2d(Vec &gsurf, Vec myH) {
  PetscErrorCode ierr;
  PetscScalar    **H;
  ierr = DAVecGetArray(grid->da2, myH, &H); CHKERRQ(ierr);
  ierr = getSurfaceValuesArray2d(gsurf, H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid->da2, myH, &H); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::getSurfaceValuesArray2d(Vec &gsurf, PetscScalar **H) {
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


PetscErrorCode  IceModelVec3::getInternalColumn(
                     const PetscInt i, const PetscInt j, PetscInt *nlevels,
                     PetscScalar **levelsOUT, PetscScalar **valsOUT) {
  
  *nlevels = (grid->p)->Mz;
  *levelsOUT = grid->zlevels;
  PetscScalar ***arr = (PetscScalar***) array;
  *valsOUT = arr[i][j];
  return 0;
}


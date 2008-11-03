// Copyright (C) 2008 Ed Bueler and Constantine Khroulev
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

// this file contains methods for derived classes IceModelVec2 and IceModelVec2Box

// methods for base class IceModelVec are in "iceModelVec.cc"

IceModelVec2::IceModelVec2() : IceModelVec() {}


PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_varname[], bool local) {

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2 with varname='%s' already allocated\n",my_varname);
  }
  PetscErrorCode ierr = create(my_grid, my_varname, local, DA_STENCIL_STAR); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode  IceModelVec2::createSameDA(IceModelVec2 imv2_source,
					   IceGrid &my_grid, const char my_varname[], bool local) {

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2 with varname='%s' already allocated\n",my_varname);
  }

  grid = &my_grid;

  da = imv2_source.da;
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

  

PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_varname[], bool local,
                                     DAStencilType my_sten) {

  grid = &my_grid;
  
  PetscInt       M, N, m, n;
  PetscErrorCode ierr;
  ierr = DAGetInfo(my_grid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate2d(my_grid.com, DA_XYPERIODIC, my_sten, N, M, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, &da); CHKERRQ(ierr);
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


// Return value of ice scalar quantity stored in an IceModelVec2.
PetscScalar     IceModelVec2::getVal(const PetscInt i, const PetscInt j) {
  checkHaveArray();
  PetscScalar **arr = (PetscScalar**) array;
  return arr[i][j];
}


// Return values on planar star stencil of ice scalar quantity stored in an IceModelVec2.
PetscErrorCode   IceModelVec2::getPlaneStar(const PetscInt i, const PetscInt j, planeStar *star) {
  PetscErrorCode ierr;
  ierr = checkHaveArray();  CHKERRQ(ierr);
  // check ownership here?
  if (!localp) {
    SETERRQ1(1,"IceModelVec2 ERROR: IceModelVec2 with varname='%s' is GLOBAL and cannot do getPlaneStar()\n",
             varname);
  }
  
  PetscScalar **arr = (PetscScalar**) array;

  star->ij = arr[i][j];
  star->ip1 = arr[i+1][j];
  star->im1 = arr[i-1][j];
  star->jp1 = arr[i][j+1];
  star->jm1 = arr[i][j-1];
  return 0;
}


PetscScalar**   IceModelVec2::arrayGet() {
  checkHaveArray();
  return (PetscScalar**) array;
}



/********* IceModelVec2Box **********/

IceModelVec2Box::IceModelVec2Box() : IceModelVec2() {}


PetscErrorCode  IceModelVec2Box::create(IceGrid &my_grid, const char my_varname[], bool local) {

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2Box with varname='%s' already allocated\n",my_varname);
  }
  PetscErrorCode ierr = IceModelVec2::create(my_grid, my_varname, local, DA_STENCIL_BOX); CHKERRQ(ierr);
  return 0;
}


// Return values on planar BOX stencil of ice scalar quantity stored in an IceModelVec2Box.
PetscErrorCode   IceModelVec2Box::getPlaneBox(const PetscInt i, const PetscInt j, planeBox *box) {
  PetscErrorCode ierr;
  ierr = checkHaveArray();  CHKERRQ(ierr);
  // check ownership here?
  if (!localp) {
    SETERRQ1(1,"IceModelVec2Box ERROR: IceModelVec2Box with varname='%s' is GLOBAL and cannot do getPlaneBox()\n",
             varname);
  }

  PetscScalar **arr = (PetscScalar**) array;

  box->ij     = arr[i][j];
  box->ip1    = arr[i+1][j];
  box->im1    = arr[i-1][j];
  box->jp1    = arr[i][j+1];
  box->jm1    = arr[i][j-1];
  box->ip1jp1 = arr[i+1][j+1];
  box->im1jp1 = arr[i-1][j+1];
  box->ip1jm1 = arr[i+1][j-1];
  box->im1jm1 = arr[i-1][j-1];

  return 0;
}


/********* IceModelVec3Bedrock **********/

IceModelVec3Bedrock::IceModelVec3Bedrock() : IceModelVec() {}


//! Allocate a DA and a Vec from information in IceGrid.
PetscErrorCode  IceModelVec3Bedrock::create(IceGrid &my_grid, 
                               const char my_varname[], bool local) {

  strcpy(varname,my_varname);

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec3Bedrock with varname='%s' already allocated\n",varname);
  }
  if (local) {
    SETERRQ1(2,"IceModelVec3Bedrock must be GLOBAL (varname='%s')\n",varname);
  }

  grid = &my_grid;
  
  PetscInt       M, N, m, n;
  PetscErrorCode ierr;
  ierr = DAGetInfo(my_grid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate3d(my_grid.com, DA_YZPERIODIC, DA_STENCIL_STAR, my_grid.Mbz, 
                    N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da); CHKERRQ(ierr);
  IOwnDA = true;

  ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);

  localp = false;
  return 0;
}


//! Set values of bedrock scalar quantity at internal levels determined by IceGrid.
/*!
Array \c valsIN must be an allocated array of \c grid->Mbz \c PetscScalar s.
 */
PetscErrorCode  IceModelVec3Bedrock::setInternalColumn(
                   const PetscInt i, const PetscInt j, PetscScalar *valsIN) {
  
  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  PetscScalar ***arr = (PetscScalar***) array;
  for (PetscInt k = 0; k < grid->Mbz; k++) {
    arr[i][j][k] = valsIN[k];
  }
  return 0;
}


//! Set values of bedrock scalar quantity: set all values in a column to the same value.
PetscErrorCode  IceModelVec3Bedrock::setToConstantColumn(
                        const PetscInt i, const PetscInt j, const PetscScalar c) {

  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  PetscScalar ***arr = (PetscScalar***) array;
  for (PetscInt k = 0; k < grid->Mbz; k++) {
    arr[i][j][k] = c;
  }
  return 0;
}


//! Return values of bedrock scalar quantity at internal levels determined by IceGrid.
/*!
Return array \c valsOUT is an allocated array of \c grid->Mbz \c PetscScalar s.
 */
PetscErrorCode  IceModelVec3Bedrock::getInternalColumn(const PetscInt i, const PetscInt j, 
                                                       PetscScalar **valsOUT) {
  
  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  PetscScalar ***arr = (PetscScalar***) array;
  *valsOUT = arr[i][j];
  return 0;
}


//! From given values, set a bedrock scalar quantity in a given column by piecewise linear interpolation.
/*!
Input arrays \c levelsIN and \c valsIN must be allocated arrays of \c nlevels scalars
(\c PetscScalar).  Upon completion, internal storage will hold values derived from 
linearly interpolating the input values.

\c levelsIN must be strictly increasing.

Piecewise linear interpolation is used and the input values must span a sufficient range
of \f$z\f$ values so that all stored values, at heights in \c zlevels, can be determined 
by interpolation; extrapolation is not allowed.  Therefore <tt>(levelsIN[0] <= -Lbz)</tt> 
and <tt>(levelsIN[nlevels-1] >= 0.0)</tt> must both be true.
 */
PetscErrorCode  IceModelVec3Bedrock::setValColumn(const PetscInt i, const PetscInt j, 
                     const PetscInt nlevels, PetscScalar *levelsIN, PetscScalar *valsIN) {

  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // check if in ownership ?

  if (levelsIN[0] > -grid->Lbz + 1.0e-3) {
    SETERRQ3(1,"levelsIN[0]=%10.9f is above base of bedrock at z=-%10.9f so *interpolation*\n"
              "   is impossible; IceModelVec3Bedrock has varname='%s';  ENDING!\n",
              levelsIN[0],grid->Lbz,varname);
  }
  if (levelsIN[nlevels - 1] < 0.0 - 1.0e-3) {
    SETERRQ2(2,"levelsIN[nlevels-1] = %10.9f is below z=0, so *interpolation* is impossible;\n"
               "   IceModelVec3Bedrock has varname='%s';  ENDING!\n",
               levelsIN[nlevels-1],varname);
  }
  for (PetscInt k=0; k < nlevels - 1; k++) {
    if (levelsIN[k] >= levelsIN[k+1]) {
      SETERRQ2(3,"levelsIN not *strictly increasing* at index %d;\n"
                 "    IceModelVec3Bedrock has varname='%s';  ENDING!\n",
                 k,varname);
    }
  }

  PetscScalar *levels;
  levels = grid->zblevels;
  PetscScalar ***arr = (PetscScalar***) array;
  
  PetscInt mcurr = 0;
  for (PetscInt k=0; k < grid->Mbz; k++) {
    while (levelsIN[mcurr+1] < levels[k]) {
      mcurr++;
    }
    const PetscScalar increment = (levels[k] - levelsIN[mcurr]) / (levelsIN[mcurr+1] - levelsIN[mcurr]);
    arr[i][j][k] = valsIN[mcurr] +  increment * (valsIN[mcurr+1] - valsIN[mcurr]);
  }
  return 0;
}


PetscErrorCode  IceModelVec3Bedrock::isLegalLevel(const PetscScalar z) {
  if (z < -grid->Lbz - 1.0e-6) {
    SETERRQ3(1,
       "level z = %10.8f is below bottom of bedrock at -Lbz = %10.8f; IceModelVec3Bedrock has varname='%s'; ENDING!\n",
       z,-grid->Lbz,varname);
  }
  if (z > 0.0 + 1.0e-6) {
    SETERRQ2(2,"level z = %10.8f is above top of bedrock at z=0; IceModelVec3Bedrock has varname='%s'; ENDING!\n",
              z,varname);
  }
  return 0;
}


//! At given levels, return values of bedrock scalar quantity in a given column using piecewise linear interpolation.
/*!
Input array \c levelsIN must be an allocated array of \c nlevels scalars (\c PetscScalar).

\c levelsIN must be strictly increasing and in the range \f$-\mathtt{grid.Lbz} <= z <= 0.0\f$.

Return array \c valsOUT must be an allocated array of \c nlevels scalars (\c PetscScalar).
Upon return, \c valsOUT will be filled with values of scalar quantity at the \f$z\f$ values in \c levelsIN.
 */
PetscErrorCode  IceModelVec3Bedrock::getValColumn(const PetscInt i, const PetscInt j, 
                     const PetscInt nlevelsIN, PetscScalar *levelsIN, PetscScalar *valsOUT) {
  
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // check if in ownership ?

  ierr = isLegalLevel(levelsIN[0]); CHKERRQ(ierr);
  ierr = isLegalLevel(levelsIN[nlevelsIN - 1]); CHKERRQ(ierr);
  for (PetscInt k=0; k < nlevelsIN - 1; k++) {
    if (levelsIN[k] >= levelsIN[k+1]) {
      SETERRQ2(1,"levelsIN not *strictly increasing* at index %d\n"
                 "    (IceModelVec3Bedrock with varname='%s')  ENDING!\n",k,varname);
    }
  }

  PetscScalar* levels = grid->zblevels;
  PetscScalar ***arr = (PetscScalar***) array;
  
  PetscInt mcurr = 0;
  for (PetscInt k = 0; k < nlevelsIN; k++) {
    while (levels[mcurr+1] < levelsIN[k]) {
      mcurr++;
    }
    const PetscScalar incr = (levelsIN[k] - levels[mcurr])
                             / (levels[mcurr+1] - levels[mcurr]);
    const PetscScalar valm = arr[i][j][mcurr];
    valsOUT[k] = valm + incr * (arr[i][j][mcurr+1] - valm);
  }

  return 0;
}

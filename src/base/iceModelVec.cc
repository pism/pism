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


IceModelVec::IceModelVec() {

  allocated = PETSC_FALSE;
  v = PETSC_NULL;
  da = PETSC_NULL;
  grid = PETSC_NULL;

  strcpy(varname_nc,"UNKNOWN NetCDF varname");
  strcpy(long_name,"UNKNOWN long_name");
  strcpy(units,"UNKNOWN units");
  strcpy(pism_intent,"UNKNOWN pism_intent");

  has_standard_name = PETSC_FALSE;
  strcpy(standard_name,"UNKNOWN standard_name");
  
  varid_nc = -9999;
};


IceModelVec::~IceModelVec() {

  destroy();
};


PetscErrorCode  IceModelVec::create(IceGrid &mygrid) {
  SETERRQ(1,"not implemented")
}


PetscErrorCode  IceModelVec::destroy() {
  if (allocated == PETSC_TRUE) {
    PetscErrorCode ierr = VecDestroy(v); CHKERRQ(ierr);
    allocated = PETSC_FALSE;
    v = PETSC_NULL;
  }
  return 0;
}


PetscErrorCode  IceModelVec::setVaridNC(const int my_varid) {
  varid_nc = my_varid;
  return 0;
}


PetscErrorCode  IceModelVec::setAttrsNC(const int my_varid, const char my_varname[],
             const char my_long_name[], const char my_units[], const char my_pism_intent[]) {
  varid_nc = my_varid;
  strcpy(varname_nc,my_varname);
  strcpy(long_name,my_long_name);
  strcpy(units,my_units);
  strcpy(pism_intent,my_pism_intent);
  has_standard_name = PETSC_FALSE;
  return 0;
}


PetscErrorCode  IceModelVec::setAttrsCFstandardNC(const int my_varid, const char my_varname[],
             const char my_long_name[], const char my_units[], const char my_pism_intent[],
             const char my_standard_name[]) {
  varid_nc = my_varid;
  strcpy(varname_nc,my_varname);
  strcpy(long_name,my_long_name);
  strcpy(units,my_units);
  strcpy(pism_intent,my_pism_intent);
  strcpy(standard_name,my_standard_name);
  has_standard_name = PETSC_TRUE;
  return 0;
}


PetscErrorCode  IceModelVec::writeAttrsNC(const int ncid) {
  SETERRQ(1,"not implemented")
}


PetscErrorCode  IceModelVec::readVecNC(const int ncid, PetscTruth *exists) {
  SETERRQ(1,"not implemented")
/*
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // on processor zero, use varname string to check if it is present
  if (grid->rank == 0) {
    int stat = nc_inq_varid(ncid, varname_nc, &varid_nc); 
    exists = (stat == NC_NOERR);
  }
  // broadcast the existence flag
  ierr = MPI_Bcast(&exists, 1, MPI_INT, 0, grid->com); CHKERRQ(ierr);
  // if it exists, read it
  if (exists == PETSC_TRUE) {
    ierr = ncVarToDAVec(ncid, varid_nc, *da, v, g2, vzero); CHKERRQ(ierr);
  }
*/
}


PetscErrorCode IceModelVec::putVecNC(const int ncid, Vec g, const int *s, const int *c, int dims, 
                                          void *a_mpi, int a_size) {
  PetscErrorCode ierr;
  if (grid->rank == 0) {
    ierr = put_local_var(grid, ncid, varid_nc, NC_FLOAT, *da, v, g,
                         s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
  }
  return 0;
}


PetscErrorCode  IceModelVec::needAccessToVals() {
  SETERRQ(1,"virtual; not implemented")
}


PetscErrorCode  IceModelVec::doneAccessToVals() {
  SETERRQ(1,"virtual; not implemented")
}


PetscErrorCode  IceModelVec::checkAllocated() {
  if (allocated == PETSC_FALSE) {
    SETERRQ1(1,"IceModelVec ERROR: IceModelVec NOT allocated (long_name = %s)\n",
             long_name);
  }
  return 0;
}


PetscErrorCode  IceModelVec::checkHaveArray() {
  SETERRQ(1,"virtual; not implemented")
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
  return DALocalToLocalBegin(*da, v, INSERT_VALUES, v);
}


PetscErrorCode  IceModelVec::endGhostComm() {
  CHKERRQ(checkAllocated());
  return DALocalToLocalEnd(*da, v, INSERT_VALUES, v);
}


PetscErrorCode  IceModelVec::setToConstant(const PetscScalar c) {
  CHKERRQ(checkAllocated());
  return VecSet(v,c);
}



/****************** 3d Vecs ********************/

IceModelVec3::IceModelVec3() : IceModelVec() {
  Mz = -1;
  dz = -1.0;
  Lz = -1.0;
  levels = PETSC_NULL;
  
  array = PETSC_NULL;
};


PetscErrorCode  IceModelVec3::destroy() {
  if (Mz > 0) {
    delete [] levels;
    levels = PETSC_NULL;
    Mz = -1;
    dz = -1.0;
    Lz = -1.0;
  }
  PetscErrorCode ierr = IceModelVec::destroy(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::create(IceGrid &mygrid) {
  if (allocated == PETSC_TRUE) {
    SETERRQ1(1,"IceModelVec3 ERROR: IceModelVec3 already allocated; has long_name = %s\n",
             long_name);
  }
  
  grid = &mygrid;
  da = &(grid->da3);
  Mz = (grid->p)->Mz;
  dz = (grid->p)->dz;
  Lz = (grid->p)->Lz;

  levels = new PetscScalar[Mz];
  for (PetscInt k=0; k < Mz; k++) {
    levels[k] = dz * ((PetscScalar) k);
  }
  
  if (PetscAbs(levels[Mz - 1] - Lz) > 1.0e-6) {
    SETERRQ(2,"IceModelVec3 ERROR: Lz does not agree with top level to within 10^-6 meters; why?");
  }
  levels[Mz - 1] = Lz;
  
  PetscErrorCode ierr = DACreateLocalVector(*da, &v); CHKERRQ(ierr);
  allocated = PETSC_TRUE;
  return 0;
}


PetscErrorCode  IceModelVec3::needAccessToVals() {
  CHKERRQ(checkAllocated());
  PetscErrorCode ierr = DAVecGetArray(*da, v, &array); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::doneAccessToVals() {
  CHKERRQ(checkAllocated());
  PetscErrorCode ierr = DAVecRestoreArray(*da, v, &array); CHKERRQ(ierr);
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
  if (z > Lz) {
    SETERRQ(2,"IceModelVec3 ERROR: level z is above top of computational grid for ice; ENDING!\n");
  }
  return 0;
}


//! Set values of scalar quantity by linear <i>interpolation</i> from given values in a given column.
/*!
Input arrays \c levelsIN and \c valsIN must be an allocated array of \c nlevels scalars
(i.e. \c PetscScalar).  Upon completion, internal storage will hold values derived from 
linearly interpolating the input values.

\c levelsIN must be strictly increasing.

Only interpolation is used.  Therefore <tt>(levelsIN[0] == 0.0)</tt> and <tt>(levelsIN[nlevels-1] >= Lz)</tt>
must both be true.
 */
PetscErrorCode  IceModelVec3::setValColumn(
                     const PetscInt i, const PetscInt j, 
                     const PetscInt nlevels, PetscScalar *levelsIN,
                     PetscScalar *valsIN) {

  if (levelsIN[0] > 0.0) {
    SETERRQ(1,"IceModelVec3 ERROR: levelsIN[0] is above base of ice at z=0 so *interpolation*\n"
              "   is impossible; ENDING!\n");
  }
  if (levelsIN[nlevels - 1] < Lz) {
    SETERRQ1(2,"IceModelVec3 ERROR: levelsIN[nlevels-1] is below top of computational domain\n"
               "   at z=Lz=%f, so *interpolation* is impossible; ENDING!\n",Lz);
  }
  for (PetscInt k=0; k < nlevels - 1; k++) {
    if (levelsIN[k] >= levelsIN[k+1]) {
      SETERRQ1(3,"IceModelVec3 ERROR: levelsIN not *strictly increasing* at index %d; ENDING!\n",k);
    }
  }

  PetscInt mcurr = 0;
  for (PetscInt k=0; k < Mz; k++) {
    while (levelsIN[mcurr+1] < levels[k]) {
      mcurr++;
    }
    if (mcurr+1 > nlevels - 1) {
      SETERRQ(4,"how did I get here?\n");
    }
    // check that [ levelsIN[mcurr], levelsIN[mcurr+1] ] is a bracket for levels[k]:
    if (!( (levelsIN[mcurr] <= levels[k]) && (levels[k] <= levelsIN[mcurr+1]) )) {
      SETERRQ2(5,"not a bracket; k=%d, mcurr=%d\n",k,mcurr);
    }
    const PetscScalar increment = (levels[k] - levelsIN[mcurr]) / (levelsIN[mcurr+1] - levelsIN[mcurr]);
    array[i][j][k] = valsIN[mcurr] +  increment * (valsIN[mcurr+1] - valsIN[mcurr]);
  }

  return 0;
}


//! Return value of scalar quantity at level z (m) above base of ice (by interpolation).
PetscScalar     IceModelVec3::getValZ(const PetscInt i, const PetscInt j, const PetscScalar z) {
  // use linear interpolation
  checkHaveArray();
  isLegalLevel(z);
  if (z == Lz)
    return array[i][j][Mz - 1];
  const PetscInt     kbz = static_cast<PetscInt>(floor(z/dz));  // k value just below z
  const PetscScalar  val_kbz = array[i][j][kbz];
  return val_kbz + ( (z - levels[kbz]) / dz ) * (array[i][j][kbz + 1] - val_kbz);
}


//! Return values on planar star stencil of scalar quantity at level z.
PetscErrorCode   IceModelVec3::getPlaneStarZ(const PetscInt i, const PetscInt j, const PetscScalar z,
                                            planeStar *star) {
  // use linear interpolation
  CHKERRQ(checkHaveArray());
  CHKERRQ(isLegalLevel(z));
  // check ownership here?
  if (z == Lz) {
    star->ij = array[i][j][Mz - 1];
    star->ip1 = array[i+1][j][Mz - 1];
    star->im1 = array[i-1][j][Mz - 1];
    star->jp1 = array[i][j+1][Mz - 1];
    star->jm1 = array[i][j-1][Mz - 1];
    return 0;
  }
  const PetscInt     kbz = static_cast<PetscInt>(floor(z/dz));  // k value just below z
  const PetscScalar  incr = ( (z - levels[kbz]) / dz );
  star->ij = array[i][j][kbz] + incr * (array[i][j][kbz + 1] - array[i][j][kbz]);
  star->ip1 = array[i+1][j][kbz] + incr * (array[i+1][j][kbz + 1] - array[i+1][j][kbz]);
  star->im1 = array[i-1][j][kbz] + incr * (array[i-1][j][kbz + 1] - array[i-1][j][kbz]);
  star->jp1 = array[i][j+1][kbz] + incr * (array[i][j+1][kbz + 1] - array[i][j+1][kbz]);
  star->jm1 = array[i][j-1][kbz] + incr * (array[i][j-1][kbz + 1] - array[i][j-1][kbz]);
  return 0;
}


//! Return values of scalar quantity at given levels (m) above base of ice (by interpolation).
/*!
Input array \c levelsIN must be an allocated array of \c nlevels scalars (i.e. \c PetscScalar).
Likewise, return array \c vals must be an allocated array of \c nlevels scalars (i.e. \c PetscScalar).
Upon return, \c vals will be filled with values of scalar quantity at the \c levels.
 */
PetscErrorCode  IceModelVec3::getValColumn(
                     const PetscInt i, const PetscInt j, 
                     const PetscInt nlevels, PetscScalar *levelsIN,
                     PetscScalar *valsOUT) {
  
  for (PetscInt k=0;  k < nlevels; k++) {
    valsOUT[k] = getValZ(i,j,levelsIN[k]);
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
  
  *nlevels = Mz;
  *levelsOUT = levels;
  *valsOUT = array[i][j];
  return 0;
}




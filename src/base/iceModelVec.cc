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

  strcpy(varname,"*****UNKNOWN**** NetCDF varname_nc");
  strcpy(long_name,"UNKNOWN long_name");
  strcpy(units,"UNKNOWN units");
  strcpy(pism_intent,"UNKNOWN pism_intent");

  has_standard_name = PETSC_FALSE;
  strcpy(standard_name,"UNKNOWN standard_name");
  
  varid_nc = -9999;
};


IceModelVec::~IceModelVec() {

//  destroy();
};


PetscErrorCode  IceModelVec::create(IceGrid &mygrid, const char my_varname[]) {
  SETERRQ(1,"not implemented");
  return 0;
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


PetscErrorCode  IceModelVec::setAttrsNC(const int my_varid,
             const char my_long_name[], const char my_units[], const char my_pism_intent[]) {
  varid_nc = my_varid;
  strcpy(long_name,my_long_name);
  strcpy(units,my_units);
  strcpy(pism_intent,my_pism_intent);
  has_standard_name = PETSC_FALSE;
  return 0;
}


PetscErrorCode  IceModelVec::setAttrsCFstandardNC(const int my_varid,
             const char my_long_name[], const char my_units[], const char my_pism_intent[],
             const char my_standard_name[]) {
  varid_nc = my_varid;
  strcpy(long_name,my_long_name);
  strcpy(units,my_units);
  strcpy(pism_intent,my_pism_intent);
  strcpy(standard_name,my_standard_name);
  has_standard_name = PETSC_TRUE;
  return 0;
}


PetscErrorCode  IceModelVec::writeAttrsNC(const int ncid) {
  SETERRQ(1,"not YET implemented");
  return 0;
}


PetscErrorCode  IceModelVec::readVecNC(const int ncid, PetscTruth *exists) {
  SETERRQ(1,"not YET implemented");
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
  return 0;
}


PetscErrorCode IceModelVec::putVecNC(const int ncid, Vec g, const int *s, const int *c, int dims, 
                                          void *a_mpi, int a_size) {
  PetscErrorCode ierr;
  ierr = put_local_var(grid, ncid, varid_nc, NC_FLOAT, *da, v, g,
                         s, c, dims, a_mpi, a_size); CHKERRQ(ierr);  
  return 0;
}


PetscErrorCode  IceModelVec::needAccessToVals() {
  SETERRQ(1,"virtual; not implemented");
  return 0;
}


PetscErrorCode  IceModelVec::doneAccessToVals() {
  SETERRQ(1,"virtual; not implemented");
  return 0;
}


PetscErrorCode  IceModelVec::checkAllocated() {
  if (allocated == PETSC_FALSE) {
    SETERRQ1(1,"IceModelVec ERROR: IceModelVec with varname='%s' NOT allocated\n",
             varname);
  }
  return 0;
}


PetscErrorCode  IceModelVec::checkHaveArray() {
  SETERRQ(1,"virtual; not implemented");
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
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(*da, v, INSERT_VALUES, v);  CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec::endGhostComm() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(*da, v, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec::setToConstant(const PetscScalar c) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = VecSet(v,c); CHKERRQ(ierr);
  return 0;
}



/****************** 3d Vecs ********************/

IceModelVec3::IceModelVec3() : IceModelVec() {
  array = PETSC_NULL;
};


PetscErrorCode  IceModelVec3::destroy() {
  PetscErrorCode ierr = IceModelVec::destroy(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::create(IceGrid &mygrid, const char my_varname[]) {
  if (allocated == PETSC_TRUE) {
    SETERRQ1(1,"IceModelVec3 with varname='%s' already allocated\n",
             varname);
  }
  
  grid = &mygrid;
  da = &(mygrid.da3);
  
  PetscErrorCode ierr = DACreateLocalVector(*da, &v); CHKERRQ(ierr);
  allocated = PETSC_TRUE;

  strcpy(varname,my_varname);
  return 0;
}


PetscErrorCode  IceModelVec3::needAccessToVals() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DAVecGetArray(*da, v, &array); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::doneAccessToVals() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(*da, v, &array); CHKERRQ(ierr);
  array = PETSC_NULL;
  return 0;
}


PetscErrorCode  IceModelVec3::checkHaveArray() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  if (array == PETSC_NULL) {
    SETERRQ1(1,"array for IceModelVec3 with varname='%s' not available\n"
               "  (REMEMBER TO RUN needAccessToVals() before access and doneAccessToVals after access)\n",
               varname);
  }
  return 0;
}


PetscErrorCode  IceModelVec3::isLegalLevel(const PetscScalar z) {
  if (z < 0.0 - 1.0e-6) {
    SETERRQ2(1,"level z = %5.4f is below base of ice (z must be nonnegative);\n"
               "  IceModelVec3 has varname='%s'; ENDING!\n",
              z,varname);
  }
  if (z > (grid->p)->Lz + 1.0e-6) {
    SETERRQ3(2,"level z = %20.18f is above top of computational grid Lz = %20.18f;\n"
               "  IceModelVec3 has varname='%s'; ENDING!\n",
              z, (grid->p)->Lz,varname);
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
  PetscInt mcurr = 0;
  for (PetscInt k=0; k < (grid->p)->Mz; k++) {
    while (levelsIN[mcurr+1] < levels[k]) {
      mcurr++;
    }
    const PetscScalar increment = (levels[k] - levelsIN[mcurr]) / (levelsIN[mcurr+1] - levelsIN[mcurr]);
    array[i][j][k] = valsIN[mcurr] +  increment * (valsIN[mcurr+1] - valsIN[mcurr]);
  }

  return 0;
}


PetscErrorCode  IceModelVec3::setToConstantColumn(const PetscInt i, const PetscInt j, 
                                                  const PetscScalar c) {

  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  // check if in ownership ?

  for (PetscInt k=0; k < (grid->p)->Mz; k++) {
    array[i][j][k] = c;
  }

  return 0;
}


//! Return value of scalar quantity at level z (m) above base of ice (by interpolation).
PetscScalar     IceModelVec3::getValZ(const PetscInt i, const PetscInt j, const PetscScalar z) {
  // use linear interpolation
  checkHaveArray();
  isLegalLevel(z);
  if (z >= (grid->p)->Lz)
    return array[i][j][(grid->p)->Mz - 1];
  else if (z <= 0.0)
    return array[i][j][0];
  const PetscScalar  dz = (grid->p)->dz;
  const PetscInt     kbz = static_cast<PetscInt>(floor( z / dz ));  // k value just below z
  const PetscScalar  val_kbz = array[i][j][kbz];
  return val_kbz + ( (z - (grid->zlevels)[kbz]) / dz ) * (array[i][j][kbz + 1] - val_kbz);
}


//! Return values on planar star stencil of scalar quantity at level z.
PetscErrorCode   IceModelVec3::getPlaneStarZ(const PetscInt i, const PetscInt j, const PetscScalar z,
                                             planeStar *star) {
  PetscErrorCode ierr;
  // use linear interpolation
  ierr = checkHaveArray();  CHKERRQ(ierr);
  ierr = isLegalLevel(z);  CHKERRQ(ierr);
  // check ownership here?

  const PetscScalar dz = (grid->p)->dz;
  PetscInt     kbz;
  PetscScalar  incr;
  if (z >= (grid->p)->Lz) {
    kbz = (grid->p)->Mz - 1;
    incr = 0.0;
  } else if (z <= 0.0) {
    kbz = 0;
    incr = 0.0;
  } else {
    kbz = static_cast<PetscInt>(floor( z / dz ));  // k value just below z
    incr = ( (z - (grid->zlevels)[kbz]) / dz );
  }

  star->ij = array[i][j][kbz] + incr * (array[i][j][kbz + 1] - array[i][j][kbz]);
  star->ip1 = array[i+1][j][kbz] + incr * (array[i+1][j][kbz + 1] - array[i+1][j][kbz]);
  star->im1 = array[i-1][j][kbz] + incr * (array[i-1][j][kbz + 1] - array[i-1][j][kbz]);
  star->jp1 = array[i][j+1][kbz] + incr * (array[i][j+1][kbz + 1] - array[i][j+1][kbz]);
  star->jm1 = array[i][j-1][kbz] + incr * (array[i][j-1][kbz + 1] - array[i][j-1][kbz]);
  return 0;
}


//! Return values on planar box stencil of scalar quantity at level z.
PetscErrorCode   IceModelVec3::getPlaneBoxZ(const PetscInt i, const PetscInt j, const PetscScalar z,
                                            planeBox *box) {
  // use linear interpolation
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = isLegalLevel(z); CHKERRQ(ierr);
  // check ownership here?

  const PetscScalar dz = (grid->p)->dz;
  PetscInt     kbz;
  PetscScalar  incr;
  if (z >= (grid->p)->Lz) {
    kbz = (grid->p)->Mz - 1;
    incr = 0.0;
  } else if (z <= 0.0) {
    kbz = 0;
    incr = 0.0;
  } else {
    kbz = static_cast<PetscInt>(floor( z / dz ));  // k value just below z
    incr = ( (z - (grid->zlevels)[kbz]) / dz );
  }

  box->ij = array[i][j][kbz] + incr * (array[i][j][kbz + 1] - array[i][j][kbz]);
  box->ip1 = array[i+1][j][kbz] + incr * (array[i+1][j][kbz + 1] - array[i+1][j][kbz]);
  box->im1 = array[i-1][j][kbz] + incr * (array[i-1][j][kbz + 1] - array[i-1][j][kbz]);
  box->jp1 = array[i][j+1][kbz] + incr * (array[i][j+1][kbz + 1] - array[i][j+1][kbz]);
  box->jm1 = array[i][j-1][kbz] + incr * (array[i][j-1][kbz + 1] - array[i][j-1][kbz]);
  box->ip1jp1 = array[i+1][j+1][kbz] + incr * (array[i+1][j+1][kbz + 1] - array[i+1][j+1][kbz]);
  box->im1jp1 = array[i-1][j+1][kbz] + incr * (array[i-1][j+1][kbz + 1] - array[i-1][j+1][kbz]);
  box->ip1jm1 = array[i+1][j-1][kbz] + incr * (array[i+1][j-1][kbz + 1] - array[i+1][j-1][kbz]);
  box->im1jm1 = array[i-1][j-1][kbz] + incr * (array[i-1][j-1][kbz + 1] - array[i-1][j-1][kbz]);

  return 0;
}


//! Return values of scalar quantity at given levels (m) above base of ice (by interpolation).
/*!
Input array \c levelsIN must be an allocated array of \c nlevels scalars (i.e. \c PetscScalar).

\c levelsIN must be strictly increasing and in the range \f$0 < z < \mathtt{grid.p->Lz}\f$.

Likewise, return array \c valsOUT must be an allocated array of \c nlevels scalars (i.e. \c PetscScalar).
Upon return, \c valsOUT will be filled with values of scalar quantity at the \f$z\f$ values in \c levelsIN.
 */
PetscErrorCode  IceModelVec3::getValColumn(
                     const PetscInt i, const PetscInt j, 
                     const PetscInt nlevelsIN, PetscScalar *levelsIN,
                     PetscScalar *valsOUT) {
  
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
  
  PetscInt mcurr = 0;
  for (PetscInt k = 0; k < nlevelsIN; k++) {
    while (levels[mcurr+1] < levelsIN[k]) {
      mcurr++;
    }
    const PetscScalar incr = (levelsIN[k] - levels[mcurr]) / (levels[mcurr+1] - levels[mcurr]);
    const PetscScalar valm = array[i][j][mcurr];
    valsOUT[k] = valm + incr * (array[i][j][mcurr+1] - valm);
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
  *valsOUT = array[i][j];
  return 0;
}


// Copyright (C) 2008, 2009 Ed Bueler and Constantine Khroulev
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

// methods for base class IceModelVec and derived class IceModelVec2
// are in "iceModelVec.cc"

IceModelVec3::IceModelVec3() : IceModelVec() {}


//! Allocate a DA and a Vec from information in IceGrid.
PetscErrorCode  IceModelVec3::create(IceGrid &mygrid, const char my_short_name[], bool local) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec3 with short_name='%s' already allocated\n",short_name);
  }
  
  grid = &mygrid;
  dims = GRID_3D;

  PetscInt       M, N, m, n;
  PetscErrorCode ierr;
  ierr = DAGetInfo(mygrid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate3d(mygrid.com, DA_YZPERIODIC, DA_STENCIL_STAR, mygrid.Mz, N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da); CHKERRQ(ierr);
  IOwnDA = true;

  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  strcpy(short_name,my_short_name);
#ifdef PISM_DEBUG
  creation_counter += 1;
#endif // PISM_DEBUG
  return 0;
}

//! Allocate a DA and a Vec from information in IceGrid; use an existing DA from an existing IceModelVec3.
PetscErrorCode  IceModelVec3::createSameDA(IceModelVec3 imv3_source,
                                           IceGrid &mygrid, const char my_short_name[], bool local) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec3 with short_name='%s' already allocated\n",short_name);
  }
  
  grid = &mygrid;
  dims = GRID_3D;

  da = imv3_source.da;
  IOwnDA = false;

  PetscErrorCode ierr;
  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  strcpy(short_name,my_short_name);
#ifdef PISM_DEBUG
  creation_counter += 1;
#endif // PISM_DEBUG
  return 0;
}

//! Defines a netcdf variable corresponding to an IceModelVec3 object. The ncid
// argument must refer to a dataset with dimensions t, x, y, z.
PetscErrorCode IceModelVec3::define_netcdf_variable(int ncid, nc_type nctype, int *varidp) {
  int stat, dimids[4], var_id;

  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "t", &dimids[0]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "y", &dimids[1]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "x", &dimids[2]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "z", &dimids[3]); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_def_var(ncid, short_name, nctype, 4, dimids, &var_id);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  stat = MPI_Bcast(&var_id, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);

  *varidp = var_id;

  return 0;
}


PetscErrorCode  IceModelVec3::beginGhostCommTransfer(IceModelVec3 imv3_source) {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec3!\n"
               "  (has short_name='%s')\n", short_name);
  }
  if (imv3_source.localp) {
    SETERRQ1(2,"source IceModelVec3 must be GLOBAL! (has short_name='%s')\n",
               imv3_source.short_name);
  }
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = imv3_source.checkAllocated(); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, imv3_source.v, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::endGhostCommTransfer(IceModelVec3 imv3_source) {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec3!\n"
               "  (has short_name='%s')\n",
               short_name);
  }
  if (imv3_source.localp) {
    SETERRQ1(2,"source IceModelVec3 must be GLOBAL! (has short_name='%s')\n",
               imv3_source.short_name);
  }
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = imv3_source.checkAllocated(); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, imv3_source.v, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::isLegalLevel(const PetscScalar z) {
  if (z < 0.0 - 1.0e-6) {
    SETERRQ2(1,"level z = %5.4f is below base of ice (z must be nonnegative);\n"
               "  IceModelVec3 has short_name='%s'; ENDING!\n",
              z,short_name);
  }
  if (z > grid->Lz + 1.0e-6) {
    SETERRQ3(2,"level z = %10.8f is above top of computational grid Lz = %10.8f;\n"
               "  IceModelVec3 has short_name='%s'; ENDING!\n",
              z, grid->Lz,short_name);
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
PetscErrorCode  IceModelVec3::setValColumnPL(
                     const PetscInt i, const PetscInt j, 
                     const PetscInt nlevels, PetscScalar *levelsIN,
                     PetscScalar *valsIN) {

  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // check if in ownership ?

  if (levelsIN[0] > 0.0 + 1.0e-3) {
    SETERRQ2(1,"levelsIN[0]=%10.9f is above base of ice at z=0 so *interpolation*\n"
              "   is impossible; IceModelVec3 has short_name='%s';  ENDING!\n",
              levelsIN[0],short_name);
  }
  if (levelsIN[nlevels - 1] < grid->Lz - 1.0e-3) {
    SETERRQ3(2,"levelsIN[nlevels-1] = %10.9f is below top of computational domain\n"
               "   at z=Lz=%10.9f, so *interpolation* is impossible;\n"
               "   IceModelVec3 has short_name='%s';  ENDING!\n",
               levelsIN[nlevels-1],grid->Lz,short_name);
  }
  for (PetscInt k=0; k < nlevels - 1; k++) {
    if (levelsIN[k] >= levelsIN[k+1]) {
      SETERRQ2(3,"levelsIN not *strictly increasing* at index %d;\n"
                 "    IceModelVec3 has short_name='%s';  ENDING!\n",
                 k,short_name);
    }
  }

  PetscScalar *levels;
  levels = grid->zlevels;
  PetscScalar ***arr = (PetscScalar***) array;
  
  PetscInt mcurr = 0;
  for (PetscInt k=0; k < grid->Mz; k++) {
    while (levelsIN[mcurr+1] < levels[k]) {
      mcurr++;
    }
    const PetscScalar increment = (levels[k] - levelsIN[mcurr])
                                  / (levelsIN[mcurr+1] - levelsIN[mcurr]);
    arr[i][j][k] = valsIN[mcurr] +  increment * (valsIN[mcurr+1] - valsIN[mcurr]);
  }

  return 0;
}


//! Set all values of scalar quantity to given a single value in a particular column.
PetscErrorCode  IceModelVec3::setColumn(
                   const PetscInt i, const PetscInt j, const PetscScalar c) {

  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  // check if in ownership ?
  PetscScalar ***arr = (PetscScalar***) array;
  for (PetscInt k=0; k < grid->Mz; k++) {
    arr[i][j][k] = c;
  }
  return 0;
}


//! Return value of scalar quantity at level z (m) above base of ice (by linear interpolation).
PetscScalar IceModelVec3::getValZ(const PetscInt i, const PetscInt j, 
                                  const PetscScalar z) {

  if (checkHaveArray() != 0) {
    PetscPrintf(PETSC_COMM_SELF, 
       "IceModelVec3 getValZ(): array was not allocated (so says\n"
       "  IceModelVec::checkHaveArray()); short_name = %s\n", short_name);
    PetscEnd();
  }
  if (isLegalLevel(z) != 0) {
    PetscPrintf(PETSC_COMM_SELF, 
       "IceModelVec3 getValZ(): isLegalLevel() says level %f was\n"
       "  not legal; short_name = %s\n", z, short_name);
    PetscEnd();
  }
  PetscScalar ***arr = (PetscScalar***) array;
  if (z >= grid->Lz)
    return arr[i][j][grid->Mz - 1];
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


//! Return values on planar star stencil of scalar quantity at level z (by linear interpolation).
PetscErrorCode   IceModelVec3::getPlaneStarZ(
                     const PetscInt i, const PetscInt j, const PetscScalar z,
                     planeStar *star) {
  PetscErrorCode ierr;
  ierr = checkHaveArray();  CHKERRQ(ierr);
  ierr = isLegalLevel(z);  CHKERRQ(ierr);
  // check ownership here?
  if (!localp) {
    SETERRQ1(1,"IceModelVec3 ERROR: IceModelVec3 with short_name='%s' is GLOBAL\n"
               "  and cannot do getPlaneStarZ()\n", short_name);
  }

  PetscInt     kbz;
  PetscScalar  incr;
  if (z >= grid->Lz) {
    kbz = grid->Mz - 1;
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
Input array \c levelsIN must be an allocated array of \c nlevels scalars 
(\c PetscScalar).

\c levelsIN must be strictly increasing and in the range 
\f$0 <= z <= \mathtt{grid.Lz}\f$.

Return array \c valsOUT must be an allocated array of \c nlevels scalars 
(\c PetscScalar).

Upon return, \c valsOUT will be filled with values of scalar quantity at 
the \f$z\f$ values in \c levelsIN.
 */
PetscErrorCode  IceModelVec3::getValColumnPL(const PetscInt i, const PetscInt j, 
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
                 "  (IceModelVec3 with short_name='%s')  ENDING!\n",k,short_name);
    }
  }

  PetscScalar* levels = grid->zlevels;
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

//! Return values of ice scalar quantity at given levels (m) above base of ice, using local quadratic interpolation.
/*!
Input array \c levelsIN must be an allocated array of \c nlevels scalars 
(\c PetscScalar).

\c levelsIN must be strictly increasing and in the range 
\f$0 <= z <= \mathtt{grid.Lz}\f$.

Return array \c valsOUT must be an allocated array of \c nlevels scalars 
(\c PetscScalar).

Upon return, \c valsOUT will be filled with values of scalar quantity 
at the \f$z\f$ values in \c levelsIN.
 */
PetscErrorCode  IceModelVec3::getValColumnQUAD(const PetscInt i, const PetscInt j, 
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
                 "    (IceModelVec3 with short_name='%s')  ENDING!\n",k,short_name);
    }
  }

  PetscScalar* levels = grid->zlevels;
  PetscScalar ***arr = (PetscScalar***) array;
  
  PetscInt mcurr = 0;
  for (PetscInt k = 0; k < nlevelsIN; k++) {
    while (levels[mcurr+1] < levelsIN[k]) {
      mcurr++;
    }
    const PetscScalar z0 = levels[mcurr],
                      f0 = arr[i][j][mcurr];
    if (mcurr >= grid->Mz - 2) {
      // just do linear interpolation at top of grid
      const PetscScalar incr = (levelsIN[k] - z0) / (levels[mcurr+1] - z0);
      valsOUT[k] = f0 + incr * (arr[i][j][mcurr+1] - f0);
    } else {
      const PetscScalar dz1 = levels[mcurr+1] - z0,
                        dz2 = levels[mcurr+2] - z0;
      const PetscScalar D1 = (arr[i][j][mcurr+1] - f0) / dz1,
                        D2 = (arr[i][j][mcurr+2] - f0) / dz2;
      const PetscScalar c = (D2 - D1) / (dz2 - dz1),
                        b = D1 - c * dz1;
      const PetscScalar s = levelsIN[k] - z0;
      valsOUT[k] = f0 + s * (b + c * s);
    }
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

PetscErrorCode  IceModelVec3::getHorSlice(IceModelVec2 gslice, const PetscScalar z) {
  PetscErrorCode ierr;
  PetscScalar    **slice_val;
  ierr = gslice.get_array(slice_val); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      slice_val[i][j] = getValZ(i,j,z);
    }
  }
  ierr = gslice.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::getSurfaceValues(IceModelVec2 &gsurf, IceModelVec2 myH) {
  PetscErrorCode ierr;
  PetscScalar    **H;
  ierr = myH.get_array(H); CHKERRQ(ierr);
  ierr = getSurfaceValues(gsurf, H); CHKERRQ(ierr);
  ierr = myH.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode  IceModelVec3::getSurfaceValues(Vec &gsurf, IceModelVec2 myH) {
  PetscErrorCode ierr;
  PetscScalar    **H, **surf_val;
  ierr = DAVecGetArray(da, gsurf, &surf_val); CHKERRQ(ierr);
  ierr = myH.get_array(H); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      surf_val[i][j] = getValZ(i,j,H[i][j]);
    }
  }
  ierr = myH.end_access(); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da, gsurf, &surf_val); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::getSurfaceValues(IceModelVec2 &gsurf, PetscScalar **H) {
  PetscErrorCode ierr;
  PetscScalar    **surf_val;
  ierr = gsurf.get_array(surf_val); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      surf_val[i][j] = getValZ(i,j,H[i][j]);
    }
  }
  ierr = gsurf.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::getInternalColumn(const PetscInt i, const PetscInt j, PetscScalar **valsPTR) {
  
  PetscScalar ***arr = (PetscScalar***) array;
  *valsPTR = arr[i][j];
  return 0;
}


PetscErrorCode  IceModelVec3::setInternalColumn(const PetscInt i, const PetscInt j, PetscScalar *valsIN) {
  
  PetscScalar ***arr = (PetscScalar***) array;
  for (PetscInt k = 0; k < grid->Mz; k++) {
    arr[i][j][k] = valsIN[k];
  }
  return 0;
}

/********* IceModelVec3Bedrock **********/

IceModelVec3Bedrock::IceModelVec3Bedrock() : IceModelVec() {}


//! Allocate a DA and a Vec from information in IceGrid.
PetscErrorCode  IceModelVec3Bedrock::create(IceGrid &my_grid, 
                               const char my_short_name[], bool local) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }
  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec3Bedrock with short_name='%s' already allocated\n",short_name);
  }
  if (local) {
    SETERRQ1(2,"IceModelVec3Bedrock must be GLOBAL (short_name='%s')\n",short_name);
  }

  strcpy(short_name,my_short_name);

  grid = &my_grid;
  dims = GRID_3D_BEDROCK;
  
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
#ifdef PISM_DEBUG
  creation_counter += 1;
#endif // PISM_DEBUG
  return 0;
}

//! Defines a netcdf variable corresponding to an IceModelVec3Bedrock object. The ncid
// argument must refer to a dataset with dimensions t, x, y, zb.
PetscErrorCode IceModelVec3Bedrock::define_netcdf_variable(int ncid, nc_type nctype, int *varidp) {
  int stat, dimids[4], var_id;

  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "t", &dimids[0]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "y", &dimids[1]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "x", &dimids[2]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "zb", &dimids[3]); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_def_var(ncid, short_name, nctype, 4, dimids, &var_id);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  stat = MPI_Bcast(&var_id, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);

  *varidp = var_id;

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
PetscErrorCode  IceModelVec3Bedrock::setColumn(
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

//   for (PetscInt k=0; k < nlevels; k++)
//     PetscPrintf(grid->com, "levels[%d] = %10.3f\n", k, levelsIN[k]);
//   PetscPrintf(grid->com, "\n");

  if (levelsIN[0] > -grid->Lbz + 1.0e-3) {
    SETERRQ3(1,"levelsIN[0]=%10.9f is above base of bedrock at z=-%10.9f so *interpolation*\n"
              "   is impossible; IceModelVec3Bedrock has short_name='%s';  ENDING!\n",
              levelsIN[0],grid->Lbz,short_name);
  }
  if (levelsIN[nlevels - 1] < 0.0 - 1.0e-3) {
    SETERRQ2(2,"levelsIN[nlevels-1] = %10.9f is below z=0, so *interpolation* is impossible;\n"
               "   IceModelVec3Bedrock has short_name='%s';  ENDING!\n",
               levelsIN[nlevels-1],short_name);
  }
  for (PetscInt k=0; k < nlevels - 1; k++) {
    if (levelsIN[k] >= levelsIN[k+1]) {
      SETERRQ2(3,"levelsIN not *strictly increasing* at index %d;\n"
                 "    IceModelVec3Bedrock has short_name='%s';  ENDING!\n",
                 k,short_name);
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
       "level z = %10.8f is below bottom of bedrock at -Lbz = %10.8f; IceModelVec3Bedrock has short_name='%s'; ENDING!\n",
       z,-grid->Lbz,short_name);
  }
  if (z > 0.0 + 1.0e-6) {
    SETERRQ2(2,"level z = %10.8f is above top of bedrock at z=0; IceModelVec3Bedrock has short_name='%s'; ENDING!\n",
              z,short_name);
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
                 "    (IceModelVec3Bedrock with short_name='%s')  ENDING!\n",k,short_name);
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

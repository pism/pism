// Copyright (C) 2008--2010 Ed Bueler and Constantine Khroulev
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

#include <sstream>
#include <cstring>
#include <cstdlib>
#include <petscda.h>
#include <netcdf.h>
#include "nc_util.hh"

#include "iceModelVec.hh"

// this file contains method for derived class IceModelVec3

// methods for base class IceModelVec and derived class IceModelVec2S
// are in "iceModelVec.cc"

IceModelVec3::IceModelVec3() : IceModelVec() {
  sounding_buffer = PETSC_NULL;
  slice_viewers = new map<string, PetscViewer>;
  sounding_viewers = new map<string, PetscViewer>;
}

IceModelVec3::~IceModelVec3() {
  if (!shallow_copy) {
    delete slice_viewers;
    delete sounding_viewers;
  }
}

IceModelVec3::IceModelVec3(const IceModelVec3 &other)
  : IceModelVec(other) {
  slice_viewers = other.slice_viewers;
  sounding_buffer = other.sounding_buffer;
  sounding_viewers = other.sounding_viewers;
  shallow_copy = true;
}

PetscErrorCode IceModelVec3::create_da(DA &result, PetscInt Mz) {
  PetscErrorCode ierr;
  PetscInt       M, N, m, n;
  const PetscInt *lx, *ly;

  ierr = DAGetInfo(grid->da2,
		   PETSC_NULL,	       // dim
		   &N, &M, PETSC_NULL, // N, M, P
		   &n, &m, PETSC_NULL, // n, m, p
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

  ierr = DAGetOwnershipRanges(grid->da2, &ly, &lx, PETSC_NULL); CHKERRQ(ierr);

  ierr = DACreate3d(grid->com, DA_YZPERIODIC, DA_STENCIL_STAR,
		    Mz, N, M,	// P, N, M
		    1, n, m,	// p, n, m
		    1, s_width,	// dof, stencil_width
                    PETSC_NULL, ly, lx, // lz, ly, lx
		    &result); CHKERRQ(ierr);

  return 0;
}

//! Allocate a DA and a Vec from information in IceGrid.
PetscErrorCode  IceModelVec3::create(IceGrid &my_grid, const char my_name[],
				     bool local, int stencil_width) {
  PetscErrorCode ierr;
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec3 with name='%s' already allocated\n",name.c_str());
  }
  
  grid = &my_grid;
  dims = GRID_3D;
  s_width = stencil_width;

  ierr = create_da(da, grid->Mz); CHKERRQ(ierr);

  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  name = my_name;

  vars[0].init(my_name, my_grid, GRID_3D);

  //  ierr = this->set(GSL_NAN); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec3::destroy() {
  PetscErrorCode ierr;
  map<string,PetscViewer>::iterator i;

  ierr = IceModelVec::destroy(); CHKERRQ(ierr);

  // soundings:
  if (sounding_viewers != NULL) {
    for (i = (*sounding_viewers).begin(); i != (*sounding_viewers).end(); ++i) {
      if ((*i).second != PETSC_NULL) {
	ierr = PetscViewerDestroy((*i).second); CHKERRQ(ierr);
      }
    }
    delete sounding_viewers;
    sounding_viewers = NULL;
  }

  if (sounding_buffer != PETSC_NULL) {
    ierr = VecDestroy(sounding_buffer); CHKERRQ(ierr);
    sounding_buffer = PETSC_NULL;
  }

  // slices:
  if (slice_viewers != NULL) {
    for (i = (*slice_viewers).begin(); i != (*slice_viewers).end(); ++i) {
      if ((*i).second != PETSC_NULL) {
	ierr = PetscViewerDestroy((*i).second); CHKERRQ(ierr);
      }
    }
    delete slice_viewers;
    slice_viewers = NULL;
  }

  return 0;
}

PetscErrorCode  IceModelVec3::beginGhostCommTransfer(IceModelVec3 &imv3_source) {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec3!\n"
               "  (has name='%s')\n", name.c_str());
  }
  if (imv3_source.localp) {
    SETERRQ1(2,"source IceModelVec3 must be GLOBAL! (has name='%s')\n",
               imv3_source.name.c_str());
  }
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = imv3_source.checkAllocated(); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da, imv3_source.v, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::endGhostCommTransfer(IceModelVec3 &imv3_source) {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec3!\n"
               "  (has name='%s')\n",
               name.c_str());
  }
  if (imv3_source.localp) {
    SETERRQ1(2,"source IceModelVec3 must be GLOBAL! (has name='%s')\n",
               imv3_source.name.c_str());
  }
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = imv3_source.checkAllocated(); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da, imv3_source.v, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::isLegalLevel(PetscScalar z) {
  if (z < 0.0 - 1.0e-6) {
    SETERRQ2(1,"level z = %5.4f is below base of ice (z must be nonnegative);\n"
               "  IceModelVec3 has name='%s'; ENDING!\n",
              z,name.c_str());
  }
  if (z > grid->Lz + 1.0e-6) {
    SETERRQ3(2,"level z = %10.8f is above top of computational grid Lz = %10.8f;\n"
               "  IceModelVec3 has name='%s'; ENDING!\n",
              z, grid->Lz,name.c_str());
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
PetscErrorCode  IceModelVec3::setValColumnPL(PetscInt i, PetscInt j, PetscScalar *source) {
#ifdef PISM_DEBUG
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // check if in ownership ?
#endif

  const PetscScalar *levels = grid->zlevels;
  const PetscScalar *levels_fine = grid->zlevels_fine;
  PetscScalar ***arr = (PetscScalar***) array;
  
  for (PetscInt k=0; k < grid->Mz; ++k) {
    PetscInt m = grid->ice_fine2storage[k];

    const PetscScalar increment = (levels[k] - levels_fine[m])
                                  / (levels_fine[m+1] - levels_fine[m]);
    arr[i][j][k] = source[m] +  increment * (source[m+1] - source[m]);
  }

  return 0;
}


//! Set all values of scalar quantity to given a single value in a particular column.
PetscErrorCode  IceModelVec3::setColumn(PetscInt i, PetscInt j, PetscScalar c) {
#ifdef PISM_DEBUG
  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  // check if in ownership ?
#endif
  PetscScalar ***arr = (PetscScalar***) array;
  for (PetscInt k=0; k < grid->Mz; k++) {
    arr[i][j][k] = c;
  }
  return 0;
}


//! Return value of scalar quantity at level z (m) above base of ice (by linear interpolation).
PetscScalar IceModelVec3::getValZ(PetscInt i, PetscInt j, PetscScalar z) {

  if (checkHaveArray() != 0) {
    PetscPrintf(PETSC_COMM_SELF, 
       "IceModelVec3 getValZ(): array was not allocated (so says\n"
       "  IceModelVec::checkHaveArray()); name = %s\n", name.c_str());
    PetscEnd();
  }
  if (isLegalLevel(z) != 0) {
    PetscPrintf(PETSC_COMM_SELF, 
       "IceModelVec3 getValZ(): isLegalLevel() says level %f was\n"
       "  not legal; name = %s\n", z, name.c_str());
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
PetscErrorCode   IceModelVec3::getPlaneStarZ(PetscInt i, PetscInt j, PetscScalar z,
					     planeStar *star) {
#ifdef PISM_DEBUG
  PetscErrorCode ierr;
  ierr = checkHaveArray();  CHKERRQ(ierr);
  ierr = isLegalLevel(z);  CHKERRQ(ierr);
  // check ownership here?
  if (!localp) {
    SETERRQ1(1,"IceModelVec3 ERROR: IceModelVec3 with name='%s' is GLOBAL\n"
               "  and cannot do getPlaneStarZ()\n", name.c_str());
  }
#endif

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

  if (kbz < grid->Mz - 1) {
    star->ij  = arr[i][j][kbz]   + incr * (arr[i][j][kbz + 1]   - arr[i][j][kbz]);
    star->ip1 = arr[i+1][j][kbz] + incr * (arr[i+1][j][kbz + 1] - arr[i+1][j][kbz]);
    star->im1 = arr[i-1][j][kbz] + incr * (arr[i-1][j][kbz + 1] - arr[i-1][j][kbz]);
    star->jp1 = arr[i][j+1][kbz] + incr * (arr[i][j+1][kbz + 1] - arr[i][j+1][kbz]);
    star->jm1 = arr[i][j-1][kbz] + incr * (arr[i][j-1][kbz + 1] - arr[i][j-1][kbz]);
  } else {
    star->ij  = arr[i][j][kbz];
    star->ip1 = arr[i+1][j][kbz];
    star->im1 = arr[i-1][j][kbz];
    star->jp1 = arr[i][j+1][kbz];
    star->jm1 = arr[i][j-1][kbz];
  }

  return 0;
}

//! Gets a map-plane star stencil directly from the storage grid.
PetscErrorCode IceModelVec3::getPlaneStar(PetscInt i, PetscInt j, PetscInt k,
						  planeStar *star) {
  PetscScalar ***arr = (PetscScalar***) array;

  star->ij  = arr[i][j][k];
  star->ip1 = arr[i+1][j][k];
  star->im1 = arr[i-1][j][k];
  star->jp1 = arr[i][j+1][k];
  star->jm1 = arr[i][j-1][k];

  return 0;
}

//! Gets a map-plane star stencil on the fine vertical grid.
PetscErrorCode IceModelVec3::getPlaneStar_fine(PetscInt i, PetscInt j, PetscInt k,
					       planeStar *star) {
  PetscInt kbz = grid->ice_storage2fine[k];

  if (kbz < grid->Mz - 1) {
    PetscScalar* levels = grid->zlevels;
    PetscScalar z = grid->zlevels_fine[k],
      incr = (z - levels[kbz]) / (levels[kbz+1] - levels[kbz]);
    PetscScalar ***arr = (PetscScalar***) array;

    star->ij  = arr[i][j][kbz]   + incr * (arr[i][j][kbz + 1]   - arr[i][j][kbz]);
    star->ip1 = arr[i+1][j][kbz] + incr * (arr[i+1][j][kbz + 1] - arr[i+1][j][kbz]);
    star->im1 = arr[i-1][j][kbz] + incr * (arr[i-1][j][kbz + 1] - arr[i-1][j][kbz]);
    star->jp1 = arr[i][j+1][kbz] + incr * (arr[i][j+1][kbz + 1] - arr[i][j+1][kbz]);
    star->jm1 = arr[i][j-1][kbz] + incr * (arr[i][j-1][kbz + 1] - arr[i][j-1][kbz]);
  } else {
    return getPlaneStar(i, j, kbz, star);
  }
  return 0;
}

//! Return values of ice scalar quantity at given levels (m) above base of ice, using piecewise linear interpolation.
/*!
Input array \c levelsIN must be an allocated array of \c nlevels scalars 
(\c PetscScalar).

\c levelsIN must be strictly increasing and positive. Exceeding Lz is allowed,
extrapolation (by the value at the top-most level) is performed in this case.

Return array \c valsOUT must be an allocated array of \c nlevels scalars 
(\c PetscScalar).

Upon return, \c valsOUT will be filled with values of scalar quantity at 
the \f$z\f$ values in \c levelsIN.
 */
PetscErrorCode IceModelVec3::getValColumnPL(PetscInt i, PetscInt j, PetscInt ks,
					    PetscScalar *result) {
  
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // check if in ownership ?

  PetscScalar *levels = grid->zlevels;
  PetscScalar *levels_fine = grid->zlevels_fine;
  PetscScalar ***arr = (PetscScalar***) array;

  for (PetscInt k = 0; k < grid->Mz_fine; k++) {
    if (k > ks) {
      result[k] = arr[i][j][grid->ice_storage2fine[k]];
      continue;
    }

    PetscInt m = grid->ice_storage2fine[k];

    // extrapolate (if necessary):
    if (m == grid->Mz - 1) {
      result[k] = arr[i][j][grid->Mz-1];
      continue;
    }

    const PetscScalar incr = (levels_fine[k] - levels[m]) / (levels[m+1] - levels[m]);
    const PetscScalar valm = arr[i][j][m];
    result[k] = valm + incr * (arr[i][j][m+1] - valm);
  }

  return 0;
}

//! Return values of ice scalar quantity on the fine computational grid, using local quadratic interpolation.
/*!

Return array \c valsOUT must be an allocated array of \c grid.Mz_fine scalars 
(\c PetscScalar).
 */
PetscErrorCode  IceModelVec3::getValColumnQUAD(PetscInt i, PetscInt j, PetscInt ks,
					       PetscScalar *result) {
  
  const PetscScalar *levels = grid->zlevels;
  const PetscScalar *levels_fine = grid->zlevels_fine;
  const PetscScalar ***arr = (const PetscScalar***) array;
  const PetscScalar *column = arr[i][j];

  for (PetscInt k = 0; k < grid->Mz_fine; k++) {
    if (k > ks) {
      result[k] = arr[i][j][grid->ice_storage2fine[k]];
      continue;
    }

    const PetscInt m = grid->ice_storage2fine[k];

    // extrapolate (if necessary):
    if (m == grid->Mz - 1) {
      result[k] = column[grid->Mz-1];
      continue;
    }

    const PetscScalar z0 = levels[m],
                      f0 = column[m];
    if (m >= grid->Mz - 2) {
				// top of the grid: just do linear interpolation
      const PetscScalar incr = (levels_fine[k] - z0) / (levels[m+1] - z0);
      result[k] = f0 + incr * (column[m+1] - f0);
    } else {			// the rest: one-sided quadratic interpolation
      const PetscScalar dz1 = levels[m+1] - z0,
                        dz2 = levels[m+2] - z0;
      const PetscScalar D1 = (column[m+1] - f0) / dz1,
                        D2 = (column[m+2] - f0) / dz2;
      const PetscScalar c = (D2 - D1) / (dz2 - dz1),
                        b = D1 - c * dz1;
      const PetscScalar s = levels_fine[k] - z0;
      result[k] = f0 + s * (b + c * s);
    }
  }

  return 0;
}


//! If the grid is equally spaced in the ice then use PL, otherwise use QUAD.
PetscErrorCode  IceModelVec3::getValColumn(PetscInt i, PetscInt j, PetscInt ks,
					   PetscScalar *result) {
  if (grid->ice_vertical_spacing == EQUAL) {
    return getValColumnPL(i, j, ks, result);
  } else {
    return getValColumnQUAD(i, j, ks, result);
  }
}


//! Copies a horizontal slice at level z of an IceModelVec3 into a Vec gslice.
PetscErrorCode  IceModelVec3::getHorSlice(Vec &gslice, PetscScalar z) {
  PetscErrorCode ierr;
  PetscScalar    **slice_val;

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid->da2, gslice, &slice_val); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      slice_val[i][j] = getValZ(i,j,z);
    }
  }
  ierr = DAVecRestoreArray(grid->da2, gslice, &slice_val); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! Copies a horizontal slice at level z of an IceModelVec3 into an IceModelVec2S gslice.
PetscErrorCode  IceModelVec3::getHorSlice(IceModelVec2S &gslice, PetscScalar z) {
  PetscErrorCode ierr;
  PetscScalar    **slice_val;

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = gslice.get_array(slice_val); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      slice_val[i][j] = getValZ(i,j,z);
    }
  }
  ierr = gslice.end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}


//! Copies the values of an IceModelVec3 at the ice surface (specified by the level myH) to an IceModelVec2S gsurf.
PetscErrorCode  IceModelVec3::getSurfaceValues(IceModelVec2S &gsurf, IceModelVec2S &myH) {
  PetscErrorCode ierr;
  PetscScalar    **H;
  ierr = myH.get_array(H); CHKERRQ(ierr);
  ierr = getSurfaceValues(gsurf, H); CHKERRQ(ierr);
  ierr = myH.end_access(); CHKERRQ(ierr);
  return 0;
}

//! Copies the values of an IceModelVec3 at the ice surface (specified by the level myH) to a Vec gsurf.
/*!
  This version is used in iMviewers.cc
 */
PetscErrorCode  IceModelVec3::getSurfaceValues(Vec &gsurf, IceModelVec2S &myH) {
  PetscErrorCode ierr;
  PetscScalar    **H, **surf_val;
  ierr = begin_access(); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid->da2, gsurf, &surf_val); CHKERRQ(ierr);
  ierr = myH.get_array(H); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      surf_val[i][j] = getValZ(i,j,H[i][j]);
    }
  }
  ierr = myH.end_access(); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid->da2, gsurf, &surf_val); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::getSurfaceValues(IceModelVec2S &gsurf, PetscScalar **H) {
  PetscErrorCode ierr;
  PetscScalar    **surf_val;

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = gsurf.get_array(surf_val); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      surf_val[i][j] = getValZ(i,j,H[i][j]);
    }
  }
  ierr = gsurf.end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode  IceModelVec3::getInternalColumn(PetscInt i, PetscInt j, PetscScalar **valsPTR) {
  
  PetscScalar ***arr = (PetscScalar***) array;
  *valsPTR = arr[i][j];
  return 0;
}


PetscErrorCode  IceModelVec3::setInternalColumn(PetscInt i, PetscInt j, PetscScalar *valsIN) {
  
  PetscScalar ***arr = (PetscScalar***) array;
  PetscErrorCode ierr = PetscMemcpy(arr[i][j], valsIN, grid->Mz*sizeof(PetscScalar));
  CHKERRQ(ierr);
  return 0;
}

/********* IceModelVec3Bedrock **********/

IceModelVec3Bedrock::IceModelVec3Bedrock() : IceModelVec() {
  sounding_buffer = PETSC_NULL;
  sounding_viewers = new map<string, PetscViewer>;
}

IceModelVec3Bedrock::~IceModelVec3Bedrock() {
  if (!shallow_copy) {
    delete sounding_viewers;
  }
}

IceModelVec3Bedrock::IceModelVec3Bedrock(const IceModelVec3Bedrock &other) : IceModelVec() {
  sounding_buffer = other.sounding_buffer;
  sounding_viewers = other.sounding_viewers;
  shallow_copy = true;
}


//! Allocate a DA and a Vec from information in IceGrid.
PetscErrorCode  IceModelVec3Bedrock::create(IceGrid &my_grid, 
                               const char my_name[], bool local) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }
  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec3Bedrock with name='%s' already allocated\n",name.c_str());
  }
  if (local) {
    SETERRQ1(2,"IceModelVec3Bedrock must be GLOBAL (name='%s')\n",name.c_str());
  }

  name = my_name;

  grid = &my_grid;
  dims = GRID_3D_BEDROCK;
  
  PetscInt       M, N, m, n;
  PetscErrorCode ierr;
  const PetscInt *lx, *ly;
  ierr = DAGetOwnershipRanges(my_grid.da2, &ly, &lx, PETSC_NULL); CHKERRQ(ierr);
  ierr = DAGetInfo(my_grid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate3d(my_grid.com, DA_YZPERIODIC, DA_STENCIL_STAR,
		    my_grid.Mbz, N, M,	// P, N, M
		    1, n, m,		// p, n, m
		    1, 1,		// dof, stencil width
                    PETSC_NULL, ly, lx, // lz, ly, lx
		    &da); CHKERRQ(ierr);

  ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);

  localp = false;

  vars[0].init(name, my_grid, GRID_3D_BEDROCK);

  //  ierr = this->set(GSL_NAN); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec3Bedrock::destroy() {
  PetscErrorCode ierr;

  ierr = IceModelVec::destroy(); CHKERRQ(ierr);

  // soundings:
  if (sounding_viewers != NULL) {
    map<string,PetscViewer>::iterator i;
    for (i = (*sounding_viewers).begin(); i != (*sounding_viewers).end(); ++i) {
      if ((*i).second != PETSC_NULL) {
	ierr = PetscViewerDestroy((*i).second); CHKERRQ(ierr);
      }
    }
    delete sounding_viewers;
    sounding_viewers = NULL;
  }

  if (sounding_buffer != PETSC_NULL) {
    ierr = VecDestroy(sounding_buffer); CHKERRQ(ierr);
    sounding_buffer = PETSC_NULL;
  }

  return 0;
}

//! Set values of bedrock scalar quantity at internal levels determined by IceGrid.
/*!
Array \c valsIN must be an allocated array of \c grid->Mbz \c PetscScalar s.
 */
PetscErrorCode  IceModelVec3Bedrock::setInternalColumn(PetscInt i, PetscInt j, PetscScalar *valsIN) {
  
  PetscErrorCode ierr = checkHaveArray();  CHKERRQ(ierr);
  PetscScalar ***arr = (PetscScalar***) array;
  for (PetscInt k = 0; k < grid->Mbz; k++) {
    arr[i][j][k] = valsIN[k];
  }
  return 0;
}


//! Set values of bedrock scalar quantity: set all values in a column to the same value.
PetscErrorCode  IceModelVec3Bedrock::setColumn(PetscInt i, PetscInt j, PetscScalar c) {

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
PetscErrorCode  IceModelVec3Bedrock::getInternalColumn(PetscInt i, PetscInt j, 
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
PetscErrorCode  IceModelVec3Bedrock::setValColumnPL(PetscInt i, PetscInt j,
						    PetscScalar *source) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // check if in ownership ?

  PetscScalar *levels = grid->zblevels;
  PetscScalar *levels_fine = grid->zblevels_fine;
  PetscScalar ***arr = (PetscScalar***) array;

  if (grid->Mbz == 1) {
    arr[i][j][0] = source[0];
    return 0;
  }

  for (PetscInt k=0; k < grid->Mbz; k++) {
    PetscInt m = grid->bed_fine2storage[k];
    
    const PetscScalar increment = (levels[k] - levels_fine[m]) / (levels_fine[m+1] - levels_fine[m]);
    arr[i][j][k] = source[m] +  increment * (source[m+1] - source[m]);
  }

  return 0;
}


PetscErrorCode  IceModelVec3Bedrock::isLegalLevel(PetscScalar z) {
  if (z < -grid->Lbz - 1.0e-6) {
    SETERRQ3(1,
       "level z = %10.8f is below bottom of bedrock at -Lbz = %10.8f; IceModelVec3Bedrock has name='%s'; ENDING!\n",
       z,-grid->Lbz,name.c_str());
  }
  if (z > 0.0 + 1.0e-6) {
    SETERRQ2(2,"level z = %10.8f is above top of bedrock at z=0; IceModelVec3Bedrock has name='%s'; ENDING!\n",
              z,name.c_str());
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
PetscErrorCode  IceModelVec3Bedrock::getValColumnPL(PetscInt i, PetscInt j, PetscScalar *result) {
  
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // check if in ownership ?

  PetscScalar *levels = grid->zblevels;
  PetscScalar *levels_fine = grid->zblevels_fine;
  PetscScalar ***arr = (PetscScalar***) array;

  if (grid->Mbz_fine == 1) {
    result[0] = arr[i][j][0];
    return 0;
  }

  for (PetscInt k = 0; k < grid->Mbz_fine; k++) {
    PetscInt m = grid->bed_storage2fine[k];

    const PetscScalar incr = (levels_fine[k] - levels[m]) / (levels[m+1] - levels[m]);
    const PetscScalar valm = arr[i][j][m];
    result[k] = valm + incr * (arr[i][j][m+1] - valm);
  }

  return 0;
}

//! \brief Return values of bedrock scalar quantity at given levels (m) below
//! the top of the bedrock, using local quadratic interpolation.
/*! Input array \c levelsIN must be an allocated array of \c nlevels scalars
  (\c PetscScalar).

  \c levelsIN must be strictly increasing and in the range
  \f$-\mathtt{grid.Lbz} <= z <= 0\f$.

  Return array \c valsOUT must be an allocated array of \c nlevels scalars (\c
  PetscScalar).

  Upon return, \c valsOUT will be filled with values of scalar quantity at the
  \f$z\f$ values in \c levelsIN.
 */
PetscErrorCode  IceModelVec3Bedrock::getValColumnQUAD(PetscInt i, PetscInt j, PetscScalar *result) {
  
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  // check if in ownership ?

  PetscScalar *levels = grid->zblevels;
  PetscScalar *levels_fine = grid->zblevels_fine;
  PetscScalar ***arr = (PetscScalar***) array;

  if (grid->Mbz_fine == 1) {
    result[0] = arr[i][j][0];
    return 0;
  }

  PetscScalar *values = new PetscScalar[grid->Mbz];

  ierr = PetscMemcpy(values, arr[i][j], grid->Mbz*sizeof(PetscScalar)); CHKERRQ(ierr);

  for (PetscInt k = 0; k < grid->Mbz_fine; k++) {
    PetscInt m = grid->bed_storage2fine[k];

    const PetscScalar z0 = levels[m],
                      f0 = values[m];
    if (m >= grid->Mbz - 2) {
      // just do linear interpolation at top of grid
      const PetscScalar incr = (levels_fine[k] - z0) / (levels[m+1] - z0);
      result[k] = f0 + incr * (values[m+1] - f0);
    } else {
      const PetscScalar dz1 = levels[m+1] - z0,
                        dz2 = levels[m+2] - z0;
      const PetscScalar D1 = (values[m+1] - f0) / dz1,
                        D2 = (values[m+2] - f0) / dz2;
      const PetscScalar c = (D2 - D1) / (dz2 - dz1),
                        b = D1 - c * dz1;
      const PetscScalar s = levels_fine[k] - z0;
      result[k] = f0 + s * (b + c * s);
    }
  }

  delete[] values;

  return 0;
}


//! If the grid is equally spaced in the bedrock then use PL, otherwise use QUAD.
PetscErrorCode  IceModelVec3Bedrock::getValColumn(PetscInt i, PetscInt j,
						  PetscScalar *result) {

  if (grid->bed_vertical_spacing == EQUAL) {
    return getValColumnPL(i, j, result);
  } else {
    return getValColumnQUAD(i, j, result);
  }
}


//! Extends an IceModelVec3 and fills all the new grid points with \c fill_value.
PetscErrorCode IceModelVec3::extend_vertically(int old_Mz, PetscScalar fill_value) {
  PetscErrorCode ierr;

  // Allocate more memory:
  ierr = extend_vertically_private(old_Mz); CHKERRQ(ierr);

  // Fill the new layer:
  PetscScalar ***a;
  ierr = DAVecGetArray(da, v, &a); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      for (PetscInt k = old_Mz; k < grid->Mz; k++)
	a[i][j][k] = fill_value;
    }
  }
  ierr = DAVecRestoreArray(da, v, &a); CHKERRQ(ierr);

  // This communicates the ghosts just to update the new levels. Since this
  // only happens when the grid is extended it should not matter.
  if (localp) {
    ierr = beginGhostComm(); CHKERRQ(ierr);
    ierr = endGhostComm(); CHKERRQ(ierr);
  }

  return 0;
}


//! Extends an IceModelVec3 and fills the new grid points with corresponding \c fill_values values.
PetscErrorCode IceModelVec3::extend_vertically(int old_Mz, IceModelVec2S &fill_values) {
  PetscErrorCode ierr;

  // Allocate more memory:
  ierr = extend_vertically_private(old_Mz); CHKERRQ(ierr);

  // Fill the new layer:
  PetscScalar ***a, **filler;
  ierr = DAVecGetArray(da, v, &a); CHKERRQ(ierr);
  ierr = fill_values.get_array(filler); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      for (PetscInt k = old_Mz; k < grid->Mz; k++)
	a[i][j][k] = filler[i][j];
    }
  }
  ierr = DAVecRestoreArray(da, v, &a); CHKERRQ(ierr);
  ierr = fill_values.end_access(); CHKERRQ(ierr);

  // This communicates the ghosts just to update the new levels. Since this
  // only happens when the grid is extended it should not matter.
  if (localp) {
    ierr = beginGhostComm(); CHKERRQ(ierr);
    ierr = endGhostComm(); CHKERRQ(ierr);
  }

  return 0;
}

//! Handles the memory allocation/deallocation and copying. Does not fill the values of the new layer.
PetscErrorCode IceModelVec3::extend_vertically_private(int old_Mz) {
  PetscErrorCode ierr;
  Vec v_new;
  DA da_new;

  // This code should match what is being done in IceModelVec3::create():

  ierr = create_da(da_new, grid->Mz); CHKERRQ(ierr);
  
  if (localp) {
    ierr = DACreateLocalVector(da_new, &v_new); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da_new, &v_new); CHKERRQ(ierr);
  }

  // Copy all the values from the old Vec to the new one:
  PetscScalar ***a_new;
  PetscScalar ***a_old;
  ierr = DAVecGetArray(da, v, &a_old); CHKERRQ(ierr);
  ierr = DAVecGetArray(da_new, v_new, &a_new); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      for (PetscInt k=0; k < old_Mz; k++)
	a_new[i][j][k] = a_old[i][j][k];
    }
  }
  ierr = DAVecRestoreArray(da, v, &a_old); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da_new, v_new, &a_new); CHKERRQ(ierr);

  // Deallocate old DA and Vec:
  ierr = VecDestroy(v); CHKERRQ(ierr);
  v = v_new;

  ierr = DADestroy(da); CHKERRQ(ierr);
  da = da_new;

  return 0;
}

PetscErrorCode IceModelVec3::view_surface(IceModelVec2S &thickness, PetscInt viewer_size) {
  PetscErrorCode ierr;
  Vec g2;

  ierr = DACreateGlobalVector(grid->da2, &g2); CHKERRQ(ierr);

  if ((*map_viewers)[name] == PETSC_NULL) {
    string title = string_attr("long_name") + " at the ice surface (" +
      string_attr("glaciological_units") + ")";

    ierr = create_viewer(viewer_size, title, (*map_viewers)[name]); CHKERRQ(ierr);
  }

  ierr = getSurfaceValues(g2, thickness); CHKERRQ(ierr);

  ierr = vars[0].to_glaciological_units(g2); CHKERRQ(ierr);

  ierr = VecView(g2, (*map_viewers)[name]); CHKERRQ(ierr);

  ierr = VecDestroy(g2); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec3::view_horizontal_slice(PetscScalar level, PetscInt viewer_size) {
  PetscErrorCode ierr;
  Vec g2;

  ierr = DACreateGlobalVector(grid->da2, &g2); CHKERRQ(ierr);

  if ((*slice_viewers)[name] == PETSC_NULL) {
    ostringstream title;
    title << string_attr("long_name") << " at " << level << " m above the base of ice, (" << 
      string_attr("glaciological_units") << ")";

    ierr = create_viewer(viewer_size, title.str(), (*slice_viewers)[name]); CHKERRQ(ierr);
  }

  ierr = getHorSlice(g2, level); CHKERRQ(ierr);

  ierr = vars[0].to_glaciological_units(g2); CHKERRQ(ierr);

  ierr = VecView(g2, (*slice_viewers)[name]); CHKERRQ(ierr);

  ierr = VecDestroy(g2); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec3::view_sounding(int i, int j, PetscInt viewer_size) {
  PetscErrorCode ierr;
  PetscScalar *ivals;
  PetscInt my_Mz = grid->Mz;

  // memory allocation:
  if (sounding_buffer == PETSC_NULL) {
    ierr = VecCreateMPI(grid->com, PETSC_DECIDE, my_Mz, &sounding_buffer); CHKERRQ(ierr);
  }

  // create the title:
  if ((*sounding_viewers)[name] == PETSC_NULL) {
    string title = string_attr("long_name") + " sounding (" + string_attr("glaciological_units") + ")";

    ierr = create_viewer(viewer_size, title, (*sounding_viewers)[name]); CHKERRQ(ierr);
  }

  // get the sounding:
  if ((i >= grid->xs) && (i < grid->xs + grid->xm) && (j >= grid->ys) && (j < grid->ys + grid->ym)) {
    PetscInt *row = new PetscInt[my_Mz];
    for (PetscInt k = 0; k < my_Mz; k++) row[k] = k;

    ierr = begin_access(); CHKERRQ(ierr);
    ierr = getInternalColumn(i, j, &ivals); CHKERRQ(ierr);
    ierr = VecSetValues(sounding_buffer, my_Mz, row, ivals, INSERT_VALUES); CHKERRQ(ierr);
    ierr = end_access(); CHKERRQ(ierr);

    delete[] row;
  }
  ierr = VecAssemblyBegin(sounding_buffer); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(sounding_buffer); CHKERRQ(ierr);

  // change units:
  ierr = vars[0].to_glaciological_units(sounding_buffer); CHKERRQ(ierr);

  ierr = VecView(sounding_buffer, (*sounding_viewers)[name]); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec3Bedrock::view_sounding(int i, int j, PetscInt viewer_size) {
  PetscErrorCode ierr;
  PetscScalar *ivals = NULL;
  PetscInt my_Mz = grid->Mbz;

  // memory allocation:
  if (sounding_buffer == PETSC_NULL) {
    ierr = VecCreateMPI(grid->com, PETSC_DECIDE, my_Mz, &sounding_buffer); CHKERRQ(ierr);
  }

  // create the title:
  if ((*sounding_viewers)[name] == PETSC_NULL) {
    string title = string_attr("long_name") + " sounding (" + string_attr("glaciological_units") + ")";

    ierr = create_viewer(viewer_size, title, (*sounding_viewers)[name]); CHKERRQ(ierr);
  }

  // get the sounding:
  if ((i >= grid->xs) && (i < grid->xs + grid->xm) && (j >= grid->ys) && (j < grid->ys + grid->ym)) {
    PetscInt *row = new PetscInt[my_Mz];
    for (PetscInt k = 0; k < my_Mz; k++) row[k] = k;

    ierr = begin_access(); CHKERRQ(ierr);
    ierr = getInternalColumn(i, j, &ivals); CHKERRQ(ierr);
    ierr = VecSetValues(sounding_buffer, my_Mz, row, ivals, INSERT_VALUES); CHKERRQ(ierr);
    ierr = end_access(); CHKERRQ(ierr);

    delete[] row;
  }
  ierr = VecAssemblyBegin(sounding_buffer); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(sounding_buffer); CHKERRQ(ierr);

  // change units:
  ierr = vars[0].to_glaciological_units(sounding_buffer); CHKERRQ(ierr);

  ierr = VecView(sounding_buffer, (*sounding_viewers)[name]); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode  IceModelVec3::has_nan() {
  PetscErrorCode ierr;
  vector<PetscReal> V;
  V.resize(grid->Mz);
  PetscReal *tmp = &V[0];

  ierr = begin_access(); CHKERRQ(ierr);
  PetscInt i, j, k;
  for (i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (j=grid->ys; j<grid->ys+grid->ym; ++j) {
      ierr = getInternalColumn(i, j, &tmp); CHKERRQ(ierr);
      for (k = 0; k < grid->Mz; k++) {
	if (gsl_isnan(tmp[k])) {
	  ierr = PetscSynchronizedPrintf(grid->com, "IceModelVec3 %s: NAN (or uninitialized) at i = %d, j = %d\n",
					 name.c_str(), i, j); CHKERRQ(ierr);
	  break;
	}
      }
    }
  }
  ierr = end_access(); CHKERRQ(ierr);

  ierr = PetscSynchronizedFlush(grid->com); CHKERRQ(ierr);

  return 0;
}

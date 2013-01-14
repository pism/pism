// Copyright (C) 2008--2013 Ed Bueler and Constantine Khroulev
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
#include <petscdmda.h>

#include "PIO.hh"
#include "iceModelVec.hh"
#include "IceGrid.hh"

// this file contains method for derived class IceModelVec3

// methods for base class IceModelVec and derived class IceModelVec2S
// are in "iceModelVec.cc"

IceModelVec3D::IceModelVec3D() : IceModelVec() {
  sounding_buffer = PETSC_NULL;
  sounding_viewers = new map<string, PetscViewer>;
}

IceModelVec3D::~IceModelVec3D() {
  if (!shallow_copy) {
    destroy();
  }
}

IceModelVec3D::IceModelVec3D(const IceModelVec3D &other)
  : IceModelVec(other) {
  sounding_buffer = other.sounding_buffer;
  sounding_viewers = other.sounding_viewers;
  shallow_copy = true;
}

//! Allocate a DA and a Vec from information in IceGrid.
PetscErrorCode  IceModelVec3D::allocate(IceGrid &my_grid, string my_name,
                                        bool local, vector<double> levels, int stencil_width) {
  PetscErrorCode ierr;
  if (!utIsInit()) {
    SETERRQ(grid->com, 1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(grid->com, 1,"IceModelVec3 with name='%s' already allocated\n",name.c_str());
  }
  
  grid = &my_grid;

  zlevels = levels;
  n_levels = (int)zlevels.size();
  da_stencil_width = stencil_width;

  ierr = grid->get_dm(this->n_levels, this->da_stencil_width, da); CHKERRQ(ierr);

  if (local) {
    ierr = DMCreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DMCreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  name = my_name;

  vars[0].init_3d(my_name, my_grid, zlevels);

  //  ierr = this->set(GSL_NAN); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec3D::destroy() {
  PetscErrorCode ierr;
  map<string,PetscViewer>::iterator i;

  // soundings:
  if (sounding_viewers != NULL) {
    for (i = (*sounding_viewers).begin(); i != (*sounding_viewers).end(); ++i) {
      if ((*i).second != PETSC_NULL) {
	ierr = PetscViewerDestroy(&(*i).second); CHKERRQ(ierr);
      }
    }
    delete sounding_viewers;
    sounding_viewers = NULL;
  }

  if (sounding_buffer != PETSC_NULL) {
    ierr = VecDestroy(&sounding_buffer); CHKERRQ(ierr);
    sounding_buffer = PETSC_NULL;
  }

  return 0;
}

PetscErrorCode  IceModelVec3D::isLegalLevel(PetscScalar z) {
  double z_min = zlevels.front(),
    z_max = zlevels.back();
  if (z < z_min - 1.0e-6 || z > z_max + 1.0e-6) {
    SETERRQ2(grid->com, 1,"level z = %5.4f is outside the valid range;\n"
               "  IceModelVec3 has name='%s'; ENDING!\n",
              z,name.c_str());
  }
  return 0;
}


//! Set values of an ice scalar quantity in a column by linear <i>interpolation</i>.
/*!
  Input array \c source and \c must contain \c grid.Mz_fine scalars
  (\c PetscScalar).  Upon completion, internal storage will hold values derived from 
  linearly interpolating the input values.
 */
PetscErrorCode  IceModelVec3::setValColumnPL(PetscInt i, PetscInt j, PetscScalar *source) {
#if (PISM_DEBUG==1)
  PetscErrorCode ierr = checkAllocated(); CHKERRQ(ierr);
  check_array_indices(i, j);
#endif

  vector<double> &zlevels_fine = grid->zlevels_fine;

  PetscScalar ***arr = (PetscScalar***) array;
  
  for (PetscInt k=0; k < n_levels; ++k) {
    PetscInt m = grid->ice_fine2storage[k];

    const PetscScalar increment = (zlevels[k] - zlevels_fine[m])
                                  / (zlevels_fine[m+1] - zlevels_fine[m]);
    arr[i][j][k] = source[m] +  increment * (source[m+1] - source[m]);
  }

  return 0;
}


//! Set all values of scalar quantity to given a single value in a particular column.
PetscErrorCode IceModelVec3D::setColumn(PetscInt i, PetscInt j, PetscScalar c) {
  PetscErrorCode ierr;
#if (PISM_DEBUG==1)
  ierr = checkHaveArray();  CHKERRQ(ierr);
  check_array_indices(i, j);
#endif

  PetscScalar ***arr = (PetscScalar***) array;

  if (c == 0.0) {
    ierr = PetscMemzero(arr[i][j], n_levels * sizeof(PetscScalar)); CHKERRQ(ierr);
  } else {
    for (PetscInt k=0; k < n_levels; k++) {
      arr[i][j][k] = c;
    }
  }
  return 0;
}


//! Return value of scalar quantity at level z (m) above base of ice (by linear interpolation).
PetscScalar IceModelVec3D::getValZ(PetscInt i, PetscInt j, PetscScalar z) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j);

  if (checkHaveArray() != 0) {
    PetscPrintf(PETSC_COMM_SELF, 
       "IceModelVec3 getValZ(): array was not allocated (so says\n"
       "  IceModelVec::checkHaveArray()); name = %s\n", name.c_str());
    PISMEnd();
  }
  if (isLegalLevel(z) != 0) {
    PetscPrintf(PETSC_COMM_SELF, 
       "IceModelVec3 getValZ(): isLegalLevel() says level %f was\n"
       "  not legal; name = %s\n", z, name.c_str());
    PISMEnd();
  }
#endif

  PetscScalar ***arr = (PetscScalar***) array;
  if (z >= zlevels.back())
    return arr[i][j][n_levels - 1];
  else if (z <= zlevels.front())
    return arr[i][j][0];

  PetscInt mcurr = 0;
  while (zlevels[mcurr+1] < z) mcurr++;

  const PetscScalar incr = (z - zlevels[mcurr]) / (zlevels[mcurr+1] - zlevels[mcurr]);
  const PetscScalar valm = arr[i][j][mcurr];
  return valm + incr * (arr[i][j][mcurr+1] - valm);
}


//! Return values on planar star stencil of scalar quantity at level z (by linear interpolation).
PetscErrorCode   IceModelVec3::getPlaneStarZ(PetscInt i, PetscInt j, PetscScalar z,
					     planeStar<PetscScalar> *star) {
#if (PISM_DEBUG==1)
  PetscErrorCode ierr;
  ierr = checkHaveArray();  CHKERRQ(ierr);
  ierr = isLegalLevel(z);  CHKERRQ(ierr);
  // check ownership here?
  if (!localp) {
    SETERRQ1(grid->com, 1,"IceModelVec3 ERROR: IceModelVec3 with name='%s' is GLOBAL\n"
               "  and cannot do getPlaneStarZ()\n", name.c_str());
  }
  check_array_indices(i, j);
#endif

  PetscInt     kbz;
  PetscScalar  incr;
  if (z >= zlevels.back()) {
    kbz = n_levels - 1;
    incr = 0.0;
  } else if (z <= zlevels.front()) {
    kbz = 0;
    incr = 0.0;
  } else {
    kbz = 0;
    while (zlevels[kbz+1] < z)
      kbz++;

    incr = (z - zlevels[kbz]) / (zlevels[kbz+1] - zlevels[kbz]);
  }

  PetscScalar ***arr = (PetscScalar***) array;

  if (kbz < n_levels - 1) {
    star->ij  = arr[i][j][kbz]   + incr * (arr[i][j][kbz + 1]   - arr[i][j][kbz]);
    star->e = arr[i+1][j][kbz] + incr * (arr[i+1][j][kbz + 1] - arr[i+1][j][kbz]);
    star->w = arr[i-1][j][kbz] + incr * (arr[i-1][j][kbz + 1] - arr[i-1][j][kbz]);
    star->n = arr[i][j+1][kbz] + incr * (arr[i][j+1][kbz + 1] - arr[i][j+1][kbz]);
    star->s = arr[i][j-1][kbz] + incr * (arr[i][j-1][kbz + 1] - arr[i][j-1][kbz]);
  } else {
    star->ij  = arr[i][j][kbz];
    star->e = arr[i+1][j][kbz];
    star->w = arr[i-1][j][kbz];
    star->n = arr[i][j+1][kbz];
    star->s = arr[i][j-1][kbz];
  }

  return 0;
}

//! Gets a map-plane star stencil directly from the storage grid.
PetscErrorCode IceModelVec3::getPlaneStar(PetscInt i, PetscInt j, PetscInt k,
						  planeStar<PetscScalar> *star) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j);
#endif

  PetscScalar ***arr = (PetscScalar***) array;

  star->ij  = arr[i][j][k];
  star->e = arr[i+1][j][k];
  star->w = arr[i-1][j][k];
  star->n = arr[i][j+1][k];
  star->s = arr[i][j-1][k];

  return 0;
}

//! Gets a map-plane star stencil on the fine vertical grid.
PetscErrorCode IceModelVec3::getPlaneStar_fine(PetscInt i, PetscInt j, PetscInt k,
					       planeStar<PetscScalar> *star) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j);
#endif

  PetscInt kbz = grid->ice_storage2fine[k];

  if (kbz < n_levels - 1) {
    PetscScalar z = grid->zlevels_fine[k],
      incr = (z - zlevels[kbz]) / (zlevels[kbz+1] - zlevels[kbz]);
    PetscScalar ***arr = (PetscScalar***) array;

    star->ij  = arr[i][j][kbz]   + incr * (arr[i][j][kbz + 1]   - arr[i][j][kbz]);
    star->e = arr[i+1][j][kbz] + incr * (arr[i+1][j][kbz + 1] - arr[i+1][j][kbz]);
    star->w = arr[i-1][j][kbz] + incr * (arr[i-1][j][kbz + 1] - arr[i-1][j][kbz]);
    star->n = arr[i][j+1][kbz] + incr * (arr[i][j+1][kbz + 1] - arr[i][j+1][kbz]);
    star->s = arr[i][j-1][kbz] + incr * (arr[i][j-1][kbz + 1] - arr[i][j-1][kbz]);
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
#if (PISM_DEBUG==1)
  PetscErrorCode ierr = checkAllocated(); CHKERRQ(ierr);
  check_array_indices(i, j);
#endif

  vector<double> &zlevels_fine = grid->zlevels_fine;
  PetscScalar ***arr = (PetscScalar***) array;

  for (PetscInt k = 0; k < grid->Mz_fine; k++) {
    if (k > ks) {
      result[k] = arr[i][j][grid->ice_storage2fine[k]];
      continue;
    }

    PetscInt m = grid->ice_storage2fine[k];

    // extrapolate (if necessary):
    if (m == n_levels - 1) {
      result[k] = arr[i][j][n_levels-1];
      continue;
    }

    const PetscScalar incr = (zlevels_fine[k] - zlevels[m]) / (zlevels[m+1] - zlevels[m]);
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
#if (PISM_DEBUG==1)
  check_array_indices(i, j);
#endif

  vector<double> &zlevels_fine = grid->zlevels_fine;
  const PetscScalar ***arr = (const PetscScalar***) array;
  const PetscScalar *column = arr[i][j];

  for (PetscInt k = 0; k < grid->Mz_fine; k++) {
    if (k > ks) {
      result[k] = arr[i][j][grid->ice_storage2fine[k]];
      continue;
    }

    const PetscInt m = grid->ice_storage2fine[k];

    // extrapolate (if necessary):
    if (m == n_levels - 1) {
      result[k] = column[n_levels-1];
      continue;
    }

    const PetscScalar z0 = zlevels[m],
                      f0 = column[m];
    if (m >= n_levels - 2) {
				// top of the grid: just do linear interpolation
      const PetscScalar incr = (zlevels_fine[k] - z0) / (zlevels[m+1] - z0);
      result[k] = f0 + incr * (column[m+1] - f0);
    } else {			// the rest: one-sided quadratic interpolation
      const PetscScalar dz1 = zlevels[m+1] - z0,
                        dz2 = zlevels[m+2] - z0;
      const PetscScalar D1 = (column[m+1] - f0) / dz1,
                        D2 = (column[m+2] - f0) / dz2;
      const PetscScalar c = (D2 - D1) / (dz2 - dz1),
                        b = D1 - c * dz1;
      const PetscScalar s = zlevels_fine[k] - z0;
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
/*!
 * FIXME: this method is misnamed: the slice is horizontal in the PISM
 * coordinate system, not in reality.
 */
PetscErrorCode  IceModelVec3::getHorSlice(Vec &gslice, PetscScalar z) {
  PetscErrorCode ierr;
  PetscScalar    **slice_val;

  DM da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da2, gslice, &slice_val); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      slice_val[i][j] = getValZ(i,j,z);
    }
  }
  ierr = DMDAVecRestoreArray(da2, gslice, &slice_val); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! Copies a horizontal slice at level z of an IceModelVec3 into an IceModelVec2S gslice.
/*!
 * FIXME: this method is misnamed: the slice is horizontal in the PISM
 * coordinate system, not in reality.
 */
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

  DM da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da2, gsurf, &surf_val); CHKERRQ(ierr);
  ierr = myH.get_array(H); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      surf_val[i][j] = getValZ(i,j,H[i][j]);
    }
  }
  ierr = myH.end_access(); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da2, gsurf, &surf_val); CHKERRQ(ierr);
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


PetscErrorCode  IceModelVec3D::getInternalColumn(PetscInt i, PetscInt j, PetscScalar **valsPTR) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j);
#endif
  PetscScalar ***arr = (PetscScalar***) array;
  *valsPTR = arr[i][j];
  return 0;
}


PetscErrorCode  IceModelVec3D::setInternalColumn(PetscInt i, PetscInt j, PetscScalar *valsIN) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j);
#endif
  PetscScalar ***arr = (PetscScalar***) array;
  PetscErrorCode ierr = PetscMemcpy(arr[i][j], valsIN, n_levels*sizeof(PetscScalar));
  CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::create(IceGrid &my_grid, string my_name, bool local,
                                     int stencil_width) {

  PetscErrorCode ierr = IceModelVec3D::allocate(my_grid, my_name, local,
                                                my_grid.zlevels, stencil_width); CHKERRQ(ierr);

  return 0;
}

//! Extends an IceModelVec3 and fills all the new grid points with \c fill_value.
PetscErrorCode IceModelVec3::extend_vertically(int old_Mz, PetscScalar fill_value) {
  PetscErrorCode ierr;

  // Allocate more memory:
  ierr = extend_vertically_private(old_Mz); CHKERRQ(ierr);

  // Fill the new layer:
  PetscScalar ***a;
  ierr = DMDAVecGetArrayDOF(da, v, &a); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      for (PetscInt k = old_Mz; k < n_levels; k++)
	a[i][j][k] = fill_value;
    }
  }
  ierr = DMDAVecRestoreArrayDOF(da, v, &a); CHKERRQ(ierr);

  // This communicates the ghosts just to update the new levels. Since this
  // only happens when the grid is extended it should not matter.
  if (localp) {
    ierr = update_ghosts(); CHKERRQ(ierr);
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
  ierr = DMDAVecGetArrayDOF(da, v, &a); CHKERRQ(ierr);
  ierr = fill_values.get_array(filler); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      for (PetscInt k = old_Mz; k < n_levels; k++)
	a[i][j][k] = filler[i][j];
    }
  }
  ierr = DMDAVecRestoreArrayDOF(da, v, &a); CHKERRQ(ierr);
  ierr = fill_values.end_access(); CHKERRQ(ierr);

  // This communicates the ghosts just to update the new levels. Since this
  // only happens when the grid is extended it should not matter.
  if (localp) {
    ierr = update_ghosts(); CHKERRQ(ierr);
  }

  return 0;
}

//! Handles the memory allocation/deallocation and copying. Does not fill the values of the new layer.
PetscErrorCode IceModelVec3::extend_vertically_private(int old_Mz) {
  PetscErrorCode ierr;
  Vec v_new;
  DM da_new;

  // This code should match what is being done in IceModelVec3D::allocate():

  zlevels = grid->zlevels;
  n_levels = (int)zlevels.size();
  for (int i = 0; i < dof; ++i)
    vars[0].set_levels(zlevels);

  ierr = grid->get_dm(this->n_levels, this->da_stencil_width, da_new); CHKERRQ(ierr);

  if (localp) {
    ierr = DMCreateLocalVector(da_new, &v_new); CHKERRQ(ierr);
  } else {
    ierr = DMCreateGlobalVector(da_new, &v_new); CHKERRQ(ierr);
  }

  // Copy all the values from the old Vec to the new one:
  PetscScalar ***a_new;
  PetscScalar ***a_old;
  ierr = DMDAVecGetArrayDOF(da, v, &a_old); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_new, v_new, &a_new); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      for (PetscInt k=0; k < old_Mz; k++)
	a_new[i][j][k] = a_old[i][j][k];
    }
  }
  ierr = DMDAVecRestoreArrayDOF(da, v, &a_old); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da_new, v_new, &a_new); CHKERRQ(ierr);

  // Deallocate old Vec:
  ierr = VecDestroy(&v); CHKERRQ(ierr);
  v = v_new;
  da = da_new;

  // IceGrid will dispose of the old DA

  // de-allocate the sounding buffer because we'll need a bigger one
  if (sounding_buffer != NULL) {
    ierr = VecDestroy(&sounding_buffer); CHKERRQ(ierr);
    sounding_buffer = PETSC_NULL;
  }

  return 0;
}

PetscErrorCode IceModelVec3D::view_sounding(int i, int j, PetscInt viewer_size) {
  PetscErrorCode ierr;

  // create the title:
  if ((*sounding_viewers)[name] == PETSC_NULL) {
    string title = string_attr("long_name") + " sounding (" + string_attr("glaciological_units") + ")";

    ierr = grid->create_viewer(viewer_size, title, (*sounding_viewers)[name]); CHKERRQ(ierr);
  }

  ierr = view_sounding(i, j, (*sounding_viewers)[name]); CHKERRQ(ierr);

  return 0;
}

//! \brief View a sounding using an existing PETSc viewer.
PetscErrorCode IceModelVec3D::view_sounding(int i, int j, PetscViewer my_viewer) {
  PetscErrorCode ierr;
  PetscScalar *ivals;

#if (PISM_DEBUG==1)
    check_array_indices(i, j);
#endif

  const string tname = string_attr("long_name"),
    tunits = " (" + string_attr("glaciological_units") + ")",
    title = tname + tunits;

  PetscDraw draw;
  ierr = PetscViewerDrawGetDraw(my_viewer, 0, &draw); CHKERRQ(ierr);
  ierr = PetscDrawSetTitle(draw, title.c_str()); CHKERRQ(ierr);

  // memory allocation:
  if (sounding_buffer == PETSC_NULL) {
    ierr = VecCreateMPI(grid->com, PETSC_DECIDE, n_levels, &sounding_buffer); CHKERRQ(ierr);
  }

  // get the sounding:
  if ((i >= grid->xs) && (i < grid->xs + grid->xm) && (j >= grid->ys) && (j < grid->ys + grid->ym)) {
    PetscInt *row = new PetscInt[n_levels];
    for (PetscInt k = 0; k < n_levels; k++) row[k] = k;

    ierr = begin_access(); CHKERRQ(ierr);
    ierr = getInternalColumn(i, j, &ivals); CHKERRQ(ierr);
    ierr = VecSetValues(sounding_buffer, n_levels, row, ivals, INSERT_VALUES); CHKERRQ(ierr);
    ierr = end_access(); CHKERRQ(ierr);

    delete[] row;
  }
  ierr = VecAssemblyBegin(sounding_buffer); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(sounding_buffer); CHKERRQ(ierr);

  // change units:
  ierr = vars[0].to_glaciological_units(sounding_buffer); CHKERRQ(ierr);

  ierr = VecView(sounding_buffer, my_viewer); CHKERRQ(ierr);
  
  return 0;
}

//! Checks if the current IceModelVec3 has NANs and reports if it does.
/*! Up to a fixed number of messages are printed at stdout.  Returns the full
 count of NANs (which is a nonzero) on this rank. */
PetscErrorCode  IceModelVec3D::has_nan() {
  PetscErrorCode ierr;
  PetscScalar *tmp;
  PetscInt retval=0;
  const PetscInt max_print_this_rank=10;

  ierr = begin_access(); CHKERRQ(ierr);
  PetscInt i, j, k;
  for (i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (j=grid->ys; j<grid->ys+grid->ym; ++j) {
      ierr = getInternalColumn(i, j, &tmp); CHKERRQ(ierr);
      for (k = 0; k < n_levels; k++) {
	if (gsl_isnan(tmp[k])) {
	  retval++;
          if (retval <= max_print_this_rank) {
            ierr = PetscSynchronizedPrintf(grid->com, 
               "IceModelVec3 %s: NAN (or uninitialized) at i = %d, j = %d, k = %d on rank = %d\n",
               name.c_str(), i, j, k, grid->rank); CHKERRQ(ierr);
          }
	  break;
	}
      }
    }
  }
  ierr = end_access(); CHKERRQ(ierr);

  if (retval > 0) {
    ierr = PetscSynchronizedPrintf(grid->com, 
       "IceModelVec3 %s: detected %d NANs (or uninitialized) on rank = %d\n",
             name.c_str(), retval, grid->rank); CHKERRQ(ierr);
  }

  ierr = PetscSynchronizedFlush(grid->com); CHKERRQ(ierr);

  return retval;
}

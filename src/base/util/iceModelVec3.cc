// Copyright (C) 2008--2014 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include <gsl/gsl_math.h>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <petscdmda.h>

#include "PIO.hh"
#include "iceModelVec.hh"
#include "IceGrid.hh"

#include <assert.h>

// this file contains method for derived class IceModelVec3

// methods for base class IceModelVec and derived class IceModelVec2S
// are in "iceModelVec.cc"

IceModelVec3D::IceModelVec3D() : IceModelVec() {
}

IceModelVec3D::~IceModelVec3D() {
  destroy();
}

//! Allocate a DA and a Vec from information in IceGrid.
PetscErrorCode  IceModelVec3D::allocate(IceGrid &my_grid, std::string my_name,
                                        IceModelVecKind ghostedp, std::vector<double> levels,
                                        unsigned int stencil_width) {
  PetscErrorCode ierr;

  assert(v == NULL);
  
  grid = &my_grid;

  zlevels = levels;
  m_n_levels = (unsigned int)zlevels.size();
  m_da_stencil_width = stencil_width;

  ierr = grid->get_dm(this->m_n_levels, this->m_da_stencil_width, m_da); CHKERRQ(ierr);

  m_has_ghosts = (ghostedp == WITH_GHOSTS);

  if (m_has_ghosts == true) {
    ierr = DMCreateLocalVector(m_da, &v); CHKERRQ(ierr);
  } else {
    ierr = DMCreateGlobalVector(m_da, &v); CHKERRQ(ierr);
  }

  m_name = my_name;

  m_metadata.resize(m_dof, NCSpatialVariable(grid->get_unit_system()));
  m_metadata[0].init_3d(my_name, my_grid, zlevels);

  return 0;
}

PetscErrorCode  IceModelVec3D::isLegalLevel(double z) {
  double z_min = zlevels.front(),
    z_max = zlevels.back();
  if (z < z_min - 1.0e-6 || z > z_max + 1.0e-6) {
    SETERRQ2(grid->com, 1,"level z = %5.4f is outside the valid range;\n"
               "  IceModelVec3 has name='%s'; ENDING!\n",
              z,m_name.c_str());
  }
  return 0;
}


//! Set values of an ice scalar quantity in a column by linear *interpolation*.
/*!
  Input array `source` and `must` contain `grid.Mz_fine` scalars
  (`double`).  Upon completion, internal storage will hold values derived from 
  linearly interpolating the input values.
 */
PetscErrorCode  IceModelVec3::setValColumnPL(int i, int j, double *source) {
#if (PISM_DEBUG==1)
  assert(v != NULL);
  check_array_indices(i, j, 0);
#endif

  std::vector<double> &zlevels_fine = grid->zlevels_fine;

  double ***arr = (double***) array;
  
  for (unsigned int k=0; k < m_n_levels-1; ++k) {
    int m = grid->ice_fine2storage[k];

    const double increment = (zlevels[k] - zlevels_fine[m])
                                  / (zlevels_fine[m+1] - zlevels_fine[m]);
    arr[i][j][k] = source[m] +  increment * (source[m+1] - source[m]);
  }
  
  arr[i][j][m_n_levels-1] = source[grid->ice_fine2storage[m_n_levels-1]];

  return 0;
}


//! Set all values of scalar quantity to given a single value in a particular column.
PetscErrorCode IceModelVec3D::setColumn(int i, int j, double c) {
  PetscErrorCode ierr;
#if (PISM_DEBUG==1)
  assert(array != NULL);
  check_array_indices(i, j, 0);
#endif

  double ***arr = (double***) array;

  if (c == 0.0) {
    ierr = PetscMemzero(arr[i][j], m_n_levels * sizeof(double)); CHKERRQ(ierr);
  } else {
    for (unsigned int k=0; k < m_n_levels; k++) {
      arr[i][j][k] = c;
    }
  }
  return 0;
}


//! Return value of scalar quantity at level z (m) above base of ice (by linear interpolation).
double IceModelVec3D::getValZ(int i, int j, double z) {
#if (PISM_DEBUG==1)
  assert(array != NULL);
  check_array_indices(i, j, 0);

  if (isLegalLevel(z) != 0) {
    PetscPrintf(PETSC_COMM_SELF, 
       "IceModelVec3 getValZ(): isLegalLevel() says level %f was\n"
       "  not legal; name = %s\n", z, m_name.c_str());
    PISMEnd();
  }
#endif

  double ***arr = (double***) array;
  if (z >= zlevels.back())
    return arr[i][j][m_n_levels - 1];
  else if (z <= zlevels.front())
    return arr[i][j][0];

  int mcurr = 0;
  while (zlevels[mcurr+1] < z) mcurr++;

  const double incr = (z - zlevels[mcurr]) / (zlevels[mcurr+1] - zlevels[mcurr]);
  const double valm = arr[i][j][mcurr];
  return valm + incr * (arr[i][j][mcurr+1] - valm);
}


//! Return values on planar star stencil of scalar quantity at level z (by linear interpolation).
PetscErrorCode   IceModelVec3::getPlaneStarZ(int i, int j, double z,
                                             planeStar<double> *star) {
#if (PISM_DEBUG==1)
  assert(array != NULL);
  assert(m_has_ghosts == true);
  PetscErrorCode ierr = isLegalLevel(z);  CHKERRQ(ierr);
  check_array_indices(i, j, 0);
#endif

  unsigned int kbz = 0;
  double incr = 0.0;
  if (z >= zlevels.back()) {
    kbz = m_n_levels - 1;
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

  double ***arr = (double***) array;

  if (kbz < m_n_levels - 1) {
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
PetscErrorCode IceModelVec3::getPlaneStar(int i, int j, unsigned int k,
                                          planeStar<double> *star) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif

  double ***arr = (double***) array;

  star->ij  = arr[i][j][k];
  star->e = arr[i+1][j][k];
  star->w = arr[i-1][j][k];
  star->n = arr[i][j+1][k];
  star->s = arr[i][j-1][k];

  return 0;
}

//! Gets a map-plane star stencil on the fine vertical grid.
PetscErrorCode IceModelVec3::getPlaneStar_fine(int i, int j, unsigned int k,
                                               planeStar<double> *star) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif

  unsigned int kbz = grid->ice_storage2fine[k];

  if (kbz < m_n_levels - 1) {
    double z = grid->zlevels_fine[k],
      incr = (z - zlevels[kbz]) / (zlevels[kbz+1] - zlevels[kbz]);
    double ***arr = (double***) array;

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

//! \brief Return values of ice scalar quantity at given levels (m)
//! above base of ice, using piecewise linear interpolation.
/*!
 * ks is the top-most fine vertical grid level within the ice
 */
PetscErrorCode IceModelVec3::getValColumnPL(int i, int j, unsigned int ks,
                                            double *result) {
#if (PISM_DEBUG==1)
  assert(v != NULL);
  check_array_indices(i, j, 0);
#endif

  std::vector<double> &zlevels_fine = grid->zlevels_fine;
  double ***arr = (double***) array;

  for (unsigned int k = 0; k < grid->Mz_fine; k++) {
    if (k > ks) {
      result[k] = arr[i][j][grid->ice_storage2fine[k]];
      continue;
    }

    unsigned int m = grid->ice_storage2fine[k];

    // extrapolate (if necessary):
    if (m == m_n_levels - 1) {
      result[k] = arr[i][j][m_n_levels-1];
      continue;
    }

    const double incr = (zlevels_fine[k] - zlevels[m]) / (zlevels[m+1] - zlevels[m]);
    const double valm = arr[i][j][m];
    result[k] = valm + incr * (arr[i][j][m+1] - valm);
  }

  return 0;
}

//! \brief Return values of ice scalar quantity on the fine
//! computational grid, using local quadratic interpolation.
PetscErrorCode  IceModelVec3::getValColumnQUAD(int i, int j, unsigned int ks,
                                               double *result) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif

  // Assume that the fine grid is equally-spaced:
  const double dz_fine = grid->zlevels_fine[1] - grid->zlevels_fine[0];
  const double *column = static_cast<const double***>(array)[i][j];

  unsigned int k = 0, m = 0;
  for (m = 0; m < m_n_levels - 2; m++) {
    if (k > ks)
      break;

    const double
      z0 = zlevels[m],
      z1 = zlevels[m+1],
      z2 = zlevels[m+2],
      f0 = column[m],
      f1 = column[m+1],
      f2 = column[m+2];

    const double
      d1 = (f1 - f0) / (z1 - z0),
      d2 = (f2 - f0) / (z2 - z0),
      b  = (d2 - d1) / (z2 - z1),
      a  = d1 - b * (z1 - z0),
      c  = f0;

    double z_fine = k * dz_fine;
    while (z_fine < z1) {
      if (k > ks)
        break;

      const double s = z_fine - z0;

      result[k] = s * (a + b * s) + c;

      k++;
      z_fine = k * dz_fine;
    }
  } // m-loop

  // check if we got to the end of the m-loop and use linear
  // interpolation between the remaining 2 coarse levels
  if (m == m_n_levels - 2) {
    const double
      z0 = zlevels[m],
      z1 = zlevels[m+1],
      f0 = column[m],
      f1 = column[m+1],
      lambda = (f1 - f0) / (z1 - z0);

    double z_fine = k * dz_fine;
    while (z_fine < z1) {
      result[k] = f0 + lambda * (z_fine - z0);

      k++;
      z_fine = k * dz_fine;
    }
  }

  // fill the rest using constant extrapolation
  const double f0 = column[m_n_levels - 1];
  while (k <= ks) {
    result[k] = f0;
    k++;
  }

  return 0;
}


//! If the grid is equally spaced in the ice then use PL, otherwise use QUAD.
PetscErrorCode  IceModelVec3::getValColumn(int i, int j, unsigned int ks,
                                           double *result) {
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
PetscErrorCode  IceModelVec3::getHorSlice(Vec &gslice, double z) {
  PetscErrorCode ierr;
  double    **slice_val;

  DM da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da2, gslice, &slice_val); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (int j=grid->ys; j<grid->ys+grid->ym; j++) {
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
PetscErrorCode  IceModelVec3::getHorSlice(IceModelVec2S &gslice, double z) {
  PetscErrorCode ierr;
  double    **slice_val;

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = gslice.get_array(slice_val); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (int j=grid->ys; j<grid->ys+grid->ym; j++) {
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
  double    **H;
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
  double    **H, **surf_val;

  DM da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da2, gsurf, &surf_val); CHKERRQ(ierr);
  ierr = myH.get_array(H); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (int j=grid->ys; j<grid->ys+grid->ym; j++) {
      surf_val[i][j] = getValZ(i,j,H[i][j]);
    }
  }
  ierr = myH.end_access(); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da2, gsurf, &surf_val); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::getSurfaceValues(IceModelVec2S &gsurf, double **H) {
  PetscErrorCode ierr;
  double    **surf_val;

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = gsurf.get_array(surf_val); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (int j=grid->ys; j<grid->ys+grid->ym; j++) {
      surf_val[i][j] = getValZ(i,j,H[i][j]);
    }
  }
  ierr = gsurf.end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode  IceModelVec3D::getInternalColumn(int i, int j, double **valsPTR) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  double ***arr = (double***) array;
  *valsPTR = arr[i][j];
  return 0;
}


PetscErrorCode  IceModelVec3D::setInternalColumn(int i, int j, double *valsIN) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  double ***arr = (double***) array;
  PetscErrorCode ierr = PetscMemcpy(arr[i][j], valsIN, m_n_levels*sizeof(double));
  CHKERRQ(ierr);
  return 0;
}


PetscErrorCode  IceModelVec3::create(IceGrid &my_grid, std::string my_name, IceModelVecKind ghostedp,
                                     unsigned int stencil_width) {

  PetscErrorCode ierr = IceModelVec3D::allocate(my_grid, my_name, ghostedp,
                                                my_grid.zlevels, stencil_width); CHKERRQ(ierr);

  return 0;
}

//! Extends an IceModelVec3 and fills all the new grid points with `fill_value`.
PetscErrorCode IceModelVec3::extend_vertically(int old_Mz, double fill_value) {
  PetscErrorCode ierr;

  // Allocate more memory:
  ierr = extend_vertically_private(old_Mz); CHKERRQ(ierr);

  // Fill the new layer:
  double ***a;
  ierr = DMDAVecGetArrayDOF(m_da, v, &a); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (int j=grid->ys; j<grid->ys+grid->ym; j++) {
      for (unsigned int k = old_Mz; k < m_n_levels; k++)
        a[i][j][k] = fill_value;
    }
  }
  ierr = DMDAVecRestoreArrayDOF(m_da, v, &a); CHKERRQ(ierr);

  // This communicates the ghosts just to update the new levels. Since this
  // only happens when the grid is extended it should not matter.
  if (m_has_ghosts) {
    ierr = update_ghosts(); CHKERRQ(ierr);
  }

  return 0;
}


//! Extends an IceModelVec3 and fills the new grid points with corresponding `fill_values` values.
PetscErrorCode IceModelVec3::extend_vertically(int old_Mz, IceModelVec2S &fill_values) {
  PetscErrorCode ierr;

  // Allocate more memory:
  ierr = extend_vertically_private(old_Mz); CHKERRQ(ierr);

  // Fill the new layer:
  double ***a, **filler;
  ierr = DMDAVecGetArrayDOF(m_da, v, &a); CHKERRQ(ierr);
  ierr = fill_values.get_array(filler); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (int j=grid->ys; j<grid->ys+grid->ym; j++) {
      for (unsigned int k = old_Mz; k < m_n_levels; k++)
        a[i][j][k] = filler[i][j];
    }
  }
  ierr = DMDAVecRestoreArrayDOF(m_da, v, &a); CHKERRQ(ierr);
  ierr = fill_values.end_access(); CHKERRQ(ierr);

  // This communicates the ghosts just to update the new levels. Since this
  // only happens when the grid is extended it should not matter.
  if (m_has_ghosts) {
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
  m_n_levels = (unsigned int)zlevels.size();
  for (unsigned int i = 0; i < m_dof; ++i)
    m_metadata[0].set_levels(zlevels);

  ierr = grid->get_dm(this->m_n_levels, this->m_da_stencil_width, da_new); CHKERRQ(ierr);

  if (m_has_ghosts) {
    ierr = DMCreateLocalVector(da_new, &v_new); CHKERRQ(ierr);
  } else {
    ierr = DMCreateGlobalVector(da_new, &v_new); CHKERRQ(ierr);
  }

  // Copy all the values from the old Vec to the new one:
  double ***a_new;
  double ***a_old;
  ierr = DMDAVecGetArrayDOF(m_da, v, &a_old); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da_new, v_new, &a_new); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (int j=grid->ys; j<grid->ys+grid->ym; j++) {
      for (int k=0; k < old_Mz; k++)
        a_new[i][j][k] = a_old[i][j][k];
    }
  }
  ierr = DMDAVecRestoreArrayDOF(m_da, v, &a_old); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da_new, v_new, &a_new); CHKERRQ(ierr);

  // Deallocate old Vec:
  ierr = VecDestroy(&v); CHKERRQ(ierr);
  v = v_new;
  m_da = da_new;

  // IceGrid will dispose of the old DA

  return 0;
}

//! Checks if the current IceModelVec3 has NANs and reports if it does.
/*! Up to a fixed number of messages are printed at stdout.  Returns the full
 count of NANs (which is a nonzero) on this rank. */
PetscErrorCode  IceModelVec3D::has_nan() {
  PetscErrorCode ierr;
  double *tmp;
  int retval=0;
  const int max_print_this_rank=10;

  ierr = begin_access(); CHKERRQ(ierr);
  int i = 0, j = 0;
  unsigned int k = 0;
  for (i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (j=grid->ys; j<grid->ys+grid->ym; ++j) {
      ierr = getInternalColumn(i, j, &tmp); CHKERRQ(ierr);
      for (k = 0; k < m_n_levels; k++) {
        if (gsl_isnan(tmp[k])) {
          retval++;
          if (retval <= max_print_this_rank) {
            ierr = PetscSynchronizedPrintf(grid->com, 
               "IceModelVec3 %s: NAN (or uninitialized) at i = %d, j = %d, k = %d on rank = %d\n",
               m_name.c_str(), i, j, k, grid->rank); CHKERRQ(ierr);
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
             m_name.c_str(), retval, grid->rank); CHKERRQ(ierr);
  }

#if PETSC_VERSION_LT(3,5,0)
  ierr = PetscSynchronizedFlush(grid->com); CHKERRQ(ierr);
#else
  ierr = PetscSynchronizedFlush(grid->com, PETSC_STDOUT); CHKERRQ(ierr);
#endif

  return retval;
}

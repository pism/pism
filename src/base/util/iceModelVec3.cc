// Copyright (C) 2008--2015 Ed Bueler and Constantine Khroulev
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
#include "PISMConfig.hh"

#include "error_handling.hh"

#include <cassert>

namespace pism {

// this file contains method for derived class IceModelVec3

// methods for base class IceModelVec and derived class IceModelVec2S
// are in "iceModelVec.cc"

IceModelVec3D::IceModelVec3D() : IceModelVec() {
}

IceModelVec3D::~IceModelVec3D() {
  destroy();
}

IceModelVec3::IceModelVec3() {
  // empty
}

IceModelVec3::~IceModelVec3() {
  // empty
}

//! Allocate a DA and a Vec from information in IceGrid.
void IceModelVec3D::allocate(const IceGrid &my_grid, const std::string &my_name,
                             IceModelVecKind ghostedp, const std::vector<double> &levels,
                             unsigned int stencil_width) {

  assert(m_v == NULL);
  
  m_grid = &my_grid;

  zlevels = levels;
  m_n_levels = (unsigned int)zlevels.size();
  m_da_stencil_width = stencil_width;

  m_da = m_grid->get_dm(this->m_n_levels, this->m_da_stencil_width);

  m_has_ghosts = (ghostedp == WITH_GHOSTS);

  if (m_has_ghosts == true) {
    DMCreateLocalVector(*m_da, &m_v);
  } else {
    DMCreateGlobalVector(*m_da, &m_v);
  }

  m_name = my_name;

  m_metadata.push_back(NCSpatialVariable(m_grid->config.get_unit_system(),
                                         my_name, *m_grid, zlevels));
}

bool IceModelVec3D::isLegalLevel(double z) const {
  double z_min = zlevels.front(),
    z_max = zlevels.back();
  if (z < z_min - 1.0e-6 || z > z_max + 1.0e-6) {
    return false;
  }
  return true;
}

//! Set all values of scalar quantity to given a single value in a particular column.
void IceModelVec3D::setColumn(int i, int j, double c) {
  PetscErrorCode ierr;
#if (PISM_DEBUG==1)
  assert(array != NULL);
  check_array_indices(i, j, 0);
#endif

  double ***arr = (double***) array;

  if (c == 0.0) {
    ierr = PetscMemzero(arr[i][j], m_n_levels * sizeof(double));
    PISM_PETSC_CHK(ierr, "PetscMemzero");
  } else {
    for (unsigned int k=0; k < m_n_levels; k++) {
      arr[i][j][k] = c;
    }
  }
}

//! Return value of scalar quantity at level z (m) above base of ice (by linear interpolation).
double IceModelVec3D::getValZ(int i, int j, double z) const {
#if (PISM_DEBUG==1)
  assert(array != NULL);
  check_array_indices(i, j, 0);

  if (not isLegalLevel(z)) {
    throw RuntimeError::formatted("IceModelVec3 getValZ(): level %f is not legal; name = %s",
                                  z, m_name.c_str());
  }
#endif

  double ***arr = (double***) array;
  if (z >= zlevels.back()) {
    return arr[i][j][m_n_levels - 1];
  } else if (z <= zlevels.front()) {
    return arr[i][j][0];
  }

  int mcurr = 0;
  while (zlevels[mcurr+1] < z) {
    mcurr++;
  }

  const double incr = (z - zlevels[mcurr]) / (zlevels[mcurr+1] - zlevels[mcurr]);
  const double valm = arr[i][j][mcurr];
  return valm + incr * (arr[i][j][mcurr+1] - valm);
}

//! Copies a horizontal slice at level z of an IceModelVec3 into a Vec gslice.
/*!
 * FIXME: this method is misnamed: the slice is horizontal in the PISM
 * coordinate system, not in reality.
 */
void  IceModelVec3::getHorSlice(Vec &gslice, double z) const {
  double    **slice_val;

  PISMDM::Ptr da2 = m_grid->get_dm(1, m_grid->config.get("grid_max_stencil_width"));

  IceModelVec::AccessList list(*this);
  DMDAVecGetArray(*da2, gslice, &slice_val);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    slice_val[i][j] = getValZ(i,j,z);
  }
  DMDAVecRestoreArray(*da2, gslice, &slice_val);
}

//! Copies a horizontal slice at level z of an IceModelVec3 into an IceModelVec2S gslice.
/*!
 * FIXME: this method is misnamed: the slice is horizontal in the PISM
 * coordinate system, not in reality.
 */
void  IceModelVec3::getHorSlice(IceModelVec2S &gslice, double z) const {
  IceModelVec::AccessList list(*this);
  list.add(gslice);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    gslice(i, j) = getValZ(i, j, z);
  }
}


//! Copies the values of an IceModelVec3 at the ice surface (specified by the level myH) to an IceModelVec2S gsurf.
void IceModelVec3::getSurfaceValues(IceModelVec2S &surface_values,
                                    const IceModelVec2S &H) const {
  IceModelVec::AccessList list(*this);
  list.add(surface_values);
  list.add(H);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    surface_values(i, j) = getValZ(i, j, H(i, j));
  }
}

void  IceModelVec3D::getInternalColumn(int i, int j, double **valsPTR) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  double ***arr = (double***) array;
  *valsPTR = arr[i][j];
}

void  IceModelVec3D::getInternalColumn(int i, int j, const double **valsPTR) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  double ***arr = (double***) array;
  *valsPTR = arr[i][j];
}

void  IceModelVec3D::setInternalColumn(int i, int j, double *valsIN) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  double ***arr = (double***) array;
  PetscErrorCode ierr = PetscMemcpy(arr[i][j], valsIN, m_n_levels*sizeof(double));
  PISM_PETSC_CHK(ierr, "PetscMemcpy");
}

void  IceModelVec3::create(const IceGrid &my_grid, const std::string &my_name,
                           IceModelVecKind ghostedp,
                           unsigned int stencil_width) {

  IceModelVec3D::allocate(my_grid, my_name, ghostedp,
                          my_grid.z(), stencil_width);
}

} // end of namespace pism

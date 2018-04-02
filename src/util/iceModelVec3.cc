// Copyright (C) 2008--2018 Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cstdlib>
#include <cassert>

#include <memory>
using std::dynamic_pointer_cast;

#include <petscdmda.h>

#include "iceModelVec.hh"
#include "IceGrid.hh"
#include "ConfigInterface.hh"

#include "error_handling.hh"

namespace pism {

// this file contains method for derived class IceModelVec3

// methods for base class IceModelVec and derived class IceModelVec2S
// are in "iceModelVec.cc"

IceModelVec3D::IceModelVec3D()
  : IceModelVec() {
  m_bsearch_accel = gsl_interp_accel_alloc();
  if (m_bsearch_accel == NULL) {
    throw RuntimeError(PISM_ERROR_LOCATION, "Failed to allocate a GSL interpolation accelerator");
  }
}

IceModelVec3D::~IceModelVec3D() {
  gsl_interp_accel_free(m_bsearch_accel);
}

IceModelVec3::IceModelVec3() {
  // empty
}

IceModelVec3::IceModelVec3(IceGrid::ConstPtr grid, const std::string &short_name,
                           IceModelVecKind ghostedp,
                           unsigned int stencil_width) {
  create(grid, short_name, ghostedp, stencil_width);
}

IceModelVec3::~IceModelVec3() {
  // empty
}

IceModelVec3::Ptr IceModelVec3::To3DScalar(IceModelVec::Ptr input) {
  IceModelVec3::Ptr result = dynamic_pointer_cast<IceModelVec3,IceModelVec>(input);
  if (not (bool)result) {
    throw RuntimeError(PISM_ERROR_LOCATION, "dynamic cast failure");
  }
  return result;
}

//! Allocate a DA and a Vec from information in IceGrid.
void IceModelVec3D::allocate(IceGrid::ConstPtr grid, const std::string &name,
                             IceModelVecKind ghostedp, const std::vector<double> &levels,
                             unsigned int stencil_width) {
  PetscErrorCode ierr;
  m_grid = grid;

  m_zlevels = levels;
  m_da_stencil_width = stencil_width;

  m_da = m_grid->get_dm(this->m_zlevels.size(), this->m_da_stencil_width);

  m_has_ghosts = (ghostedp == WITH_GHOSTS);

  if (m_has_ghosts == true) {
    ierr = DMCreateLocalVector(*m_da, m_v.rawptr());
    PISM_CHK(ierr, "DMCreateLocalVector");
  } else {
    ierr = DMCreateGlobalVector(*m_da, m_v.rawptr());
    PISM_CHK(ierr, "DMCreateGlobalVector");
  }

  m_name = name;

  m_metadata.push_back(SpatialVariableMetadata(m_grid->ctx()->unit_system(),
                                               name, m_zlevels));
}

bool IceModelVec3D::isLegalLevel(double z) const {
  double z_min = m_zlevels.front(),
    z_max = m_zlevels.back();
  if (z < z_min - 1.0e-6 || z > z_max + 1.0e-6) {
    return false;
  }
  return true;
}

//! Set all values of scalar quantity to given a single value in a particular column.
void IceModelVec3D::set_column(int i, int j, double c) {
  PetscErrorCode ierr;
#if (PISM_DEBUG==1)
  assert(m_array != NULL);
  check_array_indices(i, j, 0);
#endif

  double ***arr = (double***) m_array;

  if (c == 0.0) {
    ierr = PetscMemzero(arr[j][i], m_zlevels.size() * sizeof(double));
    PISM_CHK(ierr, "PetscMemzero");
  } else {
    unsigned int nlevels = m_zlevels.size();
    for (unsigned int k=0; k < nlevels; k++) {
      arr[j][i][k] = c;
    }
  }
}

void IceModelVec3D::set_column(int i, int j, const std::vector<double> &data) {
  double *column = get_column(i, j);
  for (unsigned int k = 0; k < data.size(); ++k) {
    column[k] = data[k];
  }

}

const std::vector<double> IceModelVec3D::get_column_vector(int i, int j) const {
  unsigned int n = m_zlevels.size();
  std::vector<double> result(n);
  const double *data = get_column(i, j);
  for (unsigned int k = 0; k < n; ++k) {
    result[k] = data[k];
  }
  return result;
}


//! Return value of scalar quantity at level z (m) above base of ice (by linear interpolation).
double IceModelVec3D::getValZ(int i, int j, double z) const {
#if (PISM_DEBUG==1)
  assert(m_array != NULL);
  check_array_indices(i, j, 0);

  if (not isLegalLevel(z)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "IceModelVec3 getValZ(): level %f is not legal; name = %s",
                                  z, m_name.c_str());
  }
#endif

  double ***arr = (double***) m_array;
  if (z >= m_zlevels.back()) {
    unsigned int nlevels = m_zlevels.size();
    return arr[j][i][nlevels - 1];
  } else if (z <= m_zlevels.front()) {
    return arr[j][i][0];
  }

  unsigned int mcurr = gsl_interp_accel_find(m_bsearch_accel, &m_zlevels[0], m_zlevels.size(), z);

  const double incr = (z - m_zlevels[mcurr]) / (m_zlevels[mcurr+1] - m_zlevels[mcurr]);
  const double valm = arr[j][i][mcurr];
  return valm + incr * (arr[j][i][mcurr+1] - valm);
}

//! Copies a horizontal slice at level z of an IceModelVec3 into a Vec gslice.
/*!
 * FIXME: this method is misnamed: the slice is horizontal in the PISM
 * coordinate system, not in reality.
 */
void  IceModelVec3::getHorSlice(Vec &gslice, double z) const {

  petsc::DM::Ptr da2 = m_grid->get_dm(1, m_grid->ctx()->config()->get_double("grid.max_stencil_width"));

  IceModelVec::AccessList list(*this);
  petsc::DMDAVecArray slice(da2, gslice);
  double **slice_val = (double**)slice.get();

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      slice_val[j][i] = getValZ(i,j,z);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

//! Copies a horizontal slice at level z of an IceModelVec3 into an IceModelVec2S gslice.
/*!
 * FIXME: this method is misnamed: the slice is horizontal in the PISM
 * coordinate system, not in reality.
 */
void  IceModelVec3::getHorSlice(IceModelVec2S &gslice, double z) const {
  IceModelVec::AccessList list{this, &gslice};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      gslice(i, j) = getValZ(i, j, z);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


//! Copies the values of an IceModelVec3 at the ice surface (specified by the level myH) to an IceModelVec2S gsurf.
void IceModelVec3::getSurfaceValues(IceModelVec2S &surface_values,
                                    const IceModelVec2S &H) const {
  IceModelVec::AccessList list{this, &surface_values, &H};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      surface_values(i, j) = getValZ(i, j, H(i, j));
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

double* IceModelVec3D::get_column(int i, int j) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return ((double***) m_array)[j][i];
}

const double* IceModelVec3D::get_column(int i, int j) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return ((double***) m_array)[j][i];
}

void  IceModelVec3D::set_column(int i, int j, const double *valsIN) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  double ***arr = (double***) m_array;
  PetscErrorCode ierr = PetscMemcpy(arr[j][i], valsIN, m_zlevels.size()*sizeof(double));
  PISM_CHK(ierr, "PetscMemcpy");
}

void  IceModelVec3::create(IceGrid::ConstPtr grid, const std::string &name,
                           IceModelVecKind ghostedp,
                           unsigned int stencil_width) {

  IceModelVec3D::allocate(grid, name, ghostedp,
                          grid->z(), stencil_width);
}

/** Sum a 3-D vector in the Z direction to create a 2-D vector.

Note that this sums up all the values in a column, including ones
above the ice. This may or may not be what you need. Also, take a look
at IceModel::compute_ice_enthalpy(PetscScalar &result) in iMreport.cc.

As for the difference between IceModelVec2 and IceModelVec2S, the
former can store fields with more than 1 "degree of freedom" per grid
point (such as 2D fields on the "staggered" grid, with the first
degree of freedom corresponding to the i-offset and second to
j-offset).

IceModelVec2S is just IceModelVec2 with "dof == 1", and
IceModelVec2V is IceModelVec2 with "dof == 2". (Plus some extra
methods, of course.)

Either one of IceModelVec2 and IceModelVec2S would work in this
case.

Computes output = A*output + B*sum_columns(input) + C

@see https://github.com/pism/pism/issues/229 */
void IceModelVec3::sumColumns(IceModelVec2S &output, double A, double B) const {
  const unsigned int Mz = m_grid->Mz();

  AccessList access{this, &output};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double *column = this->get_column(i,j);

    double scalar_sum = 0.0;
    for (unsigned int k = 0; k < Mz; ++k) {
      scalar_sum += column[k];
    }
    output(i,j) = A * output(i,j) + B * scalar_sum;
  }
}

} // end of namespace pism

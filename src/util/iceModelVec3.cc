// Copyright (C) 2008--2018, 2020, 2021 Ed Bueler and Constantine Khroulev
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
#include "pism/util/Context.hh"
#include "pism/util/IceModelVec_impl.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism {

// this file contains method for derived class IceModelVec3

// methods for base class IceModelVec and derived class IceModelVec2S
// are in "iceModelVec.cc"

IceModelVec3::IceModelVec3(IceGrid::ConstPtr grid, const std::string &name,
                           IceModelVecKind ghostedp, const std::vector<double> &levels,
                           unsigned int stencil_width)
  : IceModelVec() {

  {
    PetscErrorCode ierr;
    m_impl->grid = grid;

    m_impl->zlevels = levels;
    m_impl->da_stencil_width = stencil_width;

    m_impl->da = m_impl->grid->get_dm(m_impl->zlevels.size(), m_impl->da_stencil_width);

    m_impl->ghosted = (ghostedp == WITH_GHOSTS);

    if (m_impl->ghosted) {
      ierr = DMCreateLocalVector(*m_impl->da, m_impl->v.rawptr());
      PISM_CHK(ierr, "DMCreateLocalVector");
    } else {
      ierr = DMCreateGlobalVector(*m_impl->da, m_impl->v.rawptr());
      PISM_CHK(ierr, "DMCreateGlobalVector");
    }

    m_impl->name = name;

    m_impl->metadata = {SpatialVariableMetadata(m_impl->grid->ctx()->unit_system(),
                                                name, m_impl->zlevels)};

  }

  m_impl->bsearch_accel = gsl_interp_accel_alloc();
  if (m_impl->bsearch_accel == NULL) {
    throw RuntimeError(PISM_ERROR_LOCATION, "Failed to allocate a GSL interpolation accelerator");
  }

}

IceModelVec3::IceModelVec3(IceGrid::ConstPtr grid,
                           const std::string &name,
                           const std::string &z_name,
                           const std::vector<double> &z_levels,
                           const std::map<std::string, std::string> &z_attrs)
  : IceModelVec3(grid, name, WITHOUT_GHOSTS, z_levels) {

  m_impl->dof = 1;

  m_impl->metadata[0] = SpatialVariableMetadata(m_impl->grid->ctx()->unit_system(),
                                                m_impl->name, m_impl->zlevels);
  m_impl->metadata[0].get_z().set_name(z_name);

  for (auto z_attr : z_attrs) {
    m_impl->metadata[0].get_z().set_string(z_attr.first, z_attr.second);
  }
}

IceModelVec3::~IceModelVec3() {
  gsl_interp_accel_free(m_impl->bsearch_accel);
}

IceModelVec3::Ptr IceModelVec3::To3DScalar(IceModelVec::Ptr input) {
  IceModelVec3::Ptr result = dynamic_pointer_cast<IceModelVec3,IceModelVec>(input);
  if (not (bool)result) {
    throw RuntimeError(PISM_ERROR_LOCATION, "dynamic cast failure");
  }
  return result;
}

bool IceModelVec3::legal_level(double z) const {
  double z_min = m_impl->zlevels.front(),
    z_max = m_impl->zlevels.back();
  if (z < z_min - 1.0e-6 || z > z_max + 1.0e-6) {
    return false;
  }
  return true;
}

//! Set all values of scalar quantity to given a single value in a particular column.
void IceModelVec3::set_column(int i, int j, double c) {
  PetscErrorCode ierr;
#if (Pism_DEBUG==1)
  assert(m_array != NULL);
  check_array_indices(i, j, 0);
#endif

  double ***arr = (double***) m_array;

  if (c == 0.0) {
    ierr = PetscMemzero(arr[j][i], m_impl->zlevels.size() * sizeof(double));
    PISM_CHK(ierr, "PetscMemzero");
  } else {
    unsigned int nlevels = m_impl->zlevels.size();
    for (unsigned int k=0; k < nlevels; k++) {
      arr[j][i][k] = c;
    }
  }
}

void IceModelVec3::set_column(int i, int j, const std::vector<double> &data) {
  double *column = get_column(i, j);
  for (unsigned int k = 0; k < data.size(); ++k) {
    column[k] = data[k];
  }

}

const std::vector<double> IceModelVec3::get_column_vector(int i, int j) const {
  unsigned int n = m_impl->zlevels.size();
  std::vector<double> result(n);
  const double *data = get_column(i, j);
  for (unsigned int k = 0; k < n; ++k) {
    result[k] = data[k];
  }
  return result;
}


//! Return value of scalar quantity at level z (m) above base of ice (by linear interpolation).
double IceModelVec3::interpolate(int i, int j, double z) const {
#if (Pism_DEBUG==1)
  assert(m_array != NULL);
  check_array_indices(i, j, 0);

  if (not legal_level(z)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "IceModelVec3 interpolate(): level %f is not legal; name = %s",
                                  z, m_impl->name.c_str());
  }
#endif

  double ***arr = (double***) m_array;
  if (z >= m_impl->zlevels.back()) {
    unsigned int nlevels = m_impl->zlevels.size();
    return arr[j][i][nlevels - 1];
  } else if (z <= m_impl->zlevels.front()) {
    return arr[j][i][0];
  }

  unsigned int mcurr = gsl_interp_accel_find(m_impl->bsearch_accel,
                                             &m_impl->zlevels[0], m_impl->zlevels.size(), z);

  const double incr = (z - m_impl->zlevels[mcurr]) / (m_impl->zlevels[mcurr+1] - m_impl->zlevels[mcurr]);
  const double valm = arr[j][i][mcurr];
  return valm + incr * (arr[j][i][mcurr+1] - valm);
}

//! Copies a horizontal slice at level z of an IceModelVec3 into an IceModelVec2S gslice.
void IceModelVec3::extract_surface(double z, IceModelVec2S &output) const {
  IceModelVec::AccessList list{this, &output};

  ParallelSection loop(m_impl->grid->com);
  try {
    for (Points p(*m_impl->grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      output(i, j) = interpolate(i, j, z);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


//! Copies the values of an IceModelVec3 at the ice surface (specified by the level myH) to an IceModelVec2S gsurf.
void IceModelVec3::extract_surface(const IceModelVec2S &H, IceModelVec2S &output) const {
  IceModelVec::AccessList list{this, &output, &H};

  ParallelSection loop(m_impl->grid->com);
  try {
    for (Points p(*m_impl->grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      output(i, j) = interpolate(i, j, H(i, j));
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

double* IceModelVec3::get_column(int i, int j) {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return ((double***) m_array)[j][i];
}

const double* IceModelVec3::get_column(int i, int j) const {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return ((double***) m_array)[j][i];
}

void  IceModelVec3::set_column(int i, int j, const double *valsIN) {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  double ***arr = (double***) m_array;
  PetscErrorCode ierr = PetscMemcpy(arr[j][i], valsIN, m_impl->zlevels.size()*sizeof(double));
  PISM_CHK(ierr, "PetscMemcpy");
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
void IceModelVec3::sum_columns(double A, double B, IceModelVec2S &output) const {
  const unsigned int Mz = m_impl->grid->Mz();

  AccessList access{this, &output};

  for (Points p(*m_impl->grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double *column = this->get_column(i,j);

    double scalar_sum = 0.0;
    for (unsigned int k = 0; k < Mz; ++k) {
      scalar_sum += column[k];
    }
    output(i,j) = A * output(i,j) + B * scalar_sum;
  }
}

void IceModelVec3::copy_from(const IceModelVec3 &input) {
  IceModelVec::AccessList list {this, &input};

  assert(levels().size() == input.levels().size());

  ParallelSection loop(m_impl->grid->com);
  try {
    for (Points p(*m_impl->grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      this->set_column(i, j, input.get_column(i, j));
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  update_ghosts();

  inc_state_counter();
}


} // end of namespace pism

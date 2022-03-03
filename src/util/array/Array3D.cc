// Copyright (C) 2008--2018, 2020, 2021, 2022 Ed Bueler and Constantine Khroulev
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

#include "pism/util/array/Array3D.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/array/Array_impl.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace array {

// this file contains method for derived class array::Array3D

// methods for base class IceModelVec and derived class IceModelVec2S
// are in "iceModelVec.cc"

Array3D::Array3D(IceGrid::ConstPtr grid,
                 const std::string &name,
                 IceModelVecKind ghostedp,
                 const std::vector<double> &levels,
                 unsigned int stencil_width)
  : Array(grid, name, ghostedp, 1, stencil_width, levels) {
  set_begin_access_use_dof(true);
}

//! Set all values of scalar quantity to given a single value in a particular column.
void Array3D::set_column(int i, int j, double c) {
  PetscErrorCode ierr;
#if (Pism_DEBUG==1)
  assert(m_array != NULL);
  check_array_indices(i, j, 0);
#endif

  double ***arr = (double***) m_array;

  if (c == 0.0) {
    ierr = PetscMemzero(arr[j][i], levels().size() * sizeof(double));
    PISM_CHK(ierr, "PetscMemzero");
  } else {
    unsigned int nlevels = levels().size();
    for (unsigned int k=0; k < nlevels; k++) {
      arr[j][i][k] = c;
    }
  }
}

void  Array3D::set_column(int i, int j, const double *input) {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  double ***arr = (double***) m_array;
  PetscErrorCode ierr = PetscMemcpy(arr[j][i], input, m_impl->zlevels.size()*sizeof(double));
  PISM_CHK(ierr, "PetscMemcpy");
}

#if (Pism_DEBUG==1)
static bool legal_level(const std::vector<double> &levels, double z) {
  double z_min = levels.front();
  double z_max = levels.back();
  const double eps = 1.0e-6;
  return not (z < z_min - eps || z > z_max + eps);
}
#endif

//! Return value of scalar quantity at level z (m) above base of ice (by linear interpolation).
double Array3D::interpolate(int i, int j, double z) const {
  auto zs = levels();

#if (Pism_DEBUG==1)
  assert(m_array != NULL);
  check_array_indices(i, j, 0);

  if (not legal_level(zs, z)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Array3D interpolate(): level %f is not legal; name = %s",
                                  z, m_impl->name.c_str());
  }
#endif

  const auto* column = get_column(i, j);

  if (z >= zs.back()) {
    unsigned int nlevels = zs.size();
    return column[nlevels - 1];
  }

  if (z <= zs.front()) {
    return column[0];
  }

  auto mcurr = gsl_interp_accel_find(m_impl->bsearch_accel,
                                     &zs[0], zs.size(), z);

  const double incr = (z - zs[mcurr]) / (zs[mcurr+1] - zs[mcurr]);
  const double valm = column[mcurr];
  return valm + incr * (column[mcurr + 1] - valm);
}

double* Array3D::get_column(int i, int j) {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return ((double***) m_array)[j][i];
}

const double* Array3D::get_column(int i, int j) const {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return ((double***) m_array)[j][i];
}

//! Copies a horizontal slice at level z of an Array3D into an Scalar gslice.
void extract_surface(const Array3D &data, double z, Scalar &output) {
  array::AccessScope list{&data, &output};

  ParallelSection loop(output.grid()->com);
  try {
    for (Points p(*output.grid()); p; p.next()) {
      const int i = p.i(), j = p.j();
      output(i, j) = data.interpolate(i, j, z);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


//! Copies the values of an Array3D at the ice surface (specified by the level `z`) to an Scalar gsurf.
void extract_surface(const Array3D &data, const Scalar &z, Scalar &output) {
  array::AccessScope list{&data, &output, &z};

  ParallelSection loop(output.grid()->com);
  try {
    for (Points p(*output.grid()); p; p.next()) {
      const int i = p.i(), j = p.j();
      output(i, j) = data.interpolate(i, j, z(i, j));
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

/** Sum a 3-D vector in the Z direction to create a 2-D vector.

Note that this sums up all the values in a column, including ones
above the ice. This may or may not be what you need. Also, take a look
at IceModel::compute_ice_enthalpy(PetscScalar &result) in iMreport.cc.

As for the difference between IceModelVec2 and array::Scalar, the
former can store fields with more than 1 "degree of freedom" per grid
point (such as 2D fields on the "staggered" grid, with the first
degree of freedom corresponding to the i-offset and second to
j-offset).

array::Scalar is just IceModelVec2 with "dof == 1", and
IceModelVec2V is IceModelVec2 with "dof == 2". (Plus some extra
methods, of course.)

Either one of IceModelVec2 and array::Scalar would work in this
case.

Computes output = A*output + B*sum_columns(input) + C

@see https://github.com/pism/pism/issues/229 */
void sum_columns(const Array3D &data, double A, double B, Scalar &output) {
  const unsigned int Mz = data.grid()->Mz();

  AccessScope access{&data, &output};

  for (Points p(*data.grid()); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double *column = data.get_column(i,j);

    double scalar_sum = 0.0;
    for (unsigned int k = 0; k < Mz; ++k) {
      scalar_sum += column[k];
    }
    output(i,j) = A * output(i,j) + B * scalar_sum;
  }
}

void Array3D::copy_from(const Array3D &input) {
  array::AccessScope list {this, &input};

  assert(levels().size() == input.levels().size());
  assert(ndof() == input.ndof());

  // 3D arrays have more than one level and ndof() of 1, collections of fields have one
  // level and ndof() > 1
  auto N = std::max((size_t)ndof(), levels().size());

  ParallelSection loop(m_impl->grid->com);
  try {
    for (Points p(*m_impl->grid); p; p.next()) {
      const int i = p.i(), j = p.j();

#if PETSC_VERSION_LT(3, 12, 0)
      PetscMemmove(this->get_column(i, j),
                   const_cast<double*>(input.get_column(i, j)), N*sizeof(double));
#else
      PetscArraymove(this->get_column(i, j), input.get_column(i, j), N);
#endif
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  update_ghosts();

  inc_state_counter();
}

std::shared_ptr<Array3D> Array3D::duplicate() const {

  auto result = std::make_shared<Array3D>(this->grid(),
                                          this->get_name(),
                                          WITHOUT_GHOSTS,
                                          this->levels());

  result->metadata() = this->metadata();

  return result;
}

} // end of namespace array
} // end of namespace pism

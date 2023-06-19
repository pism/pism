/* Copyright (C) 2015, 2016, 2017, 2023 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace petsc {

// Wrapper around Vec (calls VecDestroy)

Vec::Vec() {
  m_value = NULL;
}

Vec::Vec(::Vec v) {
  m_value = v;
}

Vec::~Vec() {
  if (m_value != NULL) {
    PetscErrorCode ierr = VecDestroy(&m_value); CHKERRCONTINUE(ierr);
  }
}

// Wrapper around VecGetArray / VecRestoreArray

VecArray::VecArray(::Vec v)
  : m_v(v), m_array(NULL) {
  PetscErrorCode ierr = VecGetArray(m_v, &m_array);
  PISM_CHK(ierr, "VecGetArray");
}

VecArray::~VecArray() {
  PetscErrorCode ierr = VecRestoreArray(m_v, &m_array); CHKERRCONTINUE(ierr);
}

double* VecArray::get() {
  return m_array;
}

// Wrapper around VecGetArray2d / VecRestoreArray2d

VecArray2D::VecArray2D(::Vec vec, int Mx, int My)
    : m_Mx(Mx), m_My(My), m_i_offset(0), m_j_offset(0), m_v(vec) {
  PetscErrorCode ierr = VecGetArray2d(m_v, m_My, m_Mx, 0, 0, &m_array);
  PISM_CHK(ierr, "VecGetArray2d");
}

VecArray2D::VecArray2D(::Vec vec, int Mx, int My, int i0, int j0)
  : m_Mx(Mx), m_My(My), m_i_offset(i0), m_j_offset(j0), m_v(vec) {
  PetscErrorCode ierr = VecGetArray2d(m_v, m_My, m_Mx, 0, 0, &m_array);
  PISM_CHK(ierr, "VecGetArray2d");
}

VecArray2D::~VecArray2D() {
  PetscErrorCode ierr = VecRestoreArray2d(m_v, m_My, m_Mx, 0, 0, &m_array); CHKERRCONTINUE(ierr);
}

// Wrapper around DMDAVecGetArray / DMDAVecRestoreArray

DMDAVecArray::DMDAVecArray(DM::Ptr dm, ::Vec v)
  : m_dm(dm), m_v(v) {
  PetscErrorCode ierr = DMDAVecGetArray(*m_dm, m_v, &m_array);
  PISM_CHK(ierr, "DMDAVecGetArray");
}

DMDAVecArray::~DMDAVecArray() {
  PetscErrorCode ierr = DMDAVecRestoreArray(*m_dm, m_v, &m_array); CHKERRCONTINUE(ierr);
}

void* DMDAVecArray::get() {
  return m_array;
}

// Wrapper around DMDAVecGetArrayDOF / DMDAVecRestoreArrayDOF

DMDAVecArrayDOF::DMDAVecArrayDOF(DM::Ptr dm, ::Vec v)
  : m_dm(dm), m_v(v) {
  PetscErrorCode ierr = DMDAVecGetArrayDOF(*m_dm, m_v, &m_array);
  PISM_CHK(ierr, "DMDAVecGetArrayDOF");
}

DMDAVecArrayDOF::~DMDAVecArrayDOF() {
  PetscErrorCode ierr = DMDAVecRestoreArrayDOF(*m_dm, m_v, &m_array); CHKERRCONTINUE(ierr);
}

void* DMDAVecArrayDOF::get() {
  return m_array;
}

// Wrapper around DMGetGlobalVector / DMRestoreGlobalVector

TemporaryGlobalVec::TemporaryGlobalVec(DM::Ptr dm) {
  m_dm = dm;
  PetscErrorCode ierr = DMGetGlobalVector(*m_dm, &m_value);
  PISM_CHK(ierr, "DMGetGlobalVector");
}

TemporaryGlobalVec::~TemporaryGlobalVec() {
  // This takes advantage of the fact that the destructor of a derived
  // class is called before the destructor of its base class, so we
  // can set m_value to NULL and turn the destructor of the base class
  // (Vec) into a no-op.
  if (m_value != NULL) {
    PetscErrorCode ierr = DMRestoreGlobalVector(*m_dm, &m_value); CHKERRCONTINUE(ierr);
    m_value = NULL;
  }
}


} // end of namespace petsc
} // end of namespace pism

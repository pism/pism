// Copyright (C) 2012, 2014, 2015, 2017, 2023  David Maxwell and Constantine Khroulev
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

#include <cassert>

#include "pism/inverse/IPTwoBlockVec.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace inverse {

IPTwoBlockVec::IPTwoBlockVec(Vec a, Vec b) {
  PetscErrorCode ierr;

  MPI_Comm comm, comm_b;
  ierr = PetscObjectGetComm((PetscObject)a, &comm); PISM_CHK(ierr, "PetscObjectGetComm");
  ierr = PetscObjectGetComm((PetscObject)b, &comm_b); PISM_CHK(ierr, "PetscObjectGetComm");
  assert(comm == comm_b);

  // These have to be PetscInt because of PETSc calls below
  PetscInt lo_a, hi_a;
  ierr = VecGetOwnershipRange(a, &lo_a, &hi_a); PISM_CHK(ierr, "VecGetOwnershipRange");
  ierr = VecGetSize(a, &m_na_global); PISM_CHK(ierr, "VecGetSize");
  m_na_local = hi_a - lo_a;

  // These have to be PetscInt because of PETSc calls below
  PetscInt lo_b, hi_b;
  ierr = VecGetOwnershipRange(b, &lo_b, &hi_b); PISM_CHK(ierr, "VecGetOwnershipRange");
  ierr = VecGetSize(b, &m_nb_global); PISM_CHK(ierr, "VecGetSize");
  m_nb_local = hi_b - lo_b;

  petsc::IS is_a, is_b;
  // a in a
  ierr = ISCreateStride(comm, m_na_local, lo_a, 1,
                        is_a.rawptr()); PISM_CHK(ierr, "ISCreateStride");
  // a in ab
  ierr = ISCreateStride(comm, m_na_local, lo_a + lo_b, 1,
                        m_a_in_ab.rawptr()); PISM_CHK(ierr, "ISCreateStride");

  // b in b
  ierr = ISCreateStride(comm, m_nb_local, lo_b, 1,
                        is_b.rawptr()); PISM_CHK(ierr, "ISCreateStride");
  // b in ab
  ierr = ISCreateStride(comm, m_nb_local, lo_a + lo_b + m_na_local, 1,
                        m_b_in_ab.rawptr()); PISM_CHK(ierr, "ISCreateStride");

  ierr = VecCreate(comm, m_ab.rawptr()); PISM_CHK(ierr, "VecCreate");

  ierr = VecSetType(m_ab, "mpi"); PISM_CHK(ierr, "VecSetType");
  ierr = VecSetSizes(m_ab, m_na_local + m_nb_local, m_na_global + m_nb_global);
  PISM_CHK(ierr, "VecSetSizes");

  ierr = VecScatterCreate(m_ab, m_a_in_ab, a, is_a,
                          m_scatter_a.rawptr()); PISM_CHK(ierr, "VecScatterCreate");
  ierr = VecScatterCreate(m_ab, m_b_in_ab, b, is_b,
                          m_scatter_b.rawptr()); PISM_CHK(ierr, "VecScatterCreate");
}

IS IPTwoBlockVec::blockAIndexSet() {
  return m_a_in_ab;
}

IS IPTwoBlockVec::blockBIndexSet() {
  return m_b_in_ab;
}

void IPTwoBlockVec::scatter(Vec a, Vec b) {
  this->scatterToA(m_ab,a);
  this->scatterToB(m_ab,b);
}

void IPTwoBlockVec::scatterToA(Vec a) {
  this->scatterToA(m_ab,a);
}

void IPTwoBlockVec::scatterToB(Vec b) {
  this->scatterToB(m_ab,b);
}

void IPTwoBlockVec::scatter(Vec ab, Vec a, Vec b) {
  this->scatterToA(ab,a);
  this->scatterToB(ab,b);  
}

void IPTwoBlockVec::scatter_begin_end(VecScatter s, Vec a, Vec b, ScatterMode m) {
  PetscErrorCode ierr;
  ierr = VecScatterBegin(s, a, b, INSERT_VALUES, m);
  PISM_CHK(ierr, "VecScatterBegin");

  ierr = VecScatterEnd(s, a, b, INSERT_VALUES, m);
  PISM_CHK(ierr, "VecScatterEnd");
}

void IPTwoBlockVec::scatterToA(Vec ab, Vec a) {
  scatter_begin_end(m_scatter_a, ab, a, SCATTER_FORWARD);
}

void IPTwoBlockVec::scatterToB(Vec ab, Vec b) {
  scatter_begin_end(m_scatter_b, ab, b, SCATTER_FORWARD);
}

void IPTwoBlockVec::gather(Vec a, Vec b) {
  this->gatherFromA(a,m_ab);
  this->gatherFromB(b,m_ab);
}

void IPTwoBlockVec::gatherFromA(Vec a) {
  this->gatherFromA(a,m_ab);
}

void IPTwoBlockVec::gatherFromB(Vec b) {
  this->gatherFromA(b,m_ab);
}

void IPTwoBlockVec::gather(Vec a, Vec b, Vec ab) {
  this->gatherFromA(a,ab);
  this->gatherFromB(b,ab);
}

void IPTwoBlockVec::gatherFromA(Vec a, Vec ab) {
  scatter_begin_end(m_scatter_a, a, ab, SCATTER_REVERSE);
}

void IPTwoBlockVec::gatherFromB(Vec b, Vec ab) {
  scatter_begin_end(m_scatter_b, b, ab, SCATTER_REVERSE);
}

} // end of namespace inverse
} // end of namespace pism

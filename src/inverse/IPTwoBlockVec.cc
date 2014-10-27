// Copyright (C) 2012, 2014  David Maxwell
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

#include "IPTwoBlockVec.hh"
#include <assert.h>
#include "error_handling.hh"

namespace pism {

IPTwoBlockVec::IPTwoBlockVec(Vec a, Vec b) {
  PetscErrorCode ierr;
  ierr = this->construct(a,b);
  assert(ierr==0);
}

IPTwoBlockVec::~IPTwoBlockVec() {
  PetscErrorCode ierr;
  ierr = this->destruct();
  assert(ierr==0);
}

PetscErrorCode IPTwoBlockVec::construct(Vec a, Vec b)  {
  PetscErrorCode ierr;
  
  MPI_Comm comm, comm_b;
  ierr = PetscObjectGetComm((PetscObject)a,&comm);
  PISM_PETSC_CHK(ierr, "PetscObjectGetComm");
  ierr = PetscObjectGetComm((PetscObject)a,&comm_b);
  PISM_PETSC_CHK(ierr, "PetscObjectGetComm");
  assert(comm==comm_b);
  
  PetscInt lo_a, hi_a;
  ierr = VecGetOwnershipRange(a,&lo_a,&hi_a);
  PISM_PETSC_CHK(ierr, "VecGetOwnershipRange");
  ierr = VecGetSize(a,&m_na_global);
  PISM_PETSC_CHK(ierr, "VecGetSize");
  m_na_local = hi_a-lo_a;

  PetscInt lo_b, hi_b;
  ierr = VecGetOwnershipRange(b,&lo_b,&hi_b);
  PISM_PETSC_CHK(ierr, "VecGetOwnershipRange");
  ierr = VecGetSize(b,&m_nb_global);
  PISM_PETSC_CHK(ierr, "VecGetSize");
  m_nb_local = hi_b-lo_b;

  IS is_a, is_b;
  ierr = ISCreateStride(comm,m_na_local,lo_a,1,&is_a); CHKERRQ(ierr);  // a in a
  ierr = ISCreateStride(comm,m_na_local,lo_a+lo_b,1,&m_a_in_ab); CHKERRQ(ierr); // a in ab
  ierr = ISCreateStride(comm,m_nb_local,lo_b,1,&is_b); CHKERRQ(ierr); // b in b
  ierr = ISCreateStride(comm,m_nb_local,lo_a+lo_b+m_na_local,1,&m_b_in_ab); CHKERRQ(ierr);  // b in ab

  ierr = VecCreate(comm,&m_ab);
  PISM_PETSC_CHK(ierr, "VecCreate");
  ierr = VecSetType(m_ab,"mpi");
  PISM_PETSC_CHK(ierr, "VecSetType");
  ierr = VecSetSizes(m_ab,m_na_local+m_nb_local,m_na_global+m_nb_global);
  PISM_PETSC_CHK(ierr, "VecSetSizes");

  ierr = VecScatterCreate(m_ab,m_a_in_ab,a,is_a,&m_scatter_a);
  PISM_PETSC_CHK(ierr, "VecScatterCreate");
  ierr = VecScatterCreate(m_ab,m_b_in_ab,b,is_b,&m_scatter_b);
  PISM_PETSC_CHK(ierr, "VecScatterCreate");
  
  ierr = ISDestroy(&is_a); CHKERRQ(ierr);
  ierr = ISDestroy(&is_b); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IPTwoBlockVec::destruct() {
  PetscErrorCode ierr;

  ierr = VecDestroy(&m_ab);
  PISM_PETSC_CHK(ierr, "VecDestroy");
  ierr = ISDestroy(&m_a_in_ab); CHKERRQ(ierr);
  ierr = ISDestroy(&m_b_in_ab); CHKERRQ(ierr);
  ierr = VecScatterDestroy(&m_scatter_a);
  PISM_PETSC_CHK(ierr, "VecScatterDestroy");
  ierr = VecScatterDestroy(&m_scatter_b);
  PISM_PETSC_CHK(ierr, "VecScatterDestroy");
  return 0;
}

IS IPTwoBlockVec::blockAIndexSet() {
  return m_a_in_ab;
}

IS IPTwoBlockVec::blockBIndexSet() {
  return m_b_in_ab;
}

PetscErrorCode IPTwoBlockVec::scatter(Vec a, Vec b) {
  PetscErrorCode ierr;
  ierr = this->scatterToA(m_ab,a); CHKERRQ(ierr);
  ierr = this->scatterToB(m_ab,b); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode IPTwoBlockVec::scatterToA(Vec a) {
  PetscErrorCode ierr;
  ierr = this->scatterToA(m_ab,a); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode IPTwoBlockVec::scatterToB(Vec b) {
  PetscErrorCode ierr;
  ierr = this->scatterToB(m_ab,b); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IPTwoBlockVec::scatter(Vec ab, Vec a, Vec b) {
  this->scatterToA(ab,a);
  this->scatterToB(ab,b);  
  return 0;
}
PetscErrorCode IPTwoBlockVec::scatterToA(Vec ab, Vec a) {
  PetscErrorCode ierr;
  ierr = VecScatterBegin(m_scatter_a, ab, a, INSERT_VALUES, SCATTER_FORWARD);
  PISM_PETSC_CHK(ierr, "VecScatterBegin");
  ierr = VecScatterEnd(m_scatter_a, ab, a, INSERT_VALUES, SCATTER_FORWARD);
  PISM_PETSC_CHK(ierr, "VecScatterEnd");
  return 0;
}
PetscErrorCode IPTwoBlockVec::scatterToB(Vec ab, Vec b) {
  PetscErrorCode ierr;
  ierr = VecScatterBegin(m_scatter_b, ab, b, INSERT_VALUES, SCATTER_FORWARD);
  PISM_PETSC_CHK(ierr, "VecScatterBegin");
  ierr = VecScatterEnd(m_scatter_b, ab, b, INSERT_VALUES, SCATTER_FORWARD);
  PISM_PETSC_CHK(ierr, "VecScatterEnd");
  return 0;
}

PetscErrorCode IPTwoBlockVec::gather(Vec a, Vec b) {
  PetscErrorCode ierr;
  ierr = this->gatherFromA(a,m_ab); CHKERRQ(ierr);
  ierr = this->gatherFromB(b,m_ab); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode IPTwoBlockVec::gatherFromA(Vec a) {
  PetscErrorCode ierr;
  ierr = this->gatherFromA(a,m_ab); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode IPTwoBlockVec::gatherFromB(Vec b) {
  PetscErrorCode ierr;
  ierr = this->gatherFromA(b,m_ab); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IPTwoBlockVec::gather(Vec a, Vec b, Vec ab) {
  PetscErrorCode ierr;
  ierr = this->gatherFromA(a,ab); CHKERRQ(ierr);
  ierr = this->gatherFromB(b,ab); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode IPTwoBlockVec::gatherFromA(Vec a, Vec ab) {
  PetscErrorCode ierr;
  
  ierr = VecScatterBegin(m_scatter_a, a, ab, INSERT_VALUES, SCATTER_REVERSE);
  PISM_PETSC_CHK(ierr, "VecScatterBegin");
  ierr = VecScatterEnd(m_scatter_a, a, ab, INSERT_VALUES, SCATTER_REVERSE);
  PISM_PETSC_CHK(ierr, "VecScatterEnd");
  return 0;
}
PetscErrorCode IPTwoBlockVec::gatherFromB(Vec b, Vec ab) {
  PetscErrorCode ierr;
  ierr = VecScatterBegin(m_scatter_b, b, ab, INSERT_VALUES, SCATTER_REVERSE);
  PISM_PETSC_CHK(ierr, "VecScatterBegin");
  ierr = VecScatterEnd(m_scatter_b, b, ab, INSERT_VALUES, SCATTER_REVERSE);
  PISM_PETSC_CHK(ierr, "VecScatterEnd");
  return 0;
}

} // end of namespace pism

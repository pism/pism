/* Copyright (C) 2015 PISM Authors
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

#include "Vec.hh"
#include "error_handling.hh"

namespace pism {
namespace petsc {

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

TemporaryGlobalVec::TemporaryGlobalVec(DM::Ptr dm) {
  m_dm = dm;
  PetscErrorCode ierr = DMGetGlobalVector(*m_dm, &m_value);
  PISM_PETSC_CHK(ierr, "DMGetGlobalVector");
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

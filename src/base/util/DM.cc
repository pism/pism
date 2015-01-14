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

#include "DM.hh"

namespace pism {

PISMDM::PISMDM(DM dm) {
  m_dm = dm;
}

PISMDM::~PISMDM() {
  PetscErrorCode ierr = DMDestroy(&m_dm); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    // We can't do anything about this failure. We can't recover
    // from it, and it is almost certainly caused by a programming
    // error. So, we call abort().
    abort();
  }
}

DM PISMDM::get() const {
  return m_dm;
}

PISMDM::operator DM() const {
  return m_dm;
}

} // end of namespace pism

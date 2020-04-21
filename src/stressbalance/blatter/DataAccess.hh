/* Copyright (C) 2020 PISM Authors
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

#ifndef PISM_DATAACCESS_H
#define PISM_DATAACCESS_H

#include <cassert>
#include <petscdmda.h>

#include "pism/util/error_handling.hh"

namespace pism {

enum AccessType {GHOSTED, NOT_GHOSTED};
/*!
 * This template class manages access to 2D and 3D Vecs stored in a DM using
 * `PetscObjectCompose`. Performs cleanup at the end of scope.
 *
 * @param[in] da SNES DM for the solution containing 2D and 3D DMs and Vecs
 * @param[in] dim number of dimensions (2 or 3)
 * @param[in] type NOT_GHOSTED -- for setting parameters; GHOSTED -- for accessing ghosts
 *                 during residual and Jacobian evaluation
 */
template<typename T>
class DataAccess {
public:
  DataAccess(DM da, int dim, AccessType type)
    : m_local(type == GHOSTED) {
    int ierr;

    assert(dim == 2 or dim == 3);

    if (dim == 2) {
      ierr = setup(da, "2D_DM", "2D_DM_data");
    } else {
      ierr = setup(da, "3D_DM", "3D_DM_data");
    }

    if (ierr != 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "Failed to create an DataAccess instance");
    }
  }

  int setup(DM da, const char *dm_name, const char *vec_name) {
    PetscErrorCode ierr;

    m_com = MPI_COMM_SELF;
    ierr = PetscObjectGetComm((PetscObject)da, &m_com); CHKERRQ(ierr);

    ierr = PetscObjectQuery((PetscObject)da, dm_name, (PetscObject*)&m_da); CHKERRQ(ierr);

    if (!m_da) {
      SETERRQ(m_com, 1, "Failed to get the inner DM");
    }

    Vec X;
    ierr = PetscObjectQuery((PetscObject)da, vec_name, (PetscObject*)&X); CHKERRQ(ierr);

    if (!X) {
      SETERRQ(m_com, 1, "Failed to get the inner Vec");
    }

    if (m_local) {
      ierr = DMGetLocalVector(m_da, &m_x); CHKERRQ(ierr);

      ierr = DMGlobalToLocalBegin(m_da, X, INSERT_VALUES, m_x); CHKERRQ(ierr);

      ierr = DMGlobalToLocalEnd(m_da, X, INSERT_VALUES, m_x); CHKERRQ(ierr);
    } else {
      m_x = X;
    }

    ierr = DMDAVecGetArray(m_da, m_x, &m_a); CHKERRQ(ierr);

    return 0;
  }

  ~DataAccess() {
    try {
      PetscErrorCode ierr = DMDAVecRestoreArray(m_da, m_x, &m_a);
      PISM_CHK(ierr, "DMDAVecRestoreArray");

      if (m_local) {
        ierr = DMRestoreLocalVector(m_da, &m_x);
        PISM_CHK(ierr, "DMRestoreLocalVector");
      }
    } catch (...) {
      handle_fatal_errors(m_com);
    }
  }

  operator T() {
    return m_a;
  }
private:
  MPI_Comm m_com;
  bool m_local;
  DM m_da;
  Vec m_x;
  T m_a;
};

} // end of namespace pism

#endif /* PISM_DATAACCESS_H */

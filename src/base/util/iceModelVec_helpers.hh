// Copyright (C) 2011, 2013 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef _ICEMODELVEC_HELPERS_H_
#define _ICEMODELVEC_HELPERS_H_

void compute_params(IceModelVec* const x, IceModelVec* const y,
                    IceModelVec* const z, int &ghosts, bool &scatter);

//! \brief Computes result = x + alpha * y, where x, y, and z are 2D
//! IceModelVecs (scalar or vector).
/*!
 * This implementation tries to be smart about handling IceModelVecs with and
 * without ghosts and with different stencil widths.
 *
 * This template function was written to re-use this code for both
 * IceModelVec2S and IceModel2V.
 *
 * This cannot go into a protected member IceModelVec because
 * IceModelVec2S::operator() and IceModelVec2V::operator() return different
 * types.
 *
 * Note: this code uses overloaded operators (PISMVector2::operator*, etc).
 */
template<class V>
PetscErrorCode add_2d(IceModelVec* const x_in, PetscScalar alpha, IceModelVec* const y_in,
                      IceModelVec* const result) {
  PetscErrorCode ierr;

  V *x = dynamic_cast<V*>(x_in),
    *y = dynamic_cast<V*>(y_in),
    *z = dynamic_cast<V*>(result);

  if (x == NULL || y == NULL || z == NULL) {
    SETERRQ(PETSC_COMM_SELF, 1, "incompatible arguments");
  }

  int ghosts = 0;
  bool scatter = false;
  compute_params(x, y, z, ghosts, scatter);

  IceGrid *grid = z->get_grid();

  ierr = x->begin_access(); CHKERRQ(ierr);
  ierr = y->begin_access(); CHKERRQ(ierr);
  ierr = z->begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid->xs - ghosts; i < grid->xs+grid->xm + ghosts; ++i) {
    for (PetscInt j = grid->ys - ghosts; j < grid->ys+grid->ym + ghosts; ++j) {
      (*z)(i, j) = (*x)(i, j) + alpha * (*y)(i, j);
    }
  }
  ierr = z->end_access(); CHKERRQ(ierr);
  ierr = y->end_access(); CHKERRQ(ierr);
  ierr = x->end_access(); CHKERRQ(ierr);

  if (scatter) {
    ierr = z->update_ghosts(); CHKERRQ(ierr);
  }

  return 0;
}

template<class V>
PetscErrorCode copy_2d(IceModelVec* const source,
                       IceModelVec* const destination) {
  PetscErrorCode ierr;

  V *x = dynamic_cast<V*>(source),
    *z = dynamic_cast<V*>(destination);

  if (x == NULL || z == NULL) {
    SETERRQ(PETSC_COMM_SELF, 1, "incompatible arguments");
  }

  int ghosts = 0;
  bool scatter = false;
  compute_params(x, x, z, ghosts, scatter);

  IceGrid *grid = z->get_grid();

  ierr = x->begin_access(); CHKERRQ(ierr);
  ierr = z->begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid->xs - ghosts; i < grid->xs+grid->xm + ghosts; ++i) {
    for (PetscInt j = grid->ys - ghosts; j < grid->ys+grid->ym + ghosts; ++j) {
      (*z)(i, j) = (*x)(i, j);
    }
  }
  ierr = z->end_access(); CHKERRQ(ierr);
  ierr = x->end_access(); CHKERRQ(ierr);

  if (scatter) {
    ierr = z->update_ghosts(); CHKERRQ(ierr);
  }

  return 0;
}

#endif /* _ICEMODELVEC_HELPERS_H_ */

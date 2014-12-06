/* Copyright (C) 2014 PISM Authors
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

#include "Viewer.hh"

#include <petscdraw.h>
#include <cassert>
#include "error_handling.hh"

namespace pism {

Viewer::Viewer(MPI_Comm com,  const std::string &title, unsigned int target_size,
               double Lx, double Ly) {
  PetscErrorCode ierr;
  unsigned int X, Y;

  compute_size(target_size, Lx, Ly, X, Y);

  // note we reverse x <-> y; see IceGrid::allocate() for original reversal
  ierr = PetscViewerDrawOpen(com, NULL, title.c_str(),
                             PETSC_DECIDE, PETSC_DECIDE, Y, X, &m_viewer);
  PISM_PETSC_CHK(ierr, "PetscViewerDrawOpen");

  // following should be redundant, but may put up a title even under 2.3.3-p1:3 where
  // there is a no titles bug
  PetscDraw draw;
  ierr = PetscViewerDrawGetDraw(m_viewer, 0, &draw);
  PISM_PETSC_CHK(ierr, "PetscViewerDrawGetDraw");
  ierr = PetscDrawSetTitle(draw, title.c_str());
  PISM_PETSC_CHK(ierr, "PetscDrawSetTitle");
}

Viewer::~Viewer() {
  PetscViewerDestroy(&m_viewer);
}

PetscViewer Viewer::operator*() {
  return m_viewer;
}

void Viewer::compute_size(unsigned int target_size, double Lx, double Ly, unsigned int &X, unsigned int &Y) {

  assert(Lx > 0 && Ly > 0);

  // aim for smaller dimension equal to target, larger dimension larger by Ly/Lx or Lx/Ly proportion
  const double yTOx = Ly / Lx;
  if (Ly > Lx) {
    X = target_size;
    Y = (unsigned int) ((double)target_size * yTOx);
  } else {
    Y = target_size;
    X = (unsigned int) ((double)target_size / yTOx);
  }

  // if either dimension is larger than twice the target, shrink appropriately
  if (X > 2 * target_size) {
    Y = (unsigned int) ((double)(Y) * (2.0 * (double)target_size / (double)(X)));
    X = 2 * target_size;
  } else if (Y > 2 * target_size) {
    X = (unsigned int) ((double)(X) * (2.0 * (double)target_size / (double)(Y)));
    Y = 2 * target_size;
  }

  // make sure minimum dimension is sufficient to see
  if (X < 20) {
    X = 20;
  }
  if (Y < 20) {
    Y = 20;
  }
}

} // end of namespace pism

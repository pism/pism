// Copyright (C) 2011 Torsten Albrecht and Constantine Khroulev
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

#include "SSAFD.hh"

PetscErrorCode SSAFD_PIK::assemble_matrix(bool include_basal_shear, Mat A) {
  PetscErrorCode ierr;

  // put the new matrix assembly here

  return 0;
}

PetscErrorCode SSAFD_PIK::assemble_rhs(Vec rhs) {
  PetscErrorCode ierr;

  // put the new RHS assembly here

  return 0;
}


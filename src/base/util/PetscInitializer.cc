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

#include "PetscInitializer.hh"

#include <petsc.h>
#include <mpi.h>
#include <cstdio>

namespace pism {

PetscInitializer::PetscInitializer(int argc, char **argv, const char *help) {

  PetscErrorCode ierr = 0;
  PetscBool initialized = PETSC_FALSE;

  ierr = PetscInitialized(&initialized); CHKERRCONTINUE(ierr);
  if (initialized == PETSC_FALSE) {
    ierr = PetscInitialize(&argc, &argv, NULL, help); CHKERRCONTINUE(ierr);
    if (ierr != 0) {
      printf("PETSc initialization failed. Aborting...\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
}

PetscInitializer::~PetscInitializer() {
  PetscErrorCode ierr = 0;
  PetscBool initialized = PETSC_FALSE;
  ierr = PetscInitialized(&initialized); CHKERRCONTINUE(ierr);
  if (initialized == PETSC_TRUE) {
    // there is nothing we can do if this fails
    ierr = PetscFinalize(); CHKERRCONTINUE(ierr);
  }
}

} // end of namespace pism

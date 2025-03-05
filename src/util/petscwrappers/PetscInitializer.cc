/* Copyright (C) 2014, 2015, 2017, 2023, 2024 PISM Authors
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

#include "pism/util/petscwrappers/PetscInitializer.hh"

#include <petscsys.h>
#include <mpi.h>
#include <cstdio>

#include "pism/util/error_handling.hh"

#if (Pism_USE_YAC_INTERPOLATION == 1)
#include "pism/util/yaxt_wrapper.h"
#endif

namespace pism {
namespace petsc {

Initializer::Initializer(int argc, char **argv, const char *help) {

  PetscErrorCode ierr = 0;
  PetscBool initialized = PETSC_FALSE;

  ierr = PetscInitialized(&initialized);
  PISM_CHK(ierr, "PetscInitialized");

  if (initialized == PETSC_FALSE) {
    ierr = PetscInitialize(&argc, &argv, NULL, help);
    PISM_CHK(ierr, "PetscInitialize");

    if (ierr != 0) {
      printf("PETSc initialization failed. Aborting...\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

#if (Pism_USE_YAC_INTERPOLATION == 1)
  int yaxt_initialized = pism_yaxt_initialized();
  if (yaxt_initialized != 1) {
    pism_yaxt_initialize(PETSC_COMM_WORLD);
  }
#endif
}

Initializer::~Initializer() {

#if (Pism_USE_YAC_INTERPOLATION == 1)
  int yaxt_initialized = pism_yaxt_initialized();
  int yaxt_finalized = pism_yaxt_finalized();
  if (yaxt_initialized == 1 and yaxt_finalized != 1) {
    pism_yaxt_finalize();
  }
#endif

  PetscBool petsc_initialized = PETSC_FALSE;
  PetscErrorCode ierr = PetscInitialized(&petsc_initialized); CHKERRCONTINUE(ierr);

  if (petsc_initialized == PETSC_TRUE) {
    // there is nothing we can do if this fails
    ierr = PetscFinalize(); CHKERRCONTINUE(ierr);
  }
}

} // end of namespace petsc
} // end of namespace pism

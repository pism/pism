/* Copyright (C) 2026 PISM Authors
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

#include <mpi.h>
#include <petscsys.h>

#include "pism/util/pism_initialization.hh"

#include "pism/util/error_handling.hh"

#if (Pism_USE_YAC == 1)
#include "pism/util/yaxt_wrapper.h"

extern "C" {
#include "yac.h"
}
#endif

namespace pism {

//! true if PISM initialized YAC, false otherwise
static bool s_pism_yac_initialized = false;

//! true if PISM should finalize MPI, false otherwise
static bool s_pism_finalize_mpi = false;

static void pism_yac_error_handler(MPI_Comm /* unused */, const char *msg, const char *source,
                                   int line) {
  throw pism::RuntimeError::formatted(pism::ErrorLocation(source, line), "YAC error: %s", msg);
}

void initialize(int argc, char **argv, const char *help) {

  // Initialize MPI only if it was not done yet
  {
    int flag = 0;
    MPI_Initialized(&flag);

    if (flag == 0) {
      MPI_Init(&argc, &argv);
      s_pism_finalize_mpi = true;
    } else {
      s_pism_finalize_mpi = false;
    }
  }

#if (Pism_USE_YAC == 1)
  // YAXT
  {
    int yaxt_initialized = pism_yaxt_initialized();
    if (yaxt_initialized != 1) {
      pism_yaxt_initialize(PETSC_COMM_WORLD);
    }
  }

  // YAC
  if (not s_pism_yac_initialized) {
    int yac_comp_id;
    const char *start_datetime = "1850-01-01T00:00:00";
    const char *end_datetime   = "1850-12-31T00:00:00";

    yac_cinit();

    yac_set_abort_handler((yac_abort_func)pism_yac_error_handler);

    yac_cdef_calendar(YAC_YEAR_OF_365_DAYS);
    yac_cdef_datetime(start_datetime, end_datetime);

    yac_cdef_comp("pism", &yac_comp_id);
    yac_cget_comp_comm(yac_comp_id, &PETSC_COMM_WORLD);

    s_pism_yac_initialized = true;
  }
#endif

  // PETSc
  {
    PetscBool initialized = PETSC_FALSE;
    PetscErrorCode ierr   = PetscInitialized(&initialized);
    PISM_CHK(ierr, "PetscInitialized");

    if (initialized == PETSC_FALSE) {
      ierr = PetscInitialize(&argc, &argv, NULL, help);
      PISM_CHK(ierr, "PetscInitialize");

      if (ierr != 0) {
        printf("PETSc initialization failed. Aborting...\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
  }
}

void initialize(const char *help) {
  std::vector<char *> argv = { nullptr };
  initialize(0, argv.data(), help);
}

void initialize_options(const std::vector<std::string> &args) {
  int argc = (int) args.size();
  std::vector<const char*> argv(argc + 1);

  for (int i = 0; i < argc; ++i) {
    argv[i] = args[i].c_str();
  }
  argv[argc] = nullptr;

  // Note: this const_cast is needed for compatibility with earlier PETSc versions.
  PetscOptionsInsertArgs(NULL, argc, const_cast<char**>(argv.data()));
}


void finalize() {

  // PETSc
  {
    PetscBool petsc_initialized = PETSC_FALSE;
    PetscErrorCode ierr         = PetscInitialized(&petsc_initialized);
    CHKERRCONTINUE(ierr);

    if (petsc_initialized == PETSC_TRUE) {
      // there is nothing we can do if this fails
      ierr = PetscFinalize();
      CHKERRCONTINUE(ierr);
    }
  }

#if (Pism_USE_YAC == 1)
  // YAC
  if (s_pism_yac_initialized) {
    yac_cfinalize();
  }

  // YAXT
  {
    int yaxt_initialized = pism_yaxt_initialized();
    int yaxt_finalized   = pism_yaxt_finalized();
    if (yaxt_initialized == 1 and yaxt_finalized != 1) {
      pism_yaxt_finalize();
    }
  }
#endif

  // MPI
  if (s_pism_finalize_mpi) {
    MPI_Finalize();
  }
}

}

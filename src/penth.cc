// Copyright (C) 2004-2010 Ed Bueler, Andy Aschwanden, and Constantine Khroulev
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

static char help[] = "PENTH IS GONE.  DO NOT USE.\n";

#include <petsc.h>

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;

  ierr = PetscPrintf(com,
      "\n\nPENTH DOES NOT EXIST.  DO NOT TRY TO USE IT.\n\n\n"); CHKERRQ(ierr);
  ierr = PetscEnd(); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}


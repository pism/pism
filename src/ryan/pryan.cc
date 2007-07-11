// Copyright (C) 2007 Ryan Woodard

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

static char help[] =
  "Ice sheet driver for EISMINT II experiment H with random accumulation.\n";

#include <petscda.h>
#include "grid.hh"
#include "materials.hh"
#include "iceRYANModel.hh"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(1,com, "PRYAN (perturbed accumulation EISMINT II mode)\n"); CHKERRQ(ierr);

    PetscInt   flowlawNumber = 0; // use Paterson-Budd by default
    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);
    
    IceGrid    g(com, rank, size);
    IceType*   ice;
    IceRYANModel   m(g, *ice);

    ierr = m.setFromOptions(); CHKERRQ(ierr);
    ierr = m.initFromOptions(); CHKERRQ(ierr);

    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "done with run ... "); CHKERRQ(ierr);
    ierr = m.writeFiles("ryan_exper"); CHKERRQ(ierr);    
    ierr = verbPrintf(2,com, " ... done.\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

// Copyright (C) 2007, 2008 Ed Bueler and Constantine Khroulev
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

static char help[] =
  "Ice sheet driver for ISMIP-HEINO simplified geometry\n"
  "intercomparison simulations.  Not recommended as an intercomparison.\n";

#include <cstring>
#include <petscbag.h>

// following includes assume this driver is in pism/src/; note Makefile must 
//   also be modified to build "pismheino" executable:
#include "base/grid.hh"
#include "base/materials.hh"
#include "base/iceModel.hh"
#include "ismip/heino/iceHEINOModel.hh"

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
    IceGrid    g(com, rank, size);
    IceType*   ice = PETSC_NULL;
    
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "PISMHEINO (ISMIP-HEINO mode)\n"); CHKERRQ(ierr);

    ierr = userChoosesIceType(com, ice); CHKERRQ(ierr);  // allocates ice
    
    // call constructors on all three, but m will point to the one we use
    IceHEINOModel  mHEINO(g, ice);

    mHEINO.setflowlawNumber(flowlawNumber);
    ierr = mHEINO.setFromOptions(); CHKERRQ(ierr);
    ierr = mHEINO.initFromOptions(); CHKERRQ(ierr);

    ierr = mHEINO.setExecName("pismheino"); CHKERRQ(ierr);
    ierr = mHEINO.run(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "done with run ... \n"); CHKERRQ(ierr);
    ierr = mHEINO.writeFiles("heino_exper.nc"); CHKERRQ(ierr);
    
    ierr = mHEINO.simpFinalize(); CHKERRQ(ierr);

    delete ice;
    ierr = verbPrintf(2,com, "\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

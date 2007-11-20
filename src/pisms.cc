// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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
  "Ice sheet driver for EISMINT II, ISMIP-HEINO, and other simplified geometry\n"
  "ice sheet simulations.\n";

#include <cstring>
#include <petscbag.h>
#include "base/grid.hh"
#include "base/materials.hh"
#include "base/iceModel.hh"
#include "eismint/iceEISModel.hh"
#include "eismint/iceEISplModel.hh"
#include "ismip/iceHEINOModel.hh"

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
    IceType*   ice;
    PetscInt   flowlawNumber = 0; // use Paterson-Budd by default
    
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(1,com, "PISMS (simplified geometry mode)\n"); CHKERRQ(ierr);

    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);
    
    // call constructors on all three, but m will point to the one we use
    IceEISModel    mEISII(g, *ice);
    IceHEINOModel  mHEINO(g, *ice);
    IceEISplModel  mEISpl(g, *ice);
    IceModel*      m;

    PetscTruth  EISIIchosen, ISMIPchosen, EISplchosen;
    /* This option determines the single character name of EISMINT II experiments:
    "-eisII F", for example. */
    ierr = PetscOptionsHasName(PETSC_NULL, "-eisII", &EISIIchosen);
              CHKERRQ(ierr);
    /* This option chooses plastic till modification of EISMINT II experiment I. */
    ierr = PetscOptionsHasName(PETSC_NULL, "-eis2Ipl", &EISplchosen);
              CHKERRQ(ierr);
    /* This option chooses ISMIP; "-ismip H" is ISMIP-HEINO and none others are implemented */
    ierr = PetscOptionsHasName(PETSC_NULL, "-ismip", &ISMIPchosen);
              CHKERRQ(ierr);
    
    if ( (EISIIchosen == PETSC_TRUE) && (ISMIPchosen == PETSC_FALSE)
         && (EISplchosen == PETSC_FALSE) ) {
      mEISII.setFlowLawNumber(flowlawNumber);
      ierr = mEISII.setFromOptions(); CHKERRQ(ierr);
      ierr = mEISII.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mEISII;
    } else if ( (ISMIPchosen == PETSC_TRUE) && (EISIIchosen == PETSC_FALSE)
                && (EISplchosen == PETSC_FALSE) ) {
      mHEINO.setflowlawNumber(flowlawNumber);
      ierr = mHEINO.setFromOptions(); CHKERRQ(ierr);
      ierr = mHEINO.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mHEINO;
    } else if ( (EISplchosen == PETSC_TRUE) && (EISIIchosen == PETSC_FALSE)
                && (ISMIPchosen == PETSC_FALSE) ) {
      mEISpl.setFlowLawNumber(flowlawNumber);
      ierr = mEISpl.setFromOptions(); CHKERRQ(ierr);
      ierr = mEISpl.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mEISpl;
    } else {
      SETERRQ(1,
        "PISMS called with invalid, contradictory, or no simplified geometry experiment chosen");
    }

    ierr = m->run(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "done with run ... \n"); CHKERRQ(ierr);
    ierr = m->writeFiles("simp_exper"); CHKERRQ(ierr);
    
    if (ISMIPchosen == PETSC_TRUE) {
      ierr = mHEINO.simpFinalize(); CHKERRQ(ierr);
    }
    ierr = verbPrintf(2,com, " ... done.\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

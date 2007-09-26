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
    IceEISModel   mEISII(g, *ice);
    IceHEINOModel mHEINO(g, *ice);
    IceModel*     m;

    char        expername[20]; //ignored here; see IceEISModel,IceHEINOModel::initFromOptions
    PetscTruth  EISIIchosen, ISMIPchosen, ROSSchosen;
    /* This option determines the single character name of EISMINT II experiments:
    "-eisII F", for example. */
    ierr = PetscOptionsGetString(PETSC_NULL, "-eisII", expername, 1, &EISIIchosen);
              CHKERRQ(ierr);
    /* This option chooses ISMIP; "-ismip H" is ISMIP-HEINO and none others are implemented */
    ierr = PetscOptionsGetString(PETSC_NULL, "-ismip", expername, 1, &ISMIPchosen);
              CHKERRQ(ierr);
    /* This option chooses EISMINT ROSS; "-ross" */
    ierr = PetscOptionsHasName(PETSC_NULL, "-ross", &ROSSchosen); CHKERRQ(ierr);
    
    if ((EISIIchosen == PETSC_TRUE) && (ISMIPchosen == PETSC_FALSE) && (ROSSchosen == PETSC_FALSE)) {
      mEISII.setFlowLawNumber(flowlawNumber);
      ierr = mEISII.setFromOptions(); CHKERRQ(ierr);
      ierr = mEISII.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mEISII;
    } else if ((ISMIPchosen == PETSC_TRUE) && (EISIIchosen == PETSC_FALSE) && (ROSSchosen == PETSC_FALSE)) {
      mHEINO.setflowlawNumber(flowlawNumber);
      ierr = mHEINO.setFromOptions(); CHKERRQ(ierr);
      ierr = mHEINO.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mHEINO;
    } else if ((ROSSchosen == PETSC_TRUE) && (EISIIchosen == PETSC_FALSE) && (ISMIPchosen == PETSC_FALSE)) {
        ierr = PetscPrintf(com, "pisms no longer runs EISMINT ROSS; use pismd\nending ... \n"); CHKERRQ(ierr);
        PetscEnd();      
    } else {
      SETERRQ(1,"PISMS called with invalid, contradictory, or no experiment chosen");
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

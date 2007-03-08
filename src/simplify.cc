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
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"
#include "iceEISModel.hh"
#include "iceHEINOModel.hh"
#include "iceROSSModel.hh"

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
    
    ierr = PetscPrintf(com, "PISMS (simplified geometry mode)\n"); CHKERRQ(ierr);

    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);
    
    // call constructors on both, but have a pointer for either
    IceEISModel   mEISII(g, *ice);
    IceHEINOModel mHEINO(g, *ice);
    IceROSSModel  mROSS(g, *ice);
    IceModel*     m;

    char        expername[20]; //ignored here; see IceEISModel,IceHEINOModel::initFromOptions
    PetscTruth  EISIIchosen, ISMIPchosen, ROSSchosen;
    /* This option determines the single character name of EISMINT II experiments:
    "-eisII F", for example. */
    ierr = PetscOptionsGetString(PETSC_NULL, "-eisII", expername, 1, &EISIIchosen);
              CHKERRQ(ierr);
    /* This option chooses ISMIP; "-ismip H" is ISMIP-HEINO */
    ierr = PetscOptionsGetString(PETSC_NULL, "-ismip", expername, 1, &ISMIPchosen);
              CHKERRQ(ierr);
    /* This option chooses EISMINT ROSS; "-ross" */
    ierr = PetscOptionsHasName(PETSC_NULL, "-ross", &ROSSchosen); CHKERRQ(ierr);
    
    if ((EISIIchosen == PETSC_TRUE) && (ISMIPchosen == PETSC_FALSE) && (ROSSchosen == PETSC_FALSE)) {
      mEISII.setflowlawNumber(flowlawNumber);
      ierr = mEISII.setFromOptions(); CHKERRQ(ierr);
      ierr = mEISII.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mEISII;
    } else if ((ISMIPchosen == PETSC_TRUE) && (EISIIchosen == PETSC_FALSE) && (ROSSchosen == PETSC_FALSE)) {
      mHEINO.setflowlawNumber(flowlawNumber);
      ierr = mHEINO.setFromOptions(); CHKERRQ(ierr);
      ierr = mHEINO.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mHEINO;
    } else if ((ROSSchosen == PETSC_TRUE) && (EISIIchosen == PETSC_FALSE) && (ISMIPchosen == PETSC_FALSE)) {
      if (size > 1) {
        ierr = PetscPrintf(com, "EISMINT ROSS only works on a single processor!\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(com,"ending ... \n"); CHKERRQ(ierr);
        PetscEnd();
      }
      mROSS.setflowlawNumber(flowlawNumber);
      ierr = mROSS.setFromOptions(); CHKERRQ(ierr);
      ierr = mROSS.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mROSS;
    } else {
      SETERRQ(1,"PISMS called with invalid, contradictory, or no experiment chosen");
    }

    ierr = m->setSoundingFromOptions(); CHKERRQ(ierr);

    ierr = m->run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "done with run ... "); CHKERRQ(ierr);

    if (ROSSchosen == PETSC_FALSE) {
      ierr = m->writeFiles("simp_exper"); CHKERRQ(ierr);
    }
    
    if (ISMIPchosen == PETSC_TRUE) {
      ierr = mHEINO.simpFinalize(); CHKERRQ(ierr);
    }
    ierr = verbPrintf(2,com, " ... done.\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

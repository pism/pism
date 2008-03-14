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
  "Ice sheet driver for EISMINT II and MISMIP simplified geometry\n"
  "intercomparison simulations.\n";

#include <cstring>
#include <petscbag.h>
#include "base/grid.hh"
#include "base/materials.hh"
#include "base/iceModel.hh"
#include "eismint/iceEISModel.hh"
#include "eismint/iceEISplModel.hh"
#include "ismip/iceMISMIPModel.hh"

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
    MISMIPIce* mismipice = new MISMIPIce;
    
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "PISMS (simplified geometry mode)\n"); CHKERRQ(ierr);

    ierr = userChoosesIceType(com, ice); CHKERRQ(ierr);  // allocates ice
    
    // call constructors on all three, but m will point to the one we use
    IceEISModel    mEISII(g, ice);
    IceEISplModel  mEISpl(g, ice);
    IceMISMIPModel mMISMIP(g, mismipice, mismipice);
    IceModel*      m;

    PetscTruth  EISIIchosen, EISplchosen, MISMIPchosen;
    /* This option determines the single character name of EISMINT II experiments:
    "-eisII F", for example. */
    ierr = PetscOptionsHasName(PETSC_NULL, "-eisII", &EISIIchosen);
              CHKERRQ(ierr);
    /* This option chooses plastic till modification of EISMINT II experiment A or I. */
    ierr = PetscOptionsHasName(PETSC_NULL, "-eis2pl", &EISplchosen);
              CHKERRQ(ierr);
    /* This option chooses MISMIP; "-mismip N" is experiment N in MISMIP; N=1,2,3 */
    ierr = PetscOptionsHasName(PETSC_NULL, "-mismip", &MISMIPchosen);
              CHKERRQ(ierr);
    
    int  choiceSum = (int) EISIIchosen + (int) EISplchosen + (int) MISMIPchosen;
    if (choiceSum == 0) {
      SETERRQ(1,"PISMS called with no simplified geometry experiment chosen");
    } else if (choiceSum > 1) {
      SETERRQ(2,"PISMS called with more than one simplified geometry experiment chosen");
    }
    
    if (EISIIchosen == PETSC_TRUE) {
      ierr = mEISII.setFromOptions(); CHKERRQ(ierr);
      ierr = mEISII.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mEISII;
    } else if (EISplchosen == PETSC_TRUE) {
      ierr = mEISpl.setFromOptions(); CHKERRQ(ierr);
      ierr = mEISpl.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mEISpl;
    } else if (MISMIPchosen == PETSC_TRUE) {
      ierr = mMISMIP.setFromOptions(); CHKERRQ(ierr);
      ierr = mMISMIP.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mMISMIP;
    } else {
      SETERRQ(3,"PISMS: how did I get here?");
    }

//    ierr = m->testIceModelVec3(); CHKERRQ(ierr);
//    ierr = m->testIceModelVec3Bedrock(); CHKERRQ(ierr);

    ierr = m->setExecName("pisms"); CHKERRQ(ierr);
    ierr = m->run(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "done with run ... \n"); CHKERRQ(ierr);
    ierr = m->writeFiles("simp_exper"); CHKERRQ(ierr);
    
    delete ice;
    delete mismipice;
    ierr = verbPrintf(2,com, "\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

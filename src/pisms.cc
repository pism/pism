// Copyright (C) 2004-2008 Jed Brown, Ed Bueler and Constantine Khroulev
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
  "Ice sheet driver for EISMINT II, MISMIP, and other constant climate, simplified geometry\n"
  "intercomparison simulations.\n";

#include <cstring>
#include <petscbag.h>
#include "base/grid.hh"
#include "base/materials.hh"
#include "base/iceModel.hh"
#include "coupler/pccoupler.hh"
#include "eismint/iceEISModel.hh"
#include "eismint/icePSTexModel.hh"
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
    PISMConstAtmosCoupler pcac;
    PISMConstOceanCoupler pcoc;
    
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "PISMS  (simplified geometry mode)\n"); CHKERRQ(ierr);

    ierr = userChoosesIceType(com, ice); CHKERRQ(ierr);  // allocates ice
    
    // call constructors on all three, but m will point to the one we use
    IceEISModel    mEISII(g, ice);
    IcePSTexModel  mPSTex(g, ice);
    IceMISMIPModel mMISMIP(g, mismipice);
    IceModel*      m;

    PetscTruth  pddSet, EISIIchosen, PSTexchosen, MISMIPchosen;

    ierr = PetscOptionsHasName(PETSC_NULL, "-pdd", &pddSet); CHKERRQ(ierr);
    if (pddSet == PETSC_TRUE) {
      ierr = PetscPrintf(com, "PISM ERROR: -pdd is not currently allowed as option to pisms\n");
                CHKERRQ(ierr);
      PetscEnd();
    }

    /* This option determines the single character name of EISMINT II experiments:
    "-eisII F", for example. */
    ierr = PetscOptionsHasName(PETSC_NULL, "-eisII", &EISIIchosen);
              CHKERRQ(ierr);
    /* This option chooses Plastic till ice Stream with Thermocoupling experiment. */
    ierr = PetscOptionsHasName(PETSC_NULL, "-pst", &PSTexchosen);
              CHKERRQ(ierr);
    /* This option chooses MISMIP; "-mismip N" is experiment N in MISMIP; N=1,2,3 */
    ierr = PetscOptionsHasName(PETSC_NULL, "-mismip", &MISMIPchosen);
              CHKERRQ(ierr);
    
    int  choiceSum = (int) EISIIchosen + (int) PSTexchosen + (int) MISMIPchosen;
    if (choiceSum == 0) {
      ierr = PetscPrintf(com, 
         "PISM ERROR: pisms called with no simplified geometry experiment chosen.\n");
         CHKERRQ(ierr);
      PetscEnd();
    } else if (choiceSum > 1) {
      ierr = PetscPrintf(com,
         "PISM ERROR: pisms called with more than one simplified geometry experiment chosen.\n");
         CHKERRQ(ierr);
      PetscEnd();
    }
    
    if (EISIIchosen == PETSC_TRUE) {
      m = (IceModel*) &mEISII;
    } else if (PSTexchosen == PETSC_TRUE) {
      m = (IceModel*) &mPSTex;
    } else if (MISMIPchosen == PETSC_TRUE) {
      m = (IceModel*) &mMISMIP;
    } else {
      SETERRQ(3,"PISMS: how did I get here?");
    }

    // Put these in a new 'software testing' executable:
//     ierr = m->testIceModelVec3(); CHKERRQ(ierr);
//     ierr = m->testIceModelVec3Bedrock(); CHKERRQ(ierr);

    ierr = m->setFromOptions(); CHKERRQ(ierr);
    pcac.initializeFromFile = false;  // climate will always come from intercomparison formulas, for pisms
    ierr = m->attachAtmospherePCC(pcac); CHKERRQ(ierr);
    ierr = m->attachOceanPCC(pcoc); CHKERRQ(ierr);
    ierr = m->initFromOptions(); CHKERRQ(ierr);
    ierr = m->setExecName("pisms"); CHKERRQ(ierr);
    ierr = m->run(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "done with run ... \n"); CHKERRQ(ierr);
    ierr = m->writeFiles("simp_exper.nc"); CHKERRQ(ierr);
    if (MISMIPchosen == PETSC_TRUE) {
      ierr = mMISMIP.writeMISMIPFinalFiles(); CHKERRQ(ierr);
    }

    delete mismipice;
    delete ice;
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

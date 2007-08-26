// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

static char help[] =
"Ice sheet driver for SIA and SSA verification.  Uses exact solutions to various coupled\n"
"subsystems.  Currently implements tests A, B, C, D, E, F, G, H, I, L.\n\n";

#include <cstring>
#include <cstdio>
#include <petscda.h>
#include <petscbag.h>
#include "base/grid.hh"
#include "base/materials.hh"
#include "verif/iceCompModel.hh"
#include "verif/iceExactSSAModel.hh"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  MPI_Comm        com;
  PetscMPIInt     rank, size;

  PetscInitialize(&argc, &argv, PETSC_NULL, help);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);
      
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    IceGrid      g(com, rank, size);
    IceType*     ice;
    PetscInt     flowlawNumber = 1;  // use cold part of Paterson-Budd by default
    char         testname[20];
    PetscTruth   testchosen, dontReport;

    ierr = verbPrintf(1, com, "PISMV (verification mode)\n"); CHKERRQ(ierr);
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);
    
    // determine test (and whether to report error)
    ierr = PetscOptionsGetString(PETSC_NULL, "-test", testname, 1, &testchosen); CHKERRQ(ierr);
    char test = testname[0];  // only use the first letter
    if (testchosen == PETSC_FALSE)     test = 'A';           // default to test A
    if ((test >= 'a') && (test <= 'z'))    test += 'A'-'a';  // capitalize if lower    
    ierr = PetscOptionsHasName(PETSC_NULL, "-no_report", &dontReport); CHKERRQ(ierr);

    // actually construct and run one of the derived classes of IceModel
    if (test == 'I') {
      // run derived class (of IceModel) for plastic till ice stream
      IceExactSSAModel mSSA(g, *ice, test);  
      ierr = mSSA.setFromOptions(); CHKERRQ(ierr);
      ierr = mSSA.initFromOptions(); CHKERRQ(ierr);
      ierr = mSSA.diagnosticRun(); CHKERRQ(ierr);
      if (dontReport == PETSC_FALSE) {
        ierr = mSSA.reportErrors();  CHKERRQ(ierr);
      }
      ierr = mSSA.writeFiles("verify",PETSC_TRUE); CHKERRQ(ierr);
    } else { // all other tests are compensatory; 
             // if test is invalid then mComp.initFromOptions() will bonk
      // run derived class (of IceModel) for compensatory source term cases 
      // (i.e. compensatory accumulation or compensatory heating)
      ThermoGlenArrIce*   tgaice = (ThermoGlenArrIce*) ice;
      IceCompModel        mComp(g, *tgaice, test);
      ierr = mComp.setFromOptions(); CHKERRQ(ierr);
      ierr = mComp.initFromOptions(); CHKERRQ(ierr);
      ierr = mComp.run(); CHKERRQ(ierr);
      if (dontReport == PETSC_FALSE) {
        if ((flowlawNumber != 1) && ((test == 'F') || (test == 'G'))) {
            ierr = verbPrintf(1,com, 
                "pismv WARNING: flow law must be cold part of Paterson-Budd ('-law 1')\n"
                "   for reported errors in test %c to be meaningful!\n", test);
                CHKERRQ(ierr);
        }
        ierr = mComp.reportErrors();  CHKERRQ(ierr);
      }
      ierr = mComp.writeFiles("verify",PETSC_FALSE); CHKERRQ(ierr);
    }

    ierr = verbPrintf(2,com, " ... done\n"); CHKERRQ(ierr);
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

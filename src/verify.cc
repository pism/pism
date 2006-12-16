// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
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
"Ice sheet driver for SIA verification.  Uses exact solutions to various coupled\n"
"subsystems.  Currently implements tests A, B, C, D, E, F, G, H, I.\n\n";

#include <cstring>
#include <cstdio>
#include <petscda.h>
#include <petscbag.h>

#include "grid.hh"
#include "iceCompModel.hh"
#include "iceExactStreamModel.hh"

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
    IceGrid             g(com, rank, size);

    IceType*            tempice;
    PetscInt            flowlawNumber = 1;  // use cold part of Paterson-Budd by default

    ierr = PetscPrintf(com, "PISMV (verification mode)\n"); CHKERRQ(ierr);
    
    ierr = getFlowLawFromUser(com, tempice, flowlawNumber); CHKERRQ(ierr);

    // call constructors on both, but have a pointer for either
    ThermoGlenArrIce*   ice = (ThermoGlenArrIce*) tempice;
    IceCompModel        mComp(g, *ice);        // derived class of IceModel
    IceExactStreamModel mStream(g, *tempice);  // ditto
    IceModel*           m;

    //  determine which derived class to initialize based on test name
    char        testname[20];
    PetscTruth  testchosen;
    ierr = PetscOptionsGetString(PETSC_NULL, "-test", testname, 1, &testchosen); CHKERRQ(ierr);
    char temp = testname[0];
    if (testchosen == PETSC_FALSE)     temp = 'A';           // default to test A
    if ((temp >= 'a') && (temp <= 'z'))    temp += 'A'-'a';  // capitalize if lower    
    if ((temp >= 'A') && (temp <= 'H')) {
      ierr = mComp.setFromOptions(); CHKERRQ(ierr);
      ierr = mComp.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mComp;
    } else if (temp == 'I') {
      ierr = mStream.setFromOptions(); CHKERRQ(ierr);
      ierr = mStream.initFromOptions(); CHKERRQ(ierr);
      m = (IceModel*) &mStream;
    } else {
      SETERRQ(1,"(PISMV; verify.cc) ERROR: desired test NOT IMPLEMENTED\n");
    }

    ierr = m->setSoundingFromOptions();  CHKERRQ(ierr);

    ierr = m->run(); CHKERRQ(ierr);

    ierr = m->verbPrintf(2,com, "done with run\n"); CHKERRQ(ierr);
    
    /* Whether to report error at end. */
    PetscTruth dontReport;
    ierr = PetscOptionsHasName(PETSC_NULL, "-noreport", &dontReport); CHKERRQ(ierr);
    if (dontReport == PETSC_FALSE) {
      if ((temp >= 'A') && (temp <= 'H')) {
        if ((flowlawNumber != 1) && ((temp == 'F') || (temp == 'G'))) {
          ierr = PetscPrintf(com, 
                "verify WARNING: flow law must be cold part of Paterson-Budd ('-law 1')\n"
                "   for reported errors in tests F and G to be meaningful!\n"); CHKERRQ(ierr);
        }
        ierr = mComp.reportErrors();  CHKERRQ(ierr);
      } else if (temp == 'I') {
        ierr = mStream.reportErrors();  CHKERRQ(ierr);
      } else {
        SETERRQ(2,"(PISMV; verify.cc) ERROR: can not report error for desired test\n");
      }
    }

    ierr = m->writeFiles("verify"); CHKERRQ(ierr);
    ierr = m->verbPrintf(2,com, " ... done\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

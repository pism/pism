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
"Ice sheet driver for SIA verification.  Uses exact solutions from Bueler, Brown,\n"
"and Lingle (2006) ['Accuracy analysis of a numerical scheme for thermocoupled ice\n"
"sheets', in preparation] and Bueler et al (2005) ['Exact solutions and \n"
"verification of numerical models for isothermal ice sheets', J. Glaciology\n"
"51(173) 291--306].  Currently implements tests B, C, D, F, and G.\n"
"Computes compensatory accumulation M, heating Sigma_c terms at\n"
"each time step, and exact solution at each step if desired.\n\n";

#include <cstring>
#include <cstdio>
#include <petscda.h>
#include <petscbag.h>

#include "grid.hh"
#include "iceCompModel.hh"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  MPI_Comm    com;
  PetscMPIInt rank, size;
  PetscTruth  dontReport;

  PetscInitialize(&argc, &argv, PETSC_NULL, help);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);
      
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    IceGrid             g(com, rank, size);

    ThermoGlenArrIce*   ice;
    IceType*            tempice;
    PetscInt            flowlawNumber = 1;  // use cold part of Paterson-Budd by default

    ierr = PetscPrintf(com, "PISMV (verification mode)\n"); CHKERRQ(ierr);
    
    ierr = getFlowLawFromUser(com, tempice, flowlawNumber); CHKERRQ(ierr);
    if (flowlawNumber != 1) {
      ierr = PetscPrintf(com, 
              "verify WARNING: flow law must be cold part of Paterson-Budd ('-law 1')\n"
              "   for reported errors to be meaningful!\n"); CHKERRQ(ierr);
    }    
    ice = (ThermoGlenArrIce*) tempice;
    IceCompModel m(g, *ice);
    ierr = m.setFromOptions();  CHKERRQ(ierr);
    ierr = m.initFromOptions();  CHKERRQ(ierr);
    ierr = m.setSoundingFromOptions();  CHKERRQ(ierr);

    ierr = m.run(); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "done with run\n"); CHKERRQ(ierr);
    
    /* Whether to report error at end. */
    ierr = PetscOptionsHasName(PETSC_NULL, "-noreport", &dontReport); CHKERRQ(ierr);
    if (dontReport == PETSC_FALSE) {
      if (flowlawNumber != 1) {
        ierr = PetscPrintf(com, 
                "verify WARNING: flow law must be cold part of Paterson-Budd ('-law 1')\n"
                "   for reported errors to be meaningful!\n"); CHKERRQ(ierr);
      }
      ierr = m.reportErrors();  CHKERRQ(ierr);
    }

    m.writeFiles("verify"); CHKERRQ(ierr);

    ierr = PetscPrintf(com, " ... done\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

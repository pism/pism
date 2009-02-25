// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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
"Ice sheet driver for PISM (SIA and SSA) verification.  Uses exact solutions to various\n"
"  coupled subsystems.  Computes difference between exact solution and numerical solution.\n"
"  Can also just compute exact solution (-eo).\n"
"  Currently implements tests A, B, C, D, E, F, G, H, I, J, L, M.\n\n";

#include <cstring>
#include <cstdio>
#include <petscda.h>
#include <petscbag.h>
#include "base/grid.hh"
#include "base/materials.hh"
#include "coupler/pccoupler.hh"
#include "verif/iceCompModel.hh"
#include "verif/iceExactSSAModel.hh"
#include "verif/iceCalvBCModel.hh"

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
    PISMAtmosphereCoupler pac;
    PISMConstOceanCoupler pcoc;

    char         testname[20];
    PetscTruth   testchosen, dontReport;

    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(1, com, "PISMV  (verification mode)\n"); CHKERRQ(ierr);

    // determine test (and whether to report error)
    ierr = PetscOptionsGetString(PETSC_NULL, "-test", testname, 1, &testchosen); CHKERRQ(ierr);
    char test = testname[0];  // only use the first letter
    if (testchosen == PETSC_FALSE)         test = 'A';       // default to test A
    if ((test >= 'a') && (test <= 'z'))    test += 'A'-'a';  // capitalize if lower    
    ierr = PetscOptionsHasName(PETSC_NULL, "-no_report", &dontReport); CHKERRQ(ierr);

    // actually construct and run one of the derived classes of IceModel
    if (test == '0') {
      // run derived class for test M which includes new calving front stress
      //   boundary condition implementation
      ierr = verbPrintf(1,com, "!!!!!!!! USING IceCalvBCModel TO DO test M !!!!!!!!\n"); CHKERRQ(ierr);
      IceCalvBCModel mCBC(g, 'M');
      ierr = mCBC.getIceFactory().setType(ICE_ARR);CHKERRQ(ierr);
      ierr = mCBC.setExecName("pismv"); CHKERRQ(ierr);
      ierr = mCBC.attachAtmospherePCC(pac); CHKERRQ(ierr);
      ierr = mCBC.attachOceanPCC(pcoc); CHKERRQ(ierr);
      ierr = mCBC.setFromOptions(); CHKERRQ(ierr);
      ierr = mCBC.initFromOptions(); CHKERRQ(ierr);
      ierr = mCBC.diagnosticRun(); CHKERRQ(ierr);
      ierr = verbPrintf(2,com, "done with diagnostic run\n"); CHKERRQ(ierr);
      if (dontReport == PETSC_FALSE) {
        ierr = mCBC.reportErrors();  CHKERRQ(ierr);
      }
      ierr = mCBC.writeFiles("verify.nc",PETSC_TRUE); CHKERRQ(ierr);
      ierr = mCBC.writeCFfields("verify.nc"); CHKERRQ(ierr); // add three more fields
    } else if ((test == 'I') || (test == 'J') || (test == 'M')) {
      // run derived class for plastic till ice stream, or linearized ice shelf,
      //   or annular ice shelf with calving front
      IceExactSSAModel mSSA(g, test);
      if (test != 'I') {        // Correct errors in test I require CustomGlenIce
        ierr = mSSA.getIceFactory().setType(ICE_ARR);CHKERRQ(ierr);
      }
      ierr = mSSA.setExecName("pismv"); CHKERRQ(ierr);
      ierr = mSSA.attachAtmospherePCC(pac); CHKERRQ(ierr);
      ierr = mSSA.attachOceanPCC(pcoc); CHKERRQ(ierr);
      ierr = mSSA.setFromOptions(); CHKERRQ(ierr);
      ierr = mSSA.initFromOptions(); CHKERRQ(ierr);
      ierr = mSSA.diagnosticRun(); CHKERRQ(ierr);
      ierr = verbPrintf(2,com, "done with diagnostic run\n"); CHKERRQ(ierr);
      if (dontReport == PETSC_FALSE) {
        ierr = mSSA.reportErrors();  CHKERRQ(ierr);
      }
      ierr = mSSA.writeFiles("verify.nc",PETSC_TRUE); CHKERRQ(ierr);
    } else { // run derived class for compensatory source SIA solutions
             // (i.e. compensatory accumulation or compensatory heating)
      IceCompModel       mComp(g, test);
      ierr = mComp.getIceFactory().setType(ICE_ARR);CHKERRQ(ierr);
      ierr = mComp.setExecName("pismv"); CHKERRQ(ierr);
      ierr = mComp.attachAtmospherePCC(pac); CHKERRQ(ierr);
      ierr = mComp.attachOceanPCC(pcoc); CHKERRQ(ierr);
      ierr = mComp.setFromOptions(); CHKERRQ(ierr);
      ThermoGlenArrIce*   tgaice = dynamic_cast<ThermoGlenArrIce*>(mComp.getIce());
      if (!tgaice) SETERRQ(1,"Ice is actually not ThermoGlenArrIce");
      ierr = mComp.initFromOptions(); CHKERRQ(ierr);
      ierr = mComp.run(); CHKERRQ(ierr);
      ierr = verbPrintf(2,com, "done with run\n"); CHKERRQ(ierr);
      if (dontReport == PETSC_FALSE) {
        if (!IceTypeIsPatersonBuddCold(tgaice) && ((test == 'F') || (test == 'G'))) {
            ierr = verbPrintf(1,com, 
                "pismv WARNING: flow law must be cold part of Paterson-Budd ('-ice_type arr')\n"
                "   for reported errors in test %c to be meaningful!\n", test); CHKERRQ(ierr);
        }
        ierr = mComp.reportErrors();  CHKERRQ(ierr);
      }
      ierr = mComp.writeFiles("verify.nc",PETSC_FALSE); CHKERRQ(ierr);
    }
    
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

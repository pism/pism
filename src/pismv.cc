// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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
"  Currently implements tests A, B, C, D, E, F, G, H, I, J, K, L, M.\n\n";

#include <ctype.h>		// toupper
#include <cstring>
#include <cstdio>
#include <petscda.h>
#include <petscbag.h>
#include "base/grid.hh"
#include "base/materials.hh"
#include "verif/iceCompModel.hh"
#include "verif/iceExactSSAModel.hh"
#include "verif/iceCalvBCModel.hh"

#include "coupler/PISMSurface.hh"
#include "coupler/PISMOcean.hh"

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
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2, com, "PISMV %s (verification mode)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    vector<string> required;
    required.push_back("-test");
    ierr = show_usage_check_req_opts(com, "pismv", required,
      "  pismv -test x [-no_report] [-eo] [OTHER PISM & PETSc OPTIONS]\n\n"
      "where:\n"
      "  -test x     verification test (x = A|B|...|L)\n"
      "  -no_report  do not give error report at end of run\n"
      "  -eo         do not do numerical run; exact solution only\n"
      ); CHKERRQ(ierr);

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid      g(com, rank, size);

    // Initialize boundary models:
    PISMSurfaceModel *surface = new PSDummy(g, config);
    PISMOceanModel     *ocean = new POConstant(g, config);

    // determine test (and whether to report error)
    char         testname[20];
    PetscTruth   testchosen;
    ierr = PetscOptionsGetString(PETSC_NULL, "-test", testname, 1, 
                                 &testchosen); CHKERRQ(ierr);
    unsigned char test = testname[0];  // only use the first letter
    if (testchosen == PETSC_FALSE)         test = 'A';       // default to test A
    test = toupper(test);				     // capitalize

    PetscTruth   dontReport;
    ierr = check_option("-no_report", dontReport); CHKERRQ(ierr);

    // actually construct and run one of the derived classes of IceModel
    if (test == '0') {
      // run derived class for test M which includes new calving front stress
      //   boundary condition implementation
      ierr = verbPrintf(1,com, "!!!!!!!! USING IceCalvBCModel TO DO test M !!!!!!!!\n"); CHKERRQ(ierr);
      IceCalvBCModel mCBC(g, config, overrides, 'M');
      ierr = mCBC.setExecName("pismv"); CHKERRQ(ierr);
      mCBC.attach_surface_model(surface);
      mCBC.attach_ocean_model(ocean);

      ierr = mCBC.init(); CHKERRQ(ierr);
      
      ierr = mCBC.diagnosticRun(); CHKERRQ(ierr);
      ierr = verbPrintf(2,com, "done with diagnostic run\n"); CHKERRQ(ierr);
      if (dontReport == PETSC_FALSE) {
        ierr = mCBC.reportErrors();  CHKERRQ(ierr);
      }
      ierr = mCBC.writeFiles("verify.nc"); CHKERRQ(ierr);
      ierr = mCBC.writeCFfields("verify.nc"); CHKERRQ(ierr); // add three more fields
    } else if ((test == 'I') || (test == 'J') || (test == 'M')) {
      // run derived class for plastic till ice stream, or linearized ice shelf,
      //   or annular ice shelf with calving front
      IceExactSSAModel mSSA(g, config, overrides, test);

      ierr = mSSA.setExecName("pismv"); CHKERRQ(ierr);
      mSSA.attach_surface_model(surface);
      mSSA.attach_ocean_model(ocean);

      ierr = mSSA.init(); CHKERRQ(ierr);

      ierr = mSSA.diagnosticRun(); CHKERRQ(ierr);
      ierr = verbPrintf(2,com, "done with diagnostic run\n"); CHKERRQ(ierr);
      if (dontReport == PETSC_FALSE) {
        ierr = mSSA.reportErrors();  CHKERRQ(ierr);
      }
      ierr = mSSA.writeFiles("verify.nc"); CHKERRQ(ierr);
    } else { // run derived class for compensatory source SIA solutions
             // (i.e. compensatory accumulation or compensatory heating)
      IceCompModel mComp(g, config, overrides, test);
      ierr = mComp.setExecName("pismv"); CHKERRQ(ierr);
      mComp.attach_surface_model(surface);
      mComp.attach_ocean_model(ocean);

      ierr = mComp.init(); CHKERRQ(ierr);

      ierr = mComp.run(); CHKERRQ(ierr);
      ierr = verbPrintf(2,com, "done with run\n"); CHKERRQ(ierr);

      ThermoGlenArrIce*   tgaice = dynamic_cast<ThermoGlenArrIce*>(mComp.getIceFlowLaw());
      if (dontReport == PETSC_FALSE) {

        if (!IceFlowLawIsPatersonBuddCold(tgaice, config) && ((test == 'F') || (test == 'G'))) {
            ierr = verbPrintf(1,com, 
                "pismv WARNING: flow law must be cold part of Paterson-Budd ('-ice_type arr')\n"
                "   for reported errors in test %c to be meaningful!\n", test); CHKERRQ(ierr);
        }
        ierr = mComp.reportErrors();  CHKERRQ(ierr);
      }
      ierr = mComp.writeFiles("verify.nc"); CHKERRQ(ierr);
    }
    
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}


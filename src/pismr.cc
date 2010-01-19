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
  "Ice sheet driver for PISM ice sheet simulations, initialized from data.\n"
  "The basic PISM executable for evolution runs.\n";

#include <petsc.h>
#include "base/grid.hh"
#include "base/materials.hh"
#include "base/iceModel.hh"
#include "coupler/pccoupler.hh"
#include "coupler/pGreenlandAtmosCoupler.hh"

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
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    vector<string> required;
    required.clear(); // FIXME: either -boot_from or -i is required, but currently no
                      //        way to enforce (option a)|(option b)
                      // in fact, when -boot_from given then require
    ierr = show_usage_check_req_opts(com, "pismr", required,
      "  pismr {-i IN.nc|-boot_from IN.nc} [OTHER PISM & PETSc OPTIONS]\n\n"
      "where:\n"
      "  -i          input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_from  input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "notes:\n"
      "  * one of -i or -boot_from is required\n"
      "  * if -boot_from is used then in fact '-Mx A -My B -Mz C -Lz D' is also required\n"
      ); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "PISMR %s (basic evolution run mode)\n",
		      PISM_Revision); CHKERRQ(ierr);

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid g(com, rank, size);
    IceModel m(g, config, overrides);

    // Attach climate couplers:
    PISMConstAtmosCoupler     pcac;
    PISMGreenlandAtmosCoupler ppdd;
    PISMConstOceanCoupler     pcoc;
    PetscTruth  pddSet;
    ierr = check_option("-pdd", pddSet); CHKERRQ(ierr);
    if (pddSet == PETSC_TRUE) {
      ierr = verbPrintf(2,com, "pismr attaching PISMGreenlandAtmosCoupler to IceModel\n"); CHKERRQ(ierr);
      ierr = m.attachAtmospherePCC(ppdd); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,com, "pismr attaching PISMConstAtmosCoupler to IceModel\n"); CHKERRQ(ierr);
      ierr = m.attachAtmospherePCC(pcac); CHKERRQ(ierr);
    }
    ierr = m.attachOceanPCC(pcoc); CHKERRQ(ierr);

    ierr = m.setExecName("pismr"); CHKERRQ(ierr);

    ierr = m.init(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "running ...\n"); CHKERRQ(ierr);
    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);

    // provide a default output file name if no -o option is given.
    ierr = m.writeFiles("unnamed.nc"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

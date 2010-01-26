// Copyright (C) 2007-2010 Ed Bueler and Nathan Shemonski and Constantine Khroulev
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
  "Driver for PISM runs of EISMINT-Greenland intercomparison.\n";

#include <petsc.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "../coupler/pccoupler.hh"
#include "iceGRNModel.hh"

int main(int argc, char *argv[]){
  PetscErrorCode ierr;

  MPI_Comm com;
  PetscMPIInt rank, size;
  
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);
  
  { // explicit scoping does destructors before PetscFinalize() 
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    PetscTruth iset, bfset;
    ierr = check_option("-i", iset); CHKERRQ(ierr);
    ierr = check_option("-boot_from", bfset); CHKERRQ(ierr);
    string usage =
      "  pgrn {-i IN.nc|-boot_from IN.nc} [OTHER PISM & PETSc OPTIONS]\n\n"
      "where:\n"
      "  -i          input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_from  input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "notes:\n"
      "  * special executable for EISMINT-Greenland\n"
      "  * one of -i or -boot_from is required\n"
      "  * if -boot_from is used then in fact '-Mx A -My B -Mz C -Lz D' is also required\n"
      "  * generally behaves like pismr after initialization\n";
    if ((iset == PETSC_FALSE) && (bfset == PETSC_FALSE)) {
      ierr = PetscPrintf(com,
         "PISM ERROR: one of options -i,-boot_from is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "pgrn", usage.c_str()); CHKERRQ(ierr);
    } else {
      vector<string> required;  required.clear();
      ierr = show_usage_check_req_opts(com, "pgrn", required, usage.c_str()); CHKERRQ(ierr);
    }


    ierr = verbPrintf(2, com, "PGRN %s (PISM EISMINT-Greenland mode)\n",
		      PISM_Revision); CHKERRQ(ierr);

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid g(com, rank, size);
    IceGRNModel    m(g, config, overrides);
    ierr = m.setExecName("pgrn"); CHKERRQ(ierr);

    EISGREENAtmosCoupler   pegac;
    PISMConstOceanCoupler  pcoc;
    ierr = m.attachAtmospherePCC(pegac); CHKERRQ(ierr);
    ierr = m.attachOceanPCC(pcoc); CHKERRQ(ierr);
 
    ierr = m.init(); CHKERRQ(ierr);
 
    ierr = m.run(); CHKERRQ(ierr);
    ierr = verbPrintf(2, com, "done with run ... \n"); CHKERRQ(ierr);

    ierr = m.writeFiles("grn_exper.nc"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}


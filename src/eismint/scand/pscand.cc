// Copyright (C) 2009-2010 Andy Aschwandend and Ed Bueler
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
  "Ice sheet driver for Scandinavian modification of EISMINT II.\n";

#include <cstring>
#include <petscbag.h>
#include "../../base/grid.hh"
#include "../../base/materials.hh"
#include "iceScandModel.hh"

#include "../../coupler/PISMSurface.hh"
#include "../../coupler/PISMOcean.hh"

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

    PetscTruth  pddSet;
    ierr = check_option("-pdd", pddSet); CHKERRQ(ierr);
    if (pddSet == PETSC_TRUE) {
      ierr = PetscPrintf(com,
        "PISM ERROR: -pdd is not currently allowed as option to pscand\n");
        CHKERRQ(ierr);
      PetscEnd();
    }

    ierr = verbPrintf(2,com, "PSCAND %s (Scandinavian mod of EISMINT II mode)\n",
		      PISM_Revision); CHKERRQ(ierr);

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    // actually construct the IceModel
    IceGrid               g(com, rank, size);
    IceScandModel         m(g, config, overrides);

    // Create boundary models:
    PISMSurfaceModel *surface = new PSDummy(g, config);
    PISMOceanModel   *ocean   = new POConstant(g, config);;

    m.attach_ocean_model(ocean);
    m.attach_surface_model(surface);

    ierr = m.init(); CHKERRQ(ierr);

    ierr = m.setExecName("pscand"); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "running ...\n"); CHKERRQ(ierr);
    ierr = m.run(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "... done with run \n"); CHKERRQ(ierr);
    ierr = m.writeFiles("scand_exper.nc"); CHKERRQ(ierr);

  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}


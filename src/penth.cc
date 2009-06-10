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
  "TEMPORARY: USES IceEnthalpyModel.\n"
  "Ice sheet driver for PISM ice sheet simulations, initialized from data.\n"
  "The basic PISM executable for evolution runs.\n";

/*

suggested test procedure:  Use EISMINT II experiment A for only 20ka, so we do
SIA only and no sliding, just to see effect of corrected conservation of energy on flow:

pisms -eisII A -Mx 61 -My 61 -Mz 101 -quadZ -y 0.1 -o foo.nc
mpiexec -n 4 pismr -i foo.nc -y 19999.9 -o coldice.nc >> out.cold &
mpiexec -n 4 penth -i foo.nc -y 19999.9 -o polyice.nc >> out.poly &

*/

#include <petscvec.h>
#include "base/grid.hh"
#include "base/materials.hh"
#include "base/iceEnthalpyModel.hh"
#include "coupler/pccoupler.hh"

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
    IceGrid g(com, rank, size);
    PISMConstAtmosCoupler     pcac;
    PISMSnowModelAtmosCoupler ppdd;
    PISMConstOceanCoupler     pcoc;

    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(1,com, "PENTH %s (TEMPORARY ENTHALPY basic evolution run mode)\n",
		      PISM_Revision); CHKERRQ(ierr);

    IceEnthalpyModel m(g);
    ierr = m.setExecName("penth"); CHKERRQ(ierr);

    PetscTruth  pddSet;
    ierr = check_option("-pdd", pddSet); CHKERRQ(ierr);
    if (pddSet == PETSC_TRUE) {
      ierr = verbPrintf(2,com, "penth attaching PISMSnowModelAtmosCoupler to IceModel\n"); CHKERRQ(ierr);
      ierr = m.attachAtmospherePCC(ppdd); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,com, "penth attaching PISMConstAtmosCoupler to IceModel\n"); CHKERRQ(ierr);
      ierr = m.attachAtmospherePCC(pcac); CHKERRQ(ierr);
    }
    ierr = m.attachOceanPCC(pcoc); CHKERRQ(ierr);

    ierr = m.init(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "running ...\n"); CHKERRQ(ierr);
    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);

    // provide a default output file name if no -o option is given.
    ierr = m.writeFiles("unnamedEnth.nc"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}


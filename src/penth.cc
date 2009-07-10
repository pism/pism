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

one possible test procedure:  Use EISMINT II experiment A for only 20ka, so we do
SIA only and no sliding, just to see effect of corrected conservation of energy on flow
in case with nontrivial thickness:

mpiexec -n $NN pisms -temp_pa -eisII A -Mx 61 -My 61 -Mz 101 -quadZ -y 6000.0 -o nobr_estart.nc
mpiexec -n $NN pismr -temp_pa -i nobr_estart.nc -y 14000 -skip 10 -o nobr_coldice.nc
mpiexec -n $NN penth -i nobr_estart.nc -y 14000 -skip 10 -o nobr_polyice.nc

adding bedrock thermal:

mpiexec -n $NN pisms -temp_pa -eisII A -Mx 61 -My 61 -Mz 101 -Mbz 51 -quadZ -y 6000.0 -o estart.nc
mpiexec -n $NN pismr -temp_pa -i estart.nc -y 14000 -skip 10 -o coldice.nc
mpiexec -n $NN penth -i estart.nc -y 14000 -skip 10 -o polyice.nc

also, here is an example of regridding enthalpy, with 'y' flag in -regrid_vars:

mpiexec -n 2 penth -boot_from estart.nc -Mx 121 -My 121 -Mz 101 -Mbz 51 -quadZ -Lz 5000 -regrid_from polyice.nc -regrid_vars THLey -y 1000 -o finepolyice.nc >> out.finepoly &

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
    ierr = verbPrintf(1,com, "PENTH %s (development of ENTHALPY basic evolution run mode)\n",
		      PISM_Revision); CHKERRQ(ierr);

    IceEnthalpyModel m(g);
    ierr = m.setExecName("penth"); CHKERRQ(ierr);

    PetscTruth  pddSet;
    ierr = check_option("-pdd", pddSet); CHKERRQ(ierr);
    if (pddSet == PETSC_TRUE) {
      ierr = verbPrintf(2,com, 
        "penth attaching PISMSnowModelAtmosCoupler to IceEnthalpyModel\n"); CHKERRQ(ierr);
      ierr = m.attachAtmospherePCC(ppdd); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,com, 
        "penth attaching PISMConstAtmosCoupler to IceEnthalpyModel\n"); CHKERRQ(ierr);
      ierr = m.attachAtmospherePCC(pcac); CHKERRQ(ierr);
    }

    ierr = verbPrintf(2,com, 
        "penth attaching PISMConstOceanCoupler to IceEnthalpyModel\n"); CHKERRQ(ierr);
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


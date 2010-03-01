// Copyright (C) 2004-2010 Ed Bueler, Andy Aschwanden, and Constantine Khroulev
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

#include <petsc.h>
#include "base/grid.hh"
#include "base/iceEnthalpyModel.hh"

#include "coupler/PCFactory.hh"
#include "coupler/PISMAtmosphere.hh"
#include "coupler/PISMSurface.hh"
#include "coupler/PISMOcean.hh"

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

    ierr = verbPrintf(2,com, "PENTH %s (development of ENTHALPY basic evolution run mode)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    vector<string> required;
    required.clear();
    ierr = show_usage_check_req_opts(com, "penth", required,
      "  penth IS UNDER DEVELOPMENT AND WILL MERGE WITH pismr\n\n"
      "  SEE 'pismr -usage'\n"
      ); CHKERRQ(ierr);

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid g(com, rank, size);
    IceEnthalpyModel m(g, config, overrides);
    ierr = m.setExecName("penth"); CHKERRQ(ierr);

    // Initialize boundary models:
    PAFactory pa(g, config);
    PISMAtmosphereModel *atmosphere;

    PSFactory ps(g, config);
    PISMSurfaceModel *surface;

    POFactory po(g, config);
    PISMOceanModel *ocean;

    ierr = PetscOptionsBegin(com, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);
    pa.create(atmosphere);
    ps.create(surface);
    po.create(ocean);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    surface->attach_atmosphere_model(atmosphere);

    m.attach_ocean_model(ocean);
    m.attach_surface_model(surface);

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


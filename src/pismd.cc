// Copyright (C) 2007--2010 Ed Bueler and Constantine Khroulev
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

#include <petsc.h>
#include "base/grid.hh"
#include "base/iceModel.hh"
#include "eismint/iceROSSModel.hh"

#include "coupler/PCFactory.hh"
#include "coupler/PISMAtmosphere.hh"
#include "coupler/PISMSurface.hh"
#include "coupler/PISMOcean.hh"

static char help[] =
  "Driver for ice sheet, shelf, and stream simulations, for 'diagnostic'\n"
  "computation of velocity field from geometry and temperature field.\n"
  "(Also a driver for EISMINT-Ross diagnostic velocity computation in ice shelf.)\n";

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);
  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  { /* This explicit scoping forces destructors to be called before PetscFinalize() */
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "PISMD %s (diagnostic velocity computation mode)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    PetscTruth iset, bfset;
    ierr = check_option("-i", iset); CHKERRQ(ierr);
    ierr = check_option("-boot_from", bfset); CHKERRQ(ierr);
    string usage =
      "  pismd IS DEPRECATED\n\n"
      "  INTENDED REPLACEMENT IS 'pismr -y 0 -f3d '\n\n"
      "  SEE 'pismr -usage'\n";
    if ((iset == PETSC_FALSE) && (bfset == PETSC_FALSE)) {
      ierr = PetscPrintf(com,
         "PISM ERROR: one of options -i,-boot_from is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "pismd", usage.c_str()); CHKERRQ(ierr);
    } else {
      vector<string> required;  required.clear();
      ierr = show_usage_check_req_opts(com, "pismd", required, usage.c_str()); CHKERRQ(ierr);
    }

    // re this option, see  src/eismint/iceROSSModel.hh|cc and:
    //     D. MacAyeal and five others (1996). "An ice-shelf model test based on the 
    //     Ross ice shelf," Ann. Glaciol. 23, 46--51
    PetscTruth  doRoss;
    ierr = check_option("-ross", doRoss); CHKERRQ(ierr);

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid    g(com, rank, size);

    // Initialize boundary models:
    PAFactory pa(g, config);
    PISMAtmosphereModel *atmosphere;

    PSFactory ps(g, config);
    PISMSurfaceModel *surface;
    ierr = ps.set_default("constant"); CHKERRQ(ierr);

    POFactory po(g, config);
    PISMOceanModel *ocean;

    pa.create(atmosphere);
    ps.create(surface);
    po.create(ocean);

    surface->attach_atmosphere_model(atmosphere);

    IceModel *m;
    if (doRoss == PETSC_TRUE)
      m = new IceROSSModel(g, config, overrides);
    else 
      m = new IceModel(g, config, overrides);

    m->attach_surface_model(surface);
    m->attach_ocean_model(ocean);

    ierr = m->setExecName("pismd"); CHKERRQ(ierr);

    config.set_flag("force_full_diagnostics", true);

    ierr = m->init(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "computing velocity field (diagnostically) ...\n"); CHKERRQ(ierr);
    ierr = m->diagnosticRun(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "... done\n"); CHKERRQ(ierr);

    if (doRoss == PETSC_TRUE) {
      IceROSSModel* mRoss = dynamic_cast<IceROSSModel*>(m);
      if (!mRoss) { SETERRQ(1, "PISMD: ross finish files ... how did I get here?"); }
      ierr = mRoss->finishROSS(); CHKERRQ(ierr);
    }

    ierr = m->writeFiles("unnamed_diag.nc"); CHKERRQ(ierr);  // default filename if no -o
    
    delete m;
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

// Copyright (C) 2004-2012 Jed Brown, Ed Bueler and Constantine Khroulev
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
  "Ice sheet driver for EISMINT II, and other constant climate, simplified geometry\n"
  "intercomparison simulations.\n";

#include <cstring>
#include <petsc.h>
#include "IceGrid.hh"
#include "iceModel.hh"
#include "eismint/iceEISModel.hh"
#include "eismint/icePSTexModel.hh"
#include "pism_options.hh"

#include "PSDummy.hh"
#include "POConstant.hh"

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

    ierr = verbPrintf(2,com, "PISMS %s (simplified geometry mode)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    vector<string> required;
    required.clear(); // no actually required options; "-eisII A" is default
    ierr = show_usage_check_req_opts(com, "pisms", required,
      "  pisms [-eisII x|-pst -xxx] [OTHER PISM & PETSc OPTIONS]\n"
      "where major option chooses type of simplified experiment:\n"
      "  -eisII x    choose EISMINT II experiment (x = A|B|C|D|E|F|G|H|I|J|K|L)\n"
      "  -pst -xxx   choose plastic till ice stream experiment; see Bueler & Brown (2009);\n"
      "              (-xxx = -P0A|-P0I|-P1|-P2|-P3|-P4)\n"
      "notes:\n"
      "  -pdd        not allowed (because PISMConstAtmosCoupler is always used)\n"
      ); CHKERRQ(ierr);

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides, true); CHKERRQ(ierr);

    bool EISIIchosen, PSTexchosen;
    /* This option determines the single character name of EISMINT II experiments:
    "-eisII F", for example. */
    ierr = PISMOptionsIsSet("-eisII", EISIIchosen); CHKERRQ(ierr);
    /* This option chooses Plastic till ice Stream with Thermocoupling experiment. */
    ierr = PISMOptionsIsSet("-pst", PSTexchosen); CHKERRQ(ierr);

    int  choiceSum = (int) EISIIchosen + (int) PSTexchosen;
    if (choiceSum > 1) {
      ierr = PetscPrintf(com,
         "PISM ERROR: pisms called with more than one simplified geometry experiment chosen\n");
         CHKERRQ(ierr);
      PISMEnd();
    }

    // actually construct the IceModel
    IceGrid g(com, rank, size, config);

    // Initialize boundary models (climate will always come from
    // intercomparison formulas):
    PISMSurfaceModel *surface = new PSDummy(g, config);
    PISMOceanModel *ocean = new POConstant(g, config);

    IceModel *m;
    if (PSTexchosen == PETSC_TRUE) {
      m = new IcePSTexModel(g, config, overrides);
    } else {
      m = new IceEISModel(g, config, overrides);
    }

    m->attach_surface_model(surface);
    m->attach_ocean_model(ocean);
    ierr = m->setExecName("pisms"); CHKERRQ(ierr);

    ierr = m->init(); CHKERRQ(ierr);

    ierr = m->run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run \n"); CHKERRQ(ierr);
    ierr = m->writeFiles("simp_exper.nc"); CHKERRQ(ierr);

    delete m;
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}


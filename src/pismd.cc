// Copyright (C) 2007, 2008, 2009 Ed Bueler and Constantine Khroulev
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
#include "base/materials.hh"
#include "base/iceModel.hh"
#include "coupler/pccoupler.hh"
#include "eismint/iceROSSModel.hh"

static char help[] =
  "Driver for ice sheet, shelf, and stream simulations, for 'diagnostic' computation\n"
  "of velocity field from geometry and temperature field.\n"
  "(Also a driver for EISMINT-Ross diagnostic velocity computation in ice shelf.)\n";

/* 
1.  example of diagnostic computation of velocities from saved model state
file, using only SIA:

  $ pisms -eisII A -Mx 61 -My 61 -Mz 101 -y 6000 -o foo.nc
  $ pisms -eisII A -i foo.nc -y 0.00001 -f3d -o bar.nc
  $ pismd -i foo.nc -o full_foo.nc
  $ ncdiff -O full_foo.nc bar.nc pismddiff.nc  # velocities nearly equal, not quite at margin; why?

2. see example of diagnostic computation of Ross ice shelf velocities in manual

*/

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);
  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  { /* This explicit scoping forces destructors to be called before PetscFinalize() */
    IceGrid    g(com, rank, size);
    IceType*   ice = PETSC_NULL;
    PISMConstAtmosCoupler pcac;
    PISMConstOceanCoupler pcoc;
    
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(1,com, "PISMD  (diagnostic velocity computation mode)\n"); CHKERRQ(ierr);
    ierr = userChoosesIceType(com, ice); CHKERRQ(ierr);  // allocates ice

    IceModel*      m;
    IceModel       mPlain(g, ice);
    IceROSSModel   mRoss(g, ice);

    // re this option, see  src/eismint/iceROSSModel.hh|cc and:
    //     D. MacAyeal and five others (1996). "An ice-shelf model test based on the 
    //     Ross ice shelf," Ann. Glaciol. 23, 46--51
    PetscTruth  doRoss;
    ierr = PetscOptionsHasName(PETSC_NULL, "-ross", &doRoss); CHKERRQ(ierr);

    if (doRoss == PETSC_TRUE) {
      m = (IceModel*) &mRoss;
    } else 
      m = (IceModel*) &mPlain;
    ierr = m->setExecName("pismd"); CHKERRQ(ierr);

    ierr = m->attachAtmospherePCC(pcac); CHKERRQ(ierr);
    ierr = m->attachOceanPCC(pcoc); CHKERRQ(ierr);

    ierr = m->setFromOptions(); CHKERRQ(ierr);
    ierr = m->initFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "computing velocity field (diagnostically) ...\n"); CHKERRQ(ierr);
    ierr = m->diagnosticRun(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "... done\n"); CHKERRQ(ierr);

    if (doRoss == PETSC_TRUE) {
      ierr = mRoss.finishROSS(); CHKERRQ(ierr);
    }

    // provide a default base name if no -o option.
    ierr = m->writeFiles("unnamed_diag.nc",PETSC_TRUE); CHKERRQ(ierr);

    delete ice;
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

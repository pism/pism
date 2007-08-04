// Copyright (C) 2007 Ed Bueler
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
  "Driver for ice sheet, shelf, and stream simulations, for 'diagnostic' computation\n"
  "of velocity field from geometry and temperature field.\n"
  "(Also a driver for EISMINT-Ross diagnostic velocity computation in ice shelf.)\n";

/* 
1.  example of diagnostic computation of velocities from saved model state
file, using only SIA:

  $ pisms -eisII A -Mx 61 -My 61 -Mz 101 -y 6000 -o foo
  $ pisms -eisII A -if foo.nc -y 0.00001 -o bar
  $ pismd -if foo.nc -mu_sliding 0.0 -no_bmr_in_vert -d Xc03 -kd 10 -pause 10 -o full_foo
  $ ncdiff -O -v cbar full_foo.nc bar.nc cbardiff.nc  # cbar nearly equal, not quite at margin; why?
  
2. generic example of diagnostic computation of velocities from bootstrap file foo_boot.nc;  assumes
this file includes boundary conditions shelf/stream equations; at end saves "allfields.nc":

  $ pismd -bif foo_boot.nc -shelfstreamBC foo_boot.nc -Mx 61 -My 61 -Mz 101 \
    -mv -d cnmu -pause 10 -verbose

3. example of diagnostic computation of Ross ice shelf velocities:

  $ eis_ross.py --prefix=eisROSS/    # creates eis_ross.nc from EISMINT-Ross ascii files
  $ pismd -ross -bif eis_ross.nc -shelfstreamBC eis_ross.nc -Mx 147 -My 147 -Mz 11 \
    -Lz 1000 -mv -d cnmu -pause 5 -verbose 5
*/

#include <petsc.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"
#include "iceROSSModel.hh"

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
    IceType*   ice;
    PetscInt   flowlawNumber = 0; // use Paterson-Budd by default
    
    ierr = verbPrintf(1,com, "PISMD (diagnositic velocity computation mode)\n"); CHKERRQ(ierr);
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);

    IceModel*      m;
    IceModel       mPlain(g, *ice);
    IceROSSModel   mRoss(g, *ice);

    PetscTruth  doRoss, ssBCset;
    char        ssBCfile[PETSC_MAX_PATH_LEN];
    // re next option, see:
    //     D. MacAyeal and five others (1996). "An ice-shelf model test based on the 
    //     Ross ice shelf," Ann. Glaciol. 23, 46--51
    ierr = PetscOptionsHasName(PETSC_NULL, "-ross", &doRoss); CHKERRQ(ierr);
    ierr = PetscOptionsGetString(PETSC_NULL, "-shelfstreamBC", ssBCfile,
                                 PETSC_MAX_PATH_LEN, &ssBCset); CHKERRQ(ierr);
    if (doRoss == PETSC_TRUE) {
      ierr = verbPrintf(2,com, 
            "PISMD: initializing EISMINT Ross ice shelf velocity computation ... \n"); CHKERRQ(ierr);
      m = (IceModel*) &mRoss;
    } else 
      m = (IceModel*) &mPlain;
    ierr = m->setFromOptions(); CHKERRQ(ierr);
    ierr = m->initFromOptions(); CHKERRQ(ierr);

    if (ssBCset == PETSC_TRUE) {
       ierr = verbPrintf(2, com, 
             "  attempting to read mask and boundary conditions (ubar,vbar) from file %s\n",
             ssBCfile);   CHKERRQ(ierr);
       ierr = m->readShelfStreamBCFromFile_netCDF(ssBCfile); CHKERRQ(ierr);
       ierr = verbPrintf(2, com, 
             "  done reading -shelfstreamBC file and setting boundary conditions\n");
             CHKERRQ(ierr);
    }

    ierr = verbPrintf(2,com, "computing velocity field (diagnostically) ...\n"); CHKERRQ(ierr);
    ierr = m->diagnosticRun(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "... done\n"); CHKERRQ(ierr);

    // provide a default base name if no -o option.
    ierr = m->writeFiles("allfields",PETSC_TRUE); CHKERRQ(ierr);

    if (doRoss == PETSC_TRUE) {
      ierr = mRoss.finishROSS(); CHKERRQ(ierr);
    }

    ierr = verbPrintf(2,com, " ... done.\n"); CHKERRQ(ierr);
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

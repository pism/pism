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

#include <petsc.h>
#include "base/grid.hh"
#include "base/materials.hh"
#include "base/iceModel.hh"
#include "eismint/iceROSSModel.hh"

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

4. example of multiprocessor and grid-refined Ross ice shelf:

  $ eis_ross.py --prefix=eisROSS/    # creates eis_ross.nc from EISMINT-Ross ascii files
  $ fill_missing.py -i eis_ross.nc --eps=4e-8 -v ubar,vbar -o eis_ross1.nc
  $ mpiexec -n 2 pismd -ross -bif eis_ross1.nc -shelfstreamBC eis_ross1.nc \
    -Mx 293 -My 293 -Mz 5 -Lz 1000 -mv -d cnUu -pause 5 -verbose 5 -display :0

the second step fixes the boundary condition data in ubar,vbar variables in eis_ross.nc so
that it is usable under grid refinement

5. tuning run for Ross ice shelf (tuning constant hardness):

  $ eis_ross.py --prefix=eisROSS/    # creates eis_ross.nc from EISMINT-Ross ascii files
  $ pismd -ross -tune 1.6e8,0.1e8,2.3e8 -bif eis_ross.nc -shelfstreamBC eis_ross.nc \
    -Mx 147 -My 147 -Mz 11 -Lz 1000 -mv -d cnmu -pause 5 -verbose 5
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
    IceType*   ice;
    PetscInt   flowlawNumber = 0; // use Paterson-Budd by default
    
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(1,com, "PISMD (diagnostic velocity computation mode)\n"); CHKERRQ(ierr);
    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);

    IceModel*      m;
    IceModel       mPlain(g, *ice);
    IceROSSModel   mRoss(g, *ice);

    // re this option, see  src/eismint/iceROSSModel.hh|cc and:
    //     D. MacAyeal and five others (1996). "An ice-shelf model test based on the 
    //     Ross ice shelf," Ann. Glaciol. 23, 46--51
    PetscTruth  doRoss;
    ierr = PetscOptionsHasName(PETSC_NULL, "-ross", &doRoss); CHKERRQ(ierr);

    if (doRoss == PETSC_TRUE) {
      m = (IceModel*) &mRoss;
    } else 
      m = (IceModel*) &mPlain;
    ierr = m->setFromOptions(); CHKERRQ(ierr);
    ierr = m->initFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "computing velocity field (diagnostically) ...\n"); CHKERRQ(ierr);
    ierr = m->diagnosticRun(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "... done\n"); CHKERRQ(ierr);

    // provide a default base name if no -o option.
    ierr = m->writeFiles("unnamed_diag",PETSC_TRUE); CHKERRQ(ierr);

    if (doRoss == PETSC_TRUE) {
      ierr = mRoss.finishROSS(); CHKERRQ(ierr);
    }

    ierr = verbPrintf(2,com, " ... done.\n"); CHKERRQ(ierr);
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

// Copyright (C) 2007-2008 Ed Bueler and Constantine Khroulev
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

#include <petscda.h>
#include "iceModel.hh"


// Put these in a new 'software testing' executable; call by:
//     ierr = m->testIceModelVec3(); CHKERRQ(ierr);
//     ierr = m->testIceModelVec3Bedrock(); CHKERRQ(ierr);


// test IceModelVec3: assuming this is called from pisms, try
//   pisms -eisII A -y 1 -Mz 11    # no errors when grid coincides; significant otherwise
//   pisms -eisII A -y 1 -Mz 101   # no errors when grid coincides; small otherwise
//   pisms -eisII A -y 1 -Mx 5 -My 5 -Mz 501   # no errors
//   pisms -eisII A -y 1 -Mx 5 -My 5 -Mz 500   # small errors 
//                                               (appropriate; from linear interpolation)
PetscErrorCode IceModel::testIceModelVec3()    {
  PetscErrorCode ierr;
  IceModelVec3 test3;

  ierr = test3.create(grid,"testMe",true); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,"\n\ntesting IceModelVec3; setting to constant %f",
                    60402.70804); CHKERRQ(ierr);
  ierr = test3.set(60402.70804); CHKERRQ(ierr);

  ierr = test3.begin_access(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,"\n\nIceModelVec3::getValZ() says value is %f",
                    test3.getValZ(grid.xs,grid.ys,0.0) ); CHKERRQ(ierr);
  ierr = test3.end_access(); CHKERRQ(ierr);

  ierr = test3.beginGhostComm(); CHKERRQ(ierr);
  ierr = test3.endGhostComm(); CHKERRQ(ierr);

  PetscScalar levels[10] = {0.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 
                            1500.0, 2000.0, 2500.0, grid.Lz};
  PetscScalar valsIN[10], valsOUT[10];
  ierr = test3.begin_access(); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
    "\n\ntesting IceModelVec3::setValColumnPL() and getValColumnPL()\n"); CHKERRQ(ierr);
  for (PetscInt k=0; k < 10; k++) {
    valsIN[k] = sin(levels[k]/1000.0);
  }
  ierr = test3.setValColumnPL(grid.xs, grid.ys, 10, levels, valsIN); CHKERRQ(ierr);
  ierr = test3.getValColumnPL(grid.xs, grid.ys, 10, levels, valsOUT); CHKERRQ(ierr);
  for (PetscInt k=0; k < 10; k++) {
    ierr = verbPrintf(1,grid.com,
        "   k=%d:   level=%7.2f   valsIN=%7.4f   valsOUT=%7.4f   |diff|=%5.4e\n",
        k,levels[k],valsIN[k],valsOUT[k],PetscAbs(valsIN[k]-valsOUT[k]) ); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
    "\n\ntesting IceModelVec3::setValColumnQUAD() and getValColumnQUAD()\n");
    CHKERRQ(ierr);
  ierr = test3.getValColumnQUAD(grid.xs, grid.ys, 10, levels, valsOUT); CHKERRQ(ierr);
  for (PetscInt k=0; k < 10; k++) {
    ierr = verbPrintf(1,grid.com,
       "   k=%d:   level=%7.2f   valsIN=%7.4f   valsOUT=%7.4f   |diff|=%5.4e\n",
       k,levels[k],valsIN[k],valsOUT[k],PetscAbs(valsIN[k]-valsOUT[k]) ); 
       CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = test3.end_access(); CHKERRQ(ierr);
  ierr = test3.destroy(); CHKERRQ(ierr);
  return 0;
}


// test IceModelVec3Bedrock: assuming this is called from pisms, try
//   pisms -eisII A -y 1 -Mz 11 -Mbz 11  # no errors when grid coincides;
//                                       #   significant otherwise
//   pisms -eisII A -y 1 -Mz 101 -Mbz 101 # no errors because grid coincides
//   pisms -eisII A -y 1 -Mz 102 -Mbz 102 # small errors (grid doesn't coincide)
// same story here
//   pisms -eisII A -y 1 -Mz 11 -Mbz 11 -quadZ
//   pisms -eisII A -y 1 -Mz 101 -Mbz 101 -quadZ
//   pisms -eisII A -y 1 -Mz 102 -Mbz 102 -quadZ
PetscErrorCode IceModel::testIceModelVec3Bedrock()    {
  PetscErrorCode ierr;
  IceModelVec3Bedrock test3b;

  ierr = verbPrintf(1,grid.com,
    "\nbedrock elevations are grid.zblevels[0,...,grid.Mbz-1]:\n"); CHKERRQ(ierr);
  for (PetscInt k=0; k < grid.Mbz; k++) {
    ierr = verbPrintf(1,grid.com," %.3f,",grid.zblevels[k]); CHKERRQ(ierr);
  }
    ierr = verbPrintf(1,grid.com,"\n\n"); CHKERRQ(ierr);


  ierr = test3b.create(grid,"testMe",false); CHKERRQ(ierr);

  ierr = test3b.begin_access(); CHKERRQ(ierr);
  const PetscScalar  tv = 60402.70804;
  ierr = verbPrintf(1,grid.com,
             "\ntesting IceModelVec3Bedrock\n\nsetting to constant %f",
             tv); CHKERRQ(ierr);
  ierr = test3b.setColumn(grid.xs,grid.ys,tv); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,"\n\ngetInternalColumn() says ... ");
            CHKERRQ(ierr);
  PetscScalar *valscOUT;
  ierr = test3b.getInternalColumn(grid.xs,grid.ys,&valscOUT); CHKERRQ(ierr);
  PetscScalar maxdiff = 0.0;
  for (PetscInt k=0; k < grid.Mbz; k++) {
    maxdiff = PetscMax(maxdiff,PetscAbs(tv-valscOUT[k]));
  }
  ierr = test3b.end_access(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, "max error is %5.4e\n\n", maxdiff); CHKERRQ(ierr);

  ierr = test3b.begin_access(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,
             "\ntesting setInternalColumn() and getInternalColumn ... "); CHKERRQ(ierr);
  PetscScalar *valsiIN, *valsiOUT;
  valsiIN = new PetscScalar[grid.Mbz];
  for (PetscInt k=0; k < grid.Mbz; k++) {
    valsiIN[k] = sin(grid.zblevels[k]/833.42342);
  }
  ierr = test3b.setInternalColumn(grid.xs,grid.ys,valsiIN); CHKERRQ(ierr);
  ierr = test3b.getInternalColumn(grid.xs,grid.ys,&valsiOUT); CHKERRQ(ierr);
  maxdiff = 0.0;
  for (PetscInt k=0; k < grid.Mbz; k++) {
    maxdiff = PetscMax(maxdiff,PetscAbs(valsiIN[k]-valsiOUT[k]));
  }
  delete [] valsiIN;
  ierr = test3b.end_access(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, "max error is %5.4e\n\n", maxdiff); CHKERRQ(ierr);


  const PetscInt N=10;
  // make levels[] increasing and covering interval [-Lbz,0.0]
  PetscScalar levels[N] = {1.0, 0.9, 0.5, 0.4, 0.32, 0.31, 0.3, 0.2, 0.1, 0.0};
  for (PetscInt k=0; k < N; k++) {
    levels[k] = - grid.Lbz * levels[k];  
  }

  ierr = test3b.begin_access(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,
    "\ntesting setValColumn() and getValColumn():\n"); CHKERRQ(ierr);
  PetscScalar valsIN[N], valsOUT[N];
  for (PetscInt k=0; k < N; k++) {
    valsIN[k] = sin(levels[k]/1000.0);
  }
  ierr = test3b.setValColumn(grid.xs, grid.ys, N, levels, valsIN); CHKERRQ(ierr);
  ierr = test3b.getValColumn(grid.xs, grid.ys, N, levels, valsOUT); CHKERRQ(ierr);
  for (PetscInt k=0; k < N; k++) {
    ierr = verbPrintf(1,grid.com,
        "   k=%d:   level=%10.2f   valsIN=%7.4f   valsOUT=%7.4f   |diff|=%5.4e\n",
        k,levels[k],valsIN[k],valsOUT[k],PetscAbs(valsIN[k]-valsOUT[k]) ); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = test3b.end_access(); CHKERRQ(ierr);
  ierr = test3b.destroy(); CHKERRQ(ierr);
  return 0;
}




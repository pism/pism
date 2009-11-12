// Copyright (C) 2007-2009 Ed Bueler and Constantine Khroulev
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
#include "iMtests.hh"
#include "../coupler/iceModelVec2T.hh"
#include "../base/nc_util.hh"
#include <vector>

//! Set grid defaults for a particular unit test.
PetscErrorCode IceUnitModel::set_grid_defaults() {
  grid.Mx  = grid.My = 3;
  grid.Mbz = 11;
  grid.Lbz = 1000.0;
  return 0;
}

//! Initialize variables for a particular unit test.
PetscErrorCode IceUnitModel::set_vars_from_options() {
  return 0;
}

//! Run an unit test.
PetscErrorCode IceUnitModel::run() {
  PetscErrorCode ierr;
  PetscTruth flag;

  ierr = check_option("-IceModelVec3", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = testIceModelVec3(); CHKERRQ(ierr);
  }

  ierr = check_option("-IceModelVec3Bedrock", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = testIceModelVec3Bedrock(); CHKERRQ(ierr);
  }

  ierr = check_option("-IceModelVec2T", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = testIceModelVec2T(); CHKERRQ(ierr);
  }

  return 0;
}

//! Write output files.
PetscErrorCode IceUnitModel::writeFiles(const char*) {
  // We don't write anything by default.
  return 0;
}

// test IceModelVec3: assuming this is called from pisms, try
//   pisms -eisII A -y 1 -Mz 11    # no errors when grid coincides; significant otherwise
//   pisms -eisII A -y 1 -Mz 101   # no errors when grid coincides; small otherwise
//   pisms -eisII A -y 1 -Mx 5 -My 5 -Mz 501   # no errors
//   pisms -eisII A -y 1 -Mx 5 -My 5 -Mz 500   # small errors 
//                                               (appropriate; from linear interpolation)
PetscErrorCode IceUnitModel::testIceModelVec3()    {
  PetscErrorCode ierr;

  ierr = verbPrintf(1,grid.com,"\n\ntesting IceModelVec3; setting to constant %f",
                    60402.70804); CHKERRQ(ierr);
  ierr = T3.set(60402.70804); CHKERRQ(ierr);

  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,"\n\nIceModelVec3::getValZ() says value is %f",
                    T3.getValZ(grid.xs,grid.ys,0.0) ); CHKERRQ(ierr);
  ierr = T3.end_access(); CHKERRQ(ierr);

  ierr = T3.beginGhostComm(); CHKERRQ(ierr);
  ierr = T3.endGhostComm(); CHKERRQ(ierr);

  PetscScalar levels[10] = {0.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 
                            1500.0, 2000.0, 2500.0, grid.Lz};
  PetscScalar valsIN[10], valsOUT[10];
  ierr = T3.begin_access(); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
    "\n\ntesting IceModelVec3::setValColumnPL() and getValColumnPL()\n"); CHKERRQ(ierr);
  for (PetscInt k=0; k < 10; k++) {
    valsIN[k] = sin(levels[k]/1000.0);
  }
  ierr = T3.setValColumnPL(grid.xs, grid.ys, 10, levels, valsIN); CHKERRQ(ierr);
  ierr = T3.getValColumnPL(grid.xs, grid.ys, 10, levels, valsOUT); CHKERRQ(ierr);
  for (PetscInt k=0; k < 10; k++) {
    ierr = verbPrintf(1,grid.com,
        "   k=%d:   level=%7.2f   valsIN=%7.4f   valsOUT=%7.4f   |diff|=%5.4e\n",
        k,levels[k],valsIN[k],valsOUT[k],PetscAbs(valsIN[k]-valsOUT[k]) ); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
    "\n\ntesting IceModelVec3::setValColumnPL() and getValColumnQUAD()\n");
    CHKERRQ(ierr);
  ierr = T3.getValColumnQUAD(grid.xs, grid.ys, 10, levels, valsOUT); CHKERRQ(ierr);
  for (PetscInt k=0; k < 10; k++) {
    ierr = verbPrintf(1,grid.com,
       "   k=%d:   level=%7.2f   valsIN=%7.4f   valsOUT=%7.4f   |diff|=%5.4e\n",
       k,levels[k],valsIN[k],valsOUT[k],PetscAbs(valsIN[k]-valsOUT[k]) ); 
       CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = T3.end_access(); CHKERRQ(ierr);
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
PetscErrorCode IceUnitModel::testIceModelVec3Bedrock()    {
  PetscErrorCode ierr;

  ierr = verbPrintf(1,grid.com,
    "\nbedrock elevations are grid.zblevels[0,...,grid.Mbz-1]:\n"); CHKERRQ(ierr);
  for (PetscInt k=0; k < grid.Mbz; k++) {
    ierr = verbPrintf(1,grid.com," %.3f,",grid.zblevels[k]); CHKERRQ(ierr);
  }
    ierr = verbPrintf(1,grid.com,"\n\n"); CHKERRQ(ierr);

  ierr = Tb3.begin_access(); CHKERRQ(ierr);
  const PetscScalar  tv = 60402.70804;
  ierr = verbPrintf(1,grid.com,
             "\ntesting IceModelVec3Bedrock\n\nsetting to constant %f",
             tv); CHKERRQ(ierr);
  ierr = Tb3.setColumn(grid.xs,grid.ys,tv); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,"\n\ngetInternalColumn() says ... ");
            CHKERRQ(ierr);
  PetscScalar *valscOUT;
  ierr = Tb3.getInternalColumn(grid.xs,grid.ys,&valscOUT); CHKERRQ(ierr);
  PetscScalar maxdiff = 0.0;
  for (PetscInt k=0; k < grid.Mbz; k++) {
    maxdiff = PetscMax(maxdiff,PetscAbs(tv-valscOUT[k]));
  }
  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, "max error is %5.4e\n\n", maxdiff); CHKERRQ(ierr);

  ierr = Tb3.begin_access(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,
             "\ntesting setInternalColumn() and getInternalColumn ... "); CHKERRQ(ierr);
  PetscScalar *valsiIN, *valsiOUT;
  valsiIN = new PetscScalar[grid.Mbz];
  for (PetscInt k=0; k < grid.Mbz; k++) {
    valsiIN[k] = sin(grid.zblevels[k]/833.42342);
  }
  ierr = Tb3.setInternalColumn(grid.xs,grid.ys,valsiIN); CHKERRQ(ierr);
  ierr = Tb3.getInternalColumn(grid.xs,grid.ys,&valsiOUT); CHKERRQ(ierr);
  maxdiff = 0.0;
  for (PetscInt k=0; k < grid.Mbz; k++) {
    maxdiff = PetscMax(maxdiff,PetscAbs(valsiIN[k]-valsiOUT[k]));
  }
  delete [] valsiIN;
  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, "max error is %5.4e\n\n", maxdiff); CHKERRQ(ierr);


  const PetscInt N=10;
  // make levels[] increasing and covering interval [-Lbz,0.0]
  PetscScalar levels[N] = {1.0, 0.9, 0.5, 0.4, 0.32, 0.31, 0.3, 0.2, 0.1, 0.0};
  for (PetscInt k=0; k < N; k++) {
    levels[k] = - grid.Lbz * levels[k];  
  }

  ierr = Tb3.begin_access(); CHKERRQ(ierr);
  PetscScalar valsIN[N], valsOUT[N];
  for (PetscInt k=0; k < N; k++) {
    valsIN[k] = sin(levels[k]/1000.0);
  }

  ierr = Tb3.setValColumnPL(grid.xs, grid.ys, N, levels, valsIN); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
    "\ntesting setValColumn() and getValColumnPL():\n"); CHKERRQ(ierr);

  ierr = Tb3.getValColumnPL(grid.xs, grid.ys, N, levels, valsOUT); CHKERRQ(ierr);
  for (PetscInt k=0; k < N; k++) {
    ierr = verbPrintf(1,grid.com,
        "   k=%d:   level=%10.2f   valsIN=%7.4f   valsOUT=%7.4f   |diff|=%5.4e\n",
        k,levels[k],valsIN[k],valsOUT[k],PetscAbs(valsIN[k]-valsOUT[k]) ); CHKERRQ(ierr);
  }

  ierr = verbPrintf(1,grid.com,
    "\ntesting setValColumn() and getValColumnQUAD():\n"); CHKERRQ(ierr);

  ierr = Tb3.getValColumnQUAD(grid.xs, grid.ys, N, levels, valsOUT); CHKERRQ(ierr);
  for (PetscInt k=0; k < N; k++) {
    ierr = verbPrintf(1,grid.com,
        "   k=%d:   level=%10.2f   valsIN=%7.4f   valsOUT=%7.4f   |diff|=%5.4e\n",
        k,levels[k],valsIN[k],valsOUT[k],PetscAbs(valsIN[k]-valsOUT[k]) ); CHKERRQ(ierr);
  }

  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = Tb3.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceUnitModel::testIceModelVec2T() {
  PetscErrorCode ierr;
  NCTool nc(&grid);
  IceModelVec2T v;
  char filename[] = "test_IceModelVec2T.nc";

  // create a file to regrid from (that will have grid parameters compatible
  // with the current grid):
  ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attributes.write(filename); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  
  double t = 0, t_max = 50, dt = 0.35;
  while (t < t_max) {
    ierr = nc.open_for_writing(filename, true, true); CHKERRQ(ierr);
    ierr = nc.append_time(t); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = vH.set(t*t); CHKERRQ(ierr);
    ierr = vH.write(filename);
    t = t + dt;
  }

  ierr = v.create(grid, "thk", config.get("climate_forcing_buffer_size")); CHKERRQ(ierr);
  ierr = v.set_attrs("test", "IceModelVec2T test, using 'thk'", "m", ""); CHKERRQ(ierr);
  ierr = v.init(filename); CHKERRQ(ierr);

  double T = 1;
  double max_dt = 100;		// a big number
  ierr = v.max_timestep(T, max_dt); CHKERRQ(ierr);

  PetscPrintf(grid.com, "   max_dt = %f\n", max_dt);

  ierr = v.update(T, max_dt); CHKERRQ(ierr);

  int N = 11;
  vector<double> ts(N), values(N);
  dt = floor(max_dt) / (N - 1);

  for (int j = 0; j < N; ++j) {
    ts[j] = T + dt * j;
  }

  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = v.interp(grid.xs, grid.ys, N, &ts[0], &values[0]); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);

  for (int j = 0; j < N; ++j) {
    PetscPrintf(grid.com, "   value(%3.3f) = %f\n", ts[j], values[j]);
  }
  
  T += 3;
  ierr = v.update(T, max_dt); CHKERRQ(ierr);
  for (int j = 0; j < N; ++j) {
    ts[j] = T + dt * j;
  }

  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = v.interp(grid.xs, grid.ys, N, &ts[0], &values[0]); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);

  for (int j = 0; j < N; ++j) {
    PetscPrintf(grid.com, "   value(%3.3f) = %f\n", ts[j], values[j]);
  }

  char output[] = "test_output.nc";

  T = T + max_dt / 2.0;

  ierr = nc.open_for_writing(output, false, true); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attributes.write(filename); CHKERRQ(ierr);
  ierr = nc.append_time(T); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = v.update(T, 0); CHKERRQ(ierr);
  ierr = v.write(output, T, NC_DOUBLE); CHKERRQ(ierr);

  return 0;
}


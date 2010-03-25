// Copyright (C) 2007-2010 Ed Bueler and Constantine Khroulev
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
#include "../base/PISMIO.hh"
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

PetscErrorCode IceUnitModel::model_state_setup() {
 PetscErrorCode ierr;
  bool flag1, flag2;

  ierr = PISMOptionsIsSet("-dof1", flag1); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-dof2", flag2); CHKERRQ(ierr);
  if (flag1 || flag2) {
    // do nothing
  } else {
    ierr = IceModel::model_state_setup(); CHKERRQ(ierr);
  } 

  return 0;
}

PetscErrorCode IceUnitModel::createVecs() {
  PetscErrorCode ierr;
  bool flag1, flag2;

  ierr = PISMOptionsIsSet("-dof1", flag1); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-dof2", flag2); CHKERRQ(ierr);
  if (flag1 || flag2) {
    // create vH because the ocean coupler needs it:
    ierr = vH.create(grid, "thk", true); CHKERRQ(ierr);
    ierr = vH.set_attrs("model_state", "land ice thickness",
			"m", "land_ice_thickness"); CHKERRQ(ierr);
    ierr = vH.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    ierr = variables.add(vH); CHKERRQ(ierr);
  } else {
    ierr = IceModel::createVecs(); CHKERRQ(ierr);
  }


  return 0;
}

//! Run an unit test.
PetscErrorCode IceUnitModel::run() {
  PetscErrorCode ierr;
  bool flag;

  ierr = PISMOptionsIsSet("-IceModelVec3", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_IceModelVec3(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-IceModelVec3Bedrock", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_IceModelVec3Bedrock(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-IceModelVec2T", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_IceModelVec2T(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-IceModelVec2V", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_IceModelVec2V(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-dof1", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_dof1comm(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-dof2", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_dof2comm(); CHKERRQ(ierr);
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
PetscErrorCode IceUnitModel::test_IceModelVec3()    {
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

  PetscScalar *values_in = new PetscScalar[grid.Mz_fine],
    *values_out = new PetscScalar[grid.Mz_fine];

  ierr = T3.begin_access(); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
    "\n\ntesting IceModelVec3::setValColumnPL() and getValColumnPL()\n"); CHKERRQ(ierr);
  for (PetscInt k=0; k < grid.Mz_fine; k++) {
    values_in[k] = sin(grid.zlevels_fine[k]/1000.0);
  }
  ierr = T3.setValColumnPL(grid.xs, grid.ys, values_in); CHKERRQ(ierr);
  ierr = T3.getValColumnPL(grid.xs, grid.ys, grid.Mz_fine, values_out); CHKERRQ(ierr);
  for (PetscInt k=0; k < grid.Mz_fine; k++) {
    ierr = verbPrintf(1,grid.com,
        "   k=%d:   level=%7.2f   values_in=%7.4f   values_out=%7.4f   |diff|=%5.4e\n",
        k,grid.zlevels_fine[k],values_in[k],values_out[k],PetscAbs(values_in[k]-values_out[k]) ); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
    "\n\ntesting IceModelVec3::setValColumnPL() and getValColumnQUAD()\n");
  CHKERRQ(ierr);
  ierr = T3.getValColumnQUAD(grid.xs, grid.ys, grid.Mz_fine, values_out); CHKERRQ(ierr);
  for (PetscInt k=0; k < grid.Mz_fine; k++) {
    ierr = verbPrintf(1,grid.com,
       "   k=%d:   level=%7.2f   values_in=%7.4f   values_out=%7.4f   |diff|=%5.4e\n",
       k,grid.zlevels_fine[k],values_in[k],values_out[k],PetscAbs(values_in[k]-values_out[k]) ); 
       CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = T3.end_access(); CHKERRQ(ierr);

  delete[] values_in;
  delete[] values_out;

  return 0;
}


// test IceModelVec3Bedrock: assuming this is called from pisms, try
//   pisms -eisII A -y 1 -Mz 11 -Mbz 11  # no errors when grid coincides;
//                                       #   significant otherwise
//   pisms -eisII A -y 1 -Mz 101 -Mbz 101 # no errors because grid coincides
//   pisms -eisII A -y 1 -Mz 102 -Mbz 102 # small errors (grid doesn't coincide)
// same story here
//   pisms -eisII A -y 1 -Mz 11 -Mbz 11
//   pisms -eisII A -y 1 -Mz 101 -Mbz 101
//   pisms -eisII A -y 1 -Mz 102 -Mbz 102
PetscErrorCode IceUnitModel::test_IceModelVec3Bedrock()    {
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
  PetscScalar *values_in = new PetscScalar[grid.Mbz_fine],
    *values_out = new PetscScalar[grid.Mbz_fine];

  for (PetscInt k=0; k < grid.Mbz_fine; k++) {
    values_in[k] = sin(grid.zblevels_fine[k]/1000.0);
  }

  ierr = Tb3.setValColumnPL(grid.xs, grid.ys, values_in); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
    "\ntesting setValColumn() and getValColumnPL():\n"); CHKERRQ(ierr);

  ierr = Tb3.getValColumnPL(grid.xs, grid.ys, values_out); CHKERRQ(ierr);
  for (PetscInt k=0; k < grid.Mbz_fine; k++) {
    ierr = verbPrintf(1,grid.com,
        "   k=%d:   level=%10.2f   values_in=%7.4f   values_out=%7.4f   |diff|=%5.4e\n",
        k,grid.zblevels_fine[k],values_in[k],values_out[k],PetscAbs(values_in[k]-values_out[k]) ); CHKERRQ(ierr);
  }

  ierr = verbPrintf(1,grid.com,
    "\ntesting setValColumn() and getValColumnQUAD():\n"); CHKERRQ(ierr);

  ierr = Tb3.getValColumnQUAD(grid.xs, grid.ys, values_out); CHKERRQ(ierr);
  for (PetscInt k=0; k < grid.Mbz_fine; k++) {
    ierr = verbPrintf(1,grid.com,
        "   k=%d:   level=%10.2f   values_in=%7.4f   values_out=%7.4f   |diff|=%5.4e\n",
        k,grid.zblevels_fine[k],values_in[k],values_out[k],PetscAbs(values_in[k]-values_out[k]) ); CHKERRQ(ierr);
  }

  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = Tb3.end_access(); CHKERRQ(ierr);

  delete[] values_in;
  delete[] values_out;
  return 0;
}

PetscErrorCode IceUnitModel::test_IceModelVec2T() {
  PetscErrorCode ierr;
  PISMIO nc(&grid);
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

  v.set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
  ierr = v.create(grid, "thk", false); CHKERRQ(ierr);
  ierr = v.set_attrs("test", "IceModelVec2T test, using 'thk'", "m", ""); CHKERRQ(ierr);
  ierr = v.init(filename); CHKERRQ(ierr);

  double T = 1;
  double max_dt;
  max_dt = v.max_timestep(T);

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
    PetscPrintf(grid.com, "   (%3.3f)^2 ~= %f\n", ts[j], values[j]);
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
    PetscPrintf(grid.com, "   (%3.3f)^2 ~= %f\n", ts[j], values[j]);
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
  ierr = v.interp(T); CHKERRQ(ierr);
  ierr = v.write(output, NC_DOUBLE); CHKERRQ(ierr);

  T = 13;
  dt = 10;
  double average;

  ierr = v.update(T, dt); CHKERRQ(ierr);

  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = v.average(grid.xs, grid.ys, T, dt, average); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);

  PetscPrintf(grid.com,
	      "   average(%3.3f, %3.3f) = %3.3f, ideally should be 997/3 ~= 332.33333\n",
	      T, T + dt, average);

  return 0;
}

PetscErrorCode IceUnitModel::test_IceModelVec2V() {
  PetscErrorCode ierr;

  PISMIO nc(&grid);
  IceModelVec2V uvbar_ssa;

  ierr = uvbar_ssa.create(grid, "bar_ssa", true); CHKERRQ(ierr);

  ierr = uvbar_ssa.set_attrs("internal",
			      "x component of the SSA horizontal ice velocity",
			     "m s-1", "", 0); CHKERRQ(ierr);
  ierr = uvbar_ssa.set_attrs("internal",
			      "y component of the SSA horizontal ice velocity",
			     "m s-1", "", 1); CHKERRQ(ierr);
  ierr = uvbar_ssa.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  uvbar_ssa.write_in_glaciological_units = true;

  ierr = uvbar_ssa.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      uvbar_ssa(i,j).u = 3.0/secpera;
      uvbar_ssa(i,j).v = 4.0/secpera;
    }
  }
  ierr = uvbar_ssa.end_access(); CHKERRQ(ierr);

  // get and view components:
  ierr = uvbar_ssa.get_component(0, vubar); CHKERRQ(ierr);
  ierr = uvbar_ssa.get_component(1, vvbar); CHKERRQ(ierr);
  ierr = vubar.view(300); CHKERRQ(ierr);
  ierr = vvbar.view(300); CHKERRQ(ierr);
  PetscSleep(5);

  // write to a file:

  char filename[] = "test_IceModelVec2V.nc";
  // create a file to regrid from (that will have grid parameters compatible
  // with the current grid):
  ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attributes.write(filename); CHKERRQ(ierr);
  ierr = nc.append_time(0.0); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = uvbar_ssa.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  // reset:
  ierr = uvbar_ssa.set(0.0); CHKERRQ(ierr);

  // read in:
  ierr = uvbar_ssa.read(filename, 0); CHKERRQ(ierr);
  // write out:
  ierr = uvbar_ssa.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  ierr = uvbar_ssa.magnitude(vWork2d[0]); CHKERRQ(ierr);
  ierr = vWork2d[0].set_name("cbar"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs("diagnostic", 
			  "magnitude of vertically-integrated horizontal velocity of ice",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = vWork2d[0].set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vWork2d[0].write_in_glaciological_units = true;

  ierr = vWork2d[0].write(filename, NC_DOUBLE); CHKERRQ(ierr);

  return 0;
}

//! Copyright (C) PETSc authors.
//! copied from PETSc src/benchmarks/Index.c
static int BlastCache(void)
{
  int    i,ierr,n = 100000;
  PetscScalar *x,*y,*z,*a,*b;

  ierr = PetscMalloc(5*n*sizeof(PetscScalar),&x);CHKERRQ(ierr); // this will be ~4Mb...
  y = x + n;
  z = y + n;
  a = z + n;
  b = a + n;

  for (i=0; i<n; i++) {
    a[i] = (PetscScalar) i;
    y[i] = (PetscScalar) i;
    z[i] = (PetscScalar) i;
    b[i] = (PetscScalar) i;
    x[i] = (PetscScalar) i;
  }

  for (i=0; i<n; i++) {
    a[i] = 3.0*x[i] + 2.0*y[i] + 3.3*z[i] - 25.*b[i];
  }
  for (i=0; i<n; i++) {
    b[i] = 3.0*x[i] + 2.0*y[i] + 3.3*a[i] - 25.*b[i];
  }
  for (i=0; i<n; i++) {
    z[i] = 3.0*x[i] + 2.0*y[i] + 3.3*a[i] - 25.*b[i];
  }
  ierr = PetscFree(x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode IceUnitModel::test_dof1comm() {
  PetscErrorCode ierr;

  IceModelVec2S ubar, vbar;
  PetscScalar **u, **v;

  ierr = ubar.create(grid, "ubar", true); CHKERRQ(ierr);
  ierr = vbar.create(grid, "vbar", true); CHKERRQ(ierr);

  ierr = ubar.set(0.0); CHKERRQ(ierr);
  ierr = ubar.set(100.0); CHKERRQ(ierr);

  for (int k = 0; k < 100; ++k) {

    BlastCache();
    
    // set:
    ierr = ubar.get_array(u); CHKERRQ(ierr);
    ierr = vbar.get_array(v); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	u[i][j] = 0.25 * (v[i-1][j] + v[i+1][j] + v[i][j-1] + v[i][j+1]);
	v[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
      }
    }
    ierr = ubar.end_access(); CHKERRQ(ierr);
    ierr = vbar.end_access(); CHKERRQ(ierr);

    BlastCache();

    // communicate:
    ierr = ubar.beginGhostComm(); CHKERRQ(ierr);
    ierr = ubar.endGhostComm(); CHKERRQ(ierr);

    ierr = vbar.beginGhostComm(); CHKERRQ(ierr);
    ierr = vbar.endGhostComm(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceUnitModel::test_dof2comm() {
  PetscErrorCode ierr;

  IceModelVec2V uvbar_ssa;

  ierr = uvbar_ssa.create(grid, "bar_ssa", true); CHKERRQ(ierr);

  ierr = uvbar_ssa.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      uvbar_ssa(i,j).u = 0.0;
      uvbar_ssa(i,j).v = 100.0;
    }
  }
  ierr = uvbar_ssa.end_access(); CHKERRQ(ierr);

  for (int k = 0; k < 100; ++k) {

    BlastCache();

    // set:
    ierr = uvbar_ssa.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	uvbar_ssa(i,j).u = 0.25 * (uvbar_ssa(i-1,j).v + uvbar_ssa(i+1,j).v + uvbar_ssa(i,j-1).v + uvbar_ssa(i,j+1).v);
	uvbar_ssa(i,j).v = 0.25 * (uvbar_ssa(i-1,j).u + uvbar_ssa(i+1,j).u + uvbar_ssa(i,j-1).u + uvbar_ssa(i,j+1).u);
      }
    }
    ierr = uvbar_ssa.end_access(); CHKERRQ(ierr);

    BlastCache();

    // communicate:
    ierr = uvbar_ssa.beginGhostComm(); CHKERRQ(ierr);
    ierr = uvbar_ssa.endGhostComm(); CHKERRQ(ierr);
  }
  return 0;
}

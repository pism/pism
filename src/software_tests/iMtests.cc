// Copyright (C) 2007-2012 Ed Bueler and Constantine Khroulev
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

#include <vector>
#include <petscdmda.h>
#include "iMtests.hh"
#include "iceModelVec2T.hh"
#include "PIO.hh"
#include "PISMProf.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "PISMTime.hh"

//! Set grid defaults for a particular unit test.
PetscErrorCode IceUnitModel::set_grid_defaults() {
  grid.Mx  = grid.My = 3;
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

PetscErrorCode IceUnitModel::test_add_2d() {
  PetscErrorCode ierr;

  IceModelVec2V v1, v2, v3;

  ierr = v1.create(grid, "v1", false); CHKERRQ(ierr);
  ierr = v2.create(grid, "v2", true); CHKERRQ(ierr);
  ierr = v3.create(grid, "v3", true, 2); CHKERRQ(ierr);

  ierr = v1.set(0); CHKERRQ(ierr);

  ierr = v2.set(1); CHKERRQ(ierr);

  ierr = v3.set(2); CHKERRQ(ierr);

  // adding to a global vec, everything is done locally
  ierr = v1.add(1.0, v2); CHKERRQ(ierr);

  //  v1 = 1;
  ierr = v1.dump("v1.nc"); CHKERRQ(ierr);

  // adding to a local vec using a global one, have to scatter ghosts
  ierr = v2.add(1.0, v1); CHKERRQ(ierr);

  // v2 = 2;
  ierr = v2.dump("v2_1.nc"); CHKERRQ(ierr);

  // adding to a local vec using a different local vec, locally
  ierr = v2.add(1.0, v3); CHKERRQ(ierr);

  // v2 = 4
  ierr = v2.dump("v2_2.nc"); CHKERRQ(ierr);

  // adding to a local vec using a different local vec with a narrower stencil,
  // we have to scatter ghosts
  ierr = v3.add(1.0, v2); CHKERRQ(ierr);

  // v3 = 6
  ierr = v3.dump("v3.nc"); CHKERRQ(ierr);

  return 0;
}

//! Run an unit test.
PetscErrorCode IceUnitModel::run() {
  PetscErrorCode ierr;
  bool flag;

  ierr = PISMOptionsIsSet("-output", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_output(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-IceModelVec3", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_IceModelVec3(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-IceModelVec2T", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_IceModelVec2T(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-IceModelVec2V", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_IceModelVec2V(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-add_2d", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_add_2d(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-dof1", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_dof1comm(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-dof2", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_dof2comm(); CHKERRQ(ierr);
  }

  ierr = PISMOptionsIsSet("-prof", flag); CHKERRQ(ierr);
  if (flag) {
    ierr = test_pismprof(); CHKERRQ(ierr);
  }


  return 0;
}

//! Write output files.
PetscErrorCode IceUnitModel::writeFiles(string) {
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
  ierr = Enth3.set(60402.70804); CHKERRQ(ierr);

  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,"\n\nIceModelVec3::getValZ() says value is %f",
                    Enth3.getValZ(grid.xs,grid.ys,0.0) ); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);

  ierr = Enth3.beginGhostComm(); CHKERRQ(ierr);
  ierr = Enth3.endGhostComm(); CHKERRQ(ierr);

  PetscScalar *values_in = new PetscScalar[grid.Mz_fine],
    *values_out = new PetscScalar[grid.Mz_fine];

  ierr = Enth3.begin_access(); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
    "\n\ntesting IceModelVec3::setValColumnPL() and getValColumnPL()\n"); CHKERRQ(ierr);
  for (PetscInt k=0; k < grid.Mz_fine; k++) {
    values_in[k] = sin(grid.zlevels_fine[k]/1000.0);
  }
  ierr = Enth3.setValColumnPL(grid.xs, grid.ys, values_in); CHKERRQ(ierr);
  ierr = Enth3.getValColumnPL(grid.xs, grid.ys, grid.Mz_fine, values_out); CHKERRQ(ierr);
  for (PetscInt k=0; k < grid.Mz_fine; k++) {
    ierr = verbPrintf(1,grid.com,
        "   k=%d:   level=%7.2f   values_in=%7.4f   values_out=%7.4f   |diff|=%5.4e\n",
        k,grid.zlevels_fine[k],values_in[k],values_out[k],PetscAbs(values_in[k]-values_out[k]) ); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
    "\n\ntesting IceModelVec3::setValColumnPL() and getValColumnQUAD()\n");
  CHKERRQ(ierr);
  ierr = Enth3.getValColumnQUAD(grid.xs, grid.ys, grid.Mz_fine, values_out); CHKERRQ(ierr);
  for (PetscInt k=0; k < grid.Mz_fine; k++) {
    ierr = verbPrintf(1,grid.com,
       "   k=%d:   level=%7.2f   values_in=%7.4f   values_out=%7.4f   |diff|=%5.4e\n",
       k,grid.zlevels_fine[k],values_in[k],values_out[k],PetscAbs(values_in[k]-values_out[k]) ); 
       CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com,"done\n\n\n"); CHKERRQ(ierr);

  ierr = Enth3.end_access(); CHKERRQ(ierr);

  delete[] values_in;
  delete[] values_out;

  return 0;
}

PetscErrorCode IceUnitModel::test_output() {
  PetscErrorCode ierr;
  PetscScalar *E;

  PIO nc(grid.com, grid.rank, grid.config.get_string("output_format"));
  string filename = "test_output.nc";

  ierr = nc.open(filename, NC_WRITE); CHKERRQ(ierr);
  ierr = nc.def_time(config.get_string("time_dimension_name"),
                     config.get_string("calendar"),
                     grid.time->units()); CHKERRQ(ierr);
  ierr = nc.append_time(config.get_string("time_dimension_name"), grid.time->current()); CHKERRQ(ierr);
  ierr = nc.close();

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      vH(i,j) = grid.rank;

      ierr = Enth3.getInternalColumn(i, j, &E); CHKERRQ(ierr);

      for (int k = 0; k < grid.Mz; ++k) {
        E[k] = PetscAbs(grid.x[i] + grid.y[j]) + grid.zlevels[k];
      }
    }
  }

  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  bool no_2d = false, no_3d = false;
  ierr = PISMOptionsIsSet("-no_2d", no_2d); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-no_3d", no_3d); CHKERRQ(ierr);

  if (no_2d == false) {
    vH.time_independent = true;
    ierr = vH.write(filename); CHKERRQ(ierr);
  }

  if (no_3d == false) {
    Enth3.time_independent = true;
    ierr = Enth3.write(filename); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceUnitModel::test_IceModelVec2T() {
  PetscErrorCode ierr;
  PIO nc(grid.com, grid.rank, grid.config.get_string("output_format"));
  IceModelVec2T v;
  char filename[] = "test_IceModelVec2T.nc";

  // create a file to regrid from (that will have grid parameters compatible
  // with the current grid):
  ierr = nc.open(filename, NC_WRITE); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attributes.write(filename); CHKERRQ(ierr);

  double t = 0, t_max = 50, my_dt = 0.35;
  while (t < t_max) {
    ierr = nc.open(filename, NC_WRITE, true); CHKERRQ(ierr);
    ierr = nc.def_time(config.get_string("time_dimension_name"),
                       config.get_string("calendar"),
                       grid.time->units()); CHKERRQ(ierr);
    ierr = nc.append_time(config.get_string("time_dimension_name"), t); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = vH.set(t*t); CHKERRQ(ierr);
    ierr = vH.write(filename);
    t = t + my_dt;
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
  my_dt = floor(max_dt) / (N - 1);

  for (int j = 0; j < N; ++j) {
    ts[j] = T + my_dt * j;
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
    ts[j] = T + my_dt * j;
  }

  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = v.interp(grid.xs, grid.ys, N, &ts[0], &values[0]); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);

  for (int j = 0; j < N; ++j) {
    PetscPrintf(grid.com, "   (%3.3f)^2 ~= %f\n", ts[j], values[j]);
  }

  char output[] = "test_output.nc";

  T = T + max_dt / 2.0;

  ierr = nc.open(output, NC_WRITE); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attributes.write(filename); CHKERRQ(ierr);
  ierr = nc.def_time(config.get_string("time_dimension_name"),
                     config.get_string("calendar"),
                     grid.time->units()); CHKERRQ(ierr);
  ierr = nc.append_time(config.get_string("time_dimension_name"), T); CHKERRQ(ierr);

  ierr = v.update(T, 0); CHKERRQ(ierr);
  ierr = v.interp(T); CHKERRQ(ierr);
  ierr = v.write(output, NC_DOUBLE); CHKERRQ(ierr);

  T = 13;
  my_dt = 10;
  double average;

  ierr = v.update(T, my_dt); CHKERRQ(ierr);

  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = v.average(grid.xs, grid.ys, T, my_dt, average); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);

  PetscPrintf(grid.com,
	      "   average(%3.3f, %3.3f) = %3.3f, ideally should be 997/3 ~= 332.33333\n",
	      T, T + my_dt, average);

  return 0;
}

PetscErrorCode IceUnitModel::test_IceModelVec2V() {
  PetscErrorCode ierr;

  PIO nc(grid.com, grid.rank, grid.config.get_string("output_format"));
  IceModelVec2V velocity;

  ierr = velocity.create(grid, "bar", true); CHKERRQ(ierr);

  ierr = velocity.set_attrs("internal",
			      "x component of the horizontal ice velocity",
			     "m s-1", "", 0); CHKERRQ(ierr);
  ierr = velocity.set_attrs("internal",
			      "y component of the horizontal ice velocity",
			     "m s-1", "", 1); CHKERRQ(ierr);
  ierr = velocity.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  velocity.write_in_glaciological_units = true;

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      velocity(i,j).u = 3.0/secpera;
      velocity(i,j).v = 4.0/secpera;
    }
  }
  ierr = velocity.end_access(); CHKERRQ(ierr);

  // get and view components:
  ierr = velocity.get_component(0, vWork2d[0]); CHKERRQ(ierr);
  ierr = velocity.get_component(1, vWork2d[1]); CHKERRQ(ierr);
  ierr = vWork2d[0].view(300); CHKERRQ(ierr);
  ierr = vWork2d[1].view(300); CHKERRQ(ierr);
  PetscSleep(5);

  // write to a file:

  char filename[] = "test_IceModelVec2V.nc";
  // create a file to regrid from (that will have grid parameters compatible
  // with the current grid):
  ierr = nc.open(filename, NC_WRITE); CHKERRQ(ierr);
  ierr = nc.def_time(config.get_string("time_dimension_name"),
                     config.get_string("calendar"),
                     grid.time->units()); CHKERRQ(ierr);
  ierr = nc.append_time(config.get_string("time_dimension_name"), 0.0); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attributes.write(filename); CHKERRQ(ierr);

  ierr = velocity.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  // reset:
  ierr = velocity.set(0.0); CHKERRQ(ierr);

  // read in:
  ierr = velocity.read(filename, 0); CHKERRQ(ierr);
  // write out:
  ierr = velocity.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  ierr = velocity.magnitude(vWork2d[0]); CHKERRQ(ierr);
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

  IceModelVec2V velocity;

  ierr = velocity.create(grid, "bar", true); CHKERRQ(ierr);

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      velocity(i,j).u = 0.0;
      velocity(i,j).v = 100.0;
    }
  }
  ierr = velocity.end_access(); CHKERRQ(ierr);

  for (int k = 0; k < 100; ++k) {

    BlastCache();

    // set:
    ierr = velocity.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	velocity(i,j).u = 0.25 * (velocity(i-1,j).v + velocity(i+1,j).v + velocity(i,j-1).v + velocity(i,j+1).v);
	velocity(i,j).v = 0.25 * (velocity(i-1,j).u + velocity(i+1,j).u + velocity(i,j-1).u + velocity(i,j+1).u);
      }
    }
    ierr = velocity.end_access(); CHKERRQ(ierr);

    BlastCache();

    // communicate:
    ierr = velocity.beginGhostComm(); CHKERRQ(ierr);
    ierr = velocity.endGhostComm(); CHKERRQ(ierr);
  }
  return 0;
}

//! Test the PISM profiler class.
PetscErrorCode IceUnitModel::test_pismprof() {
  PetscErrorCode ierr;

  PISMProf *my_prof;

  my_prof = new PISMProf(grid.com, grid.rank, grid.size);
  my_prof->Nx = grid.Nx;
  my_prof->Ny = grid.Ny;

  int event_total, event_wait, event_sleep;
  event_total = my_prof->create("total",    "total");
  event_wait  = my_prof->create("waiting",  "waiting");
  event_sleep = my_prof->create("sleeping", "sleeping");

  my_prof->begin(event_total);
  {

    my_prof->begin(event_sleep);
    {
      fprintf(stderr, "Rank %d will sleep for %d seconds...\n", grid.rank, grid.rank);
      PetscSleep(grid.rank);
      fprintf(stderr, "Rank %d is up...\n", grid.rank);
    }
    my_prof->end(event_sleep);   

    my_prof->begin(event_wait);
    {
      ierr = my_prof->barrier(); CHKERRQ(ierr);
    }
    my_prof->end(event_wait);   

  }
  my_prof->end(event_total);   

  string output_name = "pism_out.nc";

  output_name = pism_filename_add_suffix(output_name, "-prof", "");

  ierr = my_prof->save_report(output_name); CHKERRQ(ierr); 

  delete my_prof;

  return 0;
}

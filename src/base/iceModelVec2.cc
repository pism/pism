// Copyright (C) 2008, 2009 Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cstdlib>
#include <petscda.h>
#include <netcdf.h>
#include "nc_util.hh"

#include "iceModelVec.hh"

// this file contains methods for derived classes IceModelVec2

// methods for base class IceModelVec are in "iceModelVec.cc"

IceModelVec2::IceModelVec2() : IceModelVec() {}


PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_short_name[], bool local) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2 with short_name='%s' already allocated\n",my_short_name);
  }
  PetscErrorCode ierr = create(my_grid, my_short_name, local, DA_STENCIL_BOX); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode  IceModelVec2::createSameDA(IceModelVec2 imv2_source,
					   IceGrid &my_grid, const char my_short_name[], bool local) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2 with short_name='%s' already allocated\n",my_short_name);
  }

  grid = &my_grid;
  dims = GRID_2D;

  da = imv2_source.da;
  IOwnDA = false;

  PetscErrorCode ierr;
  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  strcpy(short_name,my_short_name);
#ifdef PISM_DEBUG
  creation_counter += 1;
#endif // PISM_DEBUG

  return 0;
}

PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_short_name[], bool local,
                                     DAStencilType my_sten) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  grid = &my_grid;
  dims = GRID_2D;
  
  PetscInt       M, N, m, n;
  PetscErrorCode ierr;
  ierr = DAGetInfo(my_grid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate2d(my_grid.com, DA_XYPERIODIC, my_sten, N, M, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, &da); CHKERRQ(ierr);
  IOwnDA = true;

  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  strcpy(short_name,my_short_name);
#ifdef PISM_DEBUG
  creation_counter += 1;
#endif // PISM_DEBUG
  return 0;
}

//! Defines a netcdf variable corresponding to an IceModelVec2 object. The ncid
//! argument must refer to a dataset with dimensions t, x, y.
PetscErrorCode IceModelVec2::define_netcdf_variable(int ncid, nc_type nctype, int *varidp) {
  int stat, dimids[3], var_id;

  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "t", &dimids[0]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "x", &dimids[1]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "y", &dimids[2]); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_def_var(ncid, short_name, nctype, 3, dimids, &var_id);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  stat = MPI_Bcast(&var_id, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);

  *varidp = var_id;

  return 0;
}

PetscErrorCode IceModelVec2::get_array(PetscScalar** &a) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a = (PetscScalar**) array;
  return 0;
}


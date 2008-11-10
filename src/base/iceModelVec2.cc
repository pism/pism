// Copyright (C) 2008 Ed Bueler and Constantine Khroulev
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

// this file contains methods for derived classes IceModelVec2 and IceModelVec2Box

// methods for base class IceModelVec are in "iceModelVec.cc"

IceModelVec2::IceModelVec2() : IceModelVec() {}


PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_varname[], bool local) {

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2 with varname='%s' already allocated\n",my_varname);
  }
  PetscErrorCode ierr = create(my_grid, my_varname, local, DA_STENCIL_STAR); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode  IceModelVec2::createSameDA(IceModelVec2 imv2_source,
					   IceGrid &my_grid, const char my_varname[], bool local) {

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2 with varname='%s' already allocated\n",my_varname);
  }

  grid = &my_grid;

  da = imv2_source.da;
  IOwnDA = false;

  PetscErrorCode ierr;
  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  strcpy(varname,my_varname);
  return 0;
}

  

PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_varname[], bool local,
                                     DAStencilType my_sten) {

  grid = &my_grid;
  
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
  strcpy(varname,my_varname);
  return 0;
}

PetscErrorCode IceModelVec2::read(const char filename[], const unsigned int time) {
  PetscErrorCode ierr;
  ierr = read_from_netcdf(filename, time, 3, 1); CHKERRQ(ierr);
  return 0;
}

//! Defines a netcdf variable corresponding to an IceModelVec3 object. The ncid
// argument must refer to a dataset with dimensions t, x, y.
PetscErrorCode IceModelVec2::define_netcdf_variable(int ncid, nc_type nctype, int *varidp) {
  int stat, dimids[3], var_id;

  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "t", &dimids[0]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "x", &dimids[1]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "y", &dimids[2]); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_def_var(ncid, varname, nctype, 3, dimids, &var_id);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  stat = MPI_Bcast(&var_id, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);

  *varidp = var_id;

  return 0;
}

// Return value of ice scalar quantity stored in an IceModelVec2.
PetscScalar     IceModelVec2::getVal(const PetscInt i, const PetscInt j) {
  checkHaveArray();
  PetscScalar **arr = (PetscScalar**) array;
  return arr[i][j];
}


// Return values on planar star stencil of ice scalar quantity stored in an IceModelVec2.
PetscErrorCode   IceModelVec2::getPlaneStar(const PetscInt i, const PetscInt j, planeStar *star) {
  PetscErrorCode ierr;
  ierr = checkHaveArray();  CHKERRQ(ierr);
  // check ownership here?
  if (!localp) {
    SETERRQ1(1,"IceModelVec2 ERROR: IceModelVec2 with varname='%s' is GLOBAL and cannot do getPlaneStar()\n",
             varname);
  }
  
  PetscScalar **arr = (PetscScalar**) array;

  star->ij = arr[i][j];
  star->ip1 = arr[i+1][j];
  star->im1 = arr[i-1][j];
  star->jp1 = arr[i][j+1];
  star->jm1 = arr[i][j-1];
  return 0;
}


PetscScalar**   IceModelVec2::get_array() {
  checkHaveArray();
  return (PetscScalar**) array;
}



/********* IceModelVec2Box **********/

IceModelVec2Box::IceModelVec2Box() : IceModelVec2() {}


PetscErrorCode  IceModelVec2Box::create(IceGrid &my_grid, const char my_varname[], bool local) {

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2Box with varname='%s' already allocated\n",my_varname);
  }
  PetscErrorCode ierr = IceModelVec2::create(my_grid, my_varname, local, DA_STENCIL_BOX); CHKERRQ(ierr);
  return 0;
}


// Return values on planar BOX stencil of ice scalar quantity stored in an IceModelVec2Box.
PetscErrorCode   IceModelVec2Box::getPlaneBox(const PetscInt i, const PetscInt j, planeBox *box) {
  PetscErrorCode ierr;
  ierr = checkHaveArray();  CHKERRQ(ierr);
  // check ownership here?
  if (!localp) {
    SETERRQ1(1,"IceModelVec2Box ERROR: IceModelVec2Box with varname='%s' is GLOBAL and cannot do getPlaneBox()\n",
             varname);
  }

  PetscScalar **arr = (PetscScalar**) array;

  box->ij     = arr[i][j];
  box->ip1    = arr[i+1][j];
  box->im1    = arr[i-1][j];
  box->jp1    = arr[i][j+1];
  box->jm1    = arr[i][j-1];
  box->ip1jp1 = arr[i+1][j+1];
  box->im1jp1 = arr[i-1][j+1];
  box->ip1jm1 = arr[i+1][j-1];
  box->im1jm1 = arr[i-1][j-1];

  return 0;
}



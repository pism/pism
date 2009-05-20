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


PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_name[], bool local) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(1,"IceModelVec2 with name='%s' already allocated\n",my_name);
  }
  PetscErrorCode ierr = create(my_grid, my_name, local, DA_STENCIL_BOX); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_name[], bool local,
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

  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  strcpy(name,my_name);

  var1.init(my_name, my_grid, GRID_2D);

  return 0;
}

PetscErrorCode IceModelVec2::get_array(PetscScalar** &a) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a = (PetscScalar**) array;
  return 0;
}

//! Puts a local IceModelVec2 on processor 0.
/*!
 <ul>
 <li> onp0 and ctx should be created by calling VecScatterCreateToZero or be identical to one,
 <li> g2 is a preallocated temporary global vector,
 <li> g2natural is a preallocated temporary global vector with natural ordering.
 </ul>
*/
PetscErrorCode IceModelVec2::put_on_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);

  if (!localp)
    SETERRQ1(1, "Can't put a global IceModelVec '%s' on proc 0.", name);

  ierr =        DALocalToGlobal(da, v,  INSERT_VALUES, g2);        CHKERRQ(ierr);
  ierr = DAGlobalToNaturalBegin(grid->da2, g2, INSERT_VALUES, g2natural); CHKERRQ(ierr);
  ierr =   DAGlobalToNaturalEnd(grid->da2, g2, INSERT_VALUES, g2natural); CHKERRQ(ierr);

  ierr = VecScatterBegin(ctx, g2natural, onp0, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr =   VecScatterEnd(ctx, g2natural, onp0, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

  return 0;
}

//! Gets a local IceModelVec2 from processor 0.
/*!
 <ul>
 <li> onp0 and ctx should be created by calling VecScatterCreateToZero or be identical to one,
 <li> g2 is a preallocated temporary global vector,
 <li> g2natural is a preallocated temporary global vector with natural ordering.
 </ul>
*/
PetscErrorCode IceModelVec2::get_from_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);

  if (!localp)
    SETERRQ1(1, "Can't get a global IceModelVec '%s' from proc 0.", name);

  ierr = VecScatterBegin(ctx, onp0, g2natural, INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);
  ierr =   VecScatterEnd(ctx, onp0, g2natural, INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = DANaturalToGlobalBegin(grid->da2, g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr =   DANaturalToGlobalEnd(grid->da2, g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr =   DAGlobalToLocalBegin(da, g2,               INSERT_VALUES, v);  CHKERRQ(ierr);
  ierr =     DAGlobalToLocalEnd(da, g2,               INSERT_VALUES, v);  CHKERRQ(ierr);

  return 0;
}

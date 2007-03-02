// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <petscfix.h>
#include <petscda.h>
#include <petscbag.h>
#include "grid.hh"

// Note that the following choices may be overridden according to the input file.  
// That is, the
// data always has highest precedence.  The user should be warned when their
// command line options are overriden, but is not implemented at this time.  The
// netCDF file input respects the number of grid points, but not the physical
// extent.
static const PetscScalar HALFWIDTH_X = 1500.0e3;  // domain width is 2.0*HALFWIDTH_X
static const PetscScalar HALFWIDTH_Y = 1500.0e3;
static const PetscScalar HEIGHT_Z = 4000.0;
static const PetscScalar BEDROCK_HEIGHT_Z = 1000.0;
static const PetscInt GRIDPTS_X = 61;
static const PetscInt GRIDPTS_Y = 61;
static const PetscInt GRIDSPACES_Z = 31;


int initIceParam(MPI_Comm com, IceParam **param, PetscBag *bag) {
  PetscErrorCode ierr;
  PetscBag b;

  if (*bag != PETSC_NULL) {
    ierr = PetscPrintf(com, "initIceParam(): *bag = 0x%xp\n", param); CHKERRQ(ierr);
    return 0;
  }
  
  ierr = PetscBagCreate(com, sizeof(IceParam), &b); CHKERRQ(ierr);
  ierr = PetscBagGetData(b, (void **)param); CHKERRQ(ierr);
  *bag = b;
  IceParam &p = **param;

  ierr = PetscBagSetName(b, "IceParamBag", "Contains parameters controlling\n"
                         "preferences for simulation of ice sheets."); CHKERRQ(ierr);
  ierr = PetscBagRegisterString(b, &p.history, HISTORY_STRING_LENGTH, "\n", "history",
                                "History of commands used to generate this file.");
  CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.Lx, HALFWIDTH_X, "Lx",
                                "Half width of the ice model grid in x-direction (m).");
  CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.Ly, HALFWIDTH_Y, "Ly",
                                "Halfwidth of the ice model grid in y-direction (m).");
  CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.Lz, HEIGHT_Z, "Lz",
                                "Extent of the ice model grid in z-direction (m).");
  CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.Lbz, BEDROCK_HEIGHT_Z, "Lbz",
                                "Extent of the bedrock model grid in z-direction (m); DO NOT MODIFY; MOD OF Mbz OK.");
  CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(b, &p.Mx, GRIDPTS_X, "Mx",
                             "Number of grid points in x-direction."); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(b, &p.My, GRIDPTS_Y, "My",
                             "Number of grid points in y-direction."); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(b, &p.Mz, GRIDSPACES_Z, "Mz",
                             "Number of ice grid points in z-direction."); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(b, &p.Mbz, 0, "Mbz",
                             "Number of bedrock grid points in z-direction."); CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.dx, 2.0 * HALFWIDTH_X / (GRIDPTS_X - 1), "dx",
                                "Grid spacing in x-direction (m); DO NOT MODIFY; MOD OF Lx,Mx OK."); CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.dy, 2.0 * HALFWIDTH_Y / (GRIDPTS_Y - 1), "dy",
                                "Grid spacing in y-direction (m); DO NOT MODIFY; MOD OF Ly,My OK."); CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.dz, HEIGHT_Z / (GRIDSPACES_Z - 1), "dz",
                                "Grid spacing in z-direction (m); DO NOT MODIFY; MOD OF Lz,Mz OK."); CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.year, 0, "year",
                                "Current time in years."); CHKERRQ(ierr);

  return 0;
}


IceGrid::IceGrid(MPI_Comm c,
                 PetscMPIInt r,
                 PetscMPIInt s):
  com(c), rank(r), size(s), p(PETSC_NULL), bag(PETSC_NULL), createDA_done(PETSC_FALSE) { }


IceGrid::IceGrid(MPI_Comm c,
                 PetscMPIInt r,
                 PetscMPIInt s,
                 IceParam *prm,
                 PetscBag b):
  com(c), rank(r), size(s), p(prm), bag(b), createDA_done(PETSC_FALSE) { }


IceGrid::~IceGrid() {
  if (destroyDA() != 0) {
    PetscEnd();
  }
  if (PetscBagDestroy(bag) != 0) {
    PetscEnd();
  }
}


PetscErrorCode IceGrid::createDA() {
  PetscErrorCode ierr;
  PetscInt    M, N, m, n;
  DALocalInfo info;

  if (createDA_done == PETSC_TRUE) {
    ierr = destroyDA(); CHKERRQ(ierr);
  }
  
  ierr = DACreate2d(com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    p->My, p->Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 1,
                    PETSC_NULL, PETSC_NULL, &da2); CHKERRQ(ierr);
  ierr = DAGetInfo(da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate3d(com, DA_YZPERIODIC, DA_STENCIL_STAR, p->Mz, N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da3); CHKERRQ(ierr);
  ierr = DACreate3d(com, DA_YZPERIODIC, DA_STENCIL_STAR, (p->Mbz > 0) ? p->Mbz : 1,
                    N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da3b); CHKERRQ(ierr);

  ierr = DAGetLocalInfo(da2, &info); CHKERRQ(ierr);
  xs = info.ys; xm = info.ym;
  ys = info.xs; ym = info.xm;

  ierr = setCoordinatesDA(); CHKERRQ(ierr);

  createDA_done = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceGrid::setCoordinatesDA() {
  PetscErrorCode ierr;

  ierr = DASetUniformCoordinates(da2, -p->Ly, p->Ly, -p->Lx, p->Lx,
                                 PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DASetUniformCoordinates(da3, 0, p->Lz, -p->Ly, p->Ly, -p->Lx, p->Lx); CHKERRQ(ierr);
  ierr = DASetUniformCoordinates(da3b, ((p->Mbz == 0) ? -1.0 : -p->Lbz), 0.0,
                                 -p->Ly, p->Ly, -p->Lx, p->Lx); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode IceGrid::destroyDA() {
  PetscErrorCode ierr;

  ierr = DADestroy(da2); CHKERRQ(ierr);
  ierr = DADestroy(da3); CHKERRQ(ierr);
  ierr = DADestroy(da3b); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceGrid::rescale(PetscScalar lx, PetscScalar ly, PetscScalar lz) {
  PetscErrorCode ierr;
  
  p->Lx = lx; p->Ly = ly; p->Lz = lz;
  p->dx = 2.0 * p->Lx / (p->Mx - 1);
  p->dy = 2.0 * p->Ly / (p->My - 1);
  p->dz = p->Lz / (p->Mz - 1);
  p->Lbz = p->dz * p->Mbz;

  ierr = setCoordinatesDA(); CHKERRQ(ierr);
  
  return 0;
}

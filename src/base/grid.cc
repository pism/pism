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
#include "grid.hh"

// Note that choices herein may be overridden according to the input file.  
// That is, the data always has highest precedence.  The user should be warned when their
// command line options are overriden, but is not implemented at this time.  The
// NetCDF file input respects the number of grid points, but not the physical extent.

static const PetscScalar HALFWIDTH_X = 1500.0e3;  // domain width is 2.0*HALFWIDTH_X
static const PetscScalar HALFWIDTH_Y = 1500.0e3;
static const PetscScalar HEIGHT_Z = 4000.0;
static const PetscInt GRIDPTS_X = 61;
static const PetscInt GRIDPTS_Y = 61;
static const PetscInt GRIDPTS_Z = 31;


int initIceParam(MPI_Comm com, IceParam **param) {
  IceParam* p = *param;

  p->Lx = HALFWIDTH_X;
  p->Ly = HALFWIDTH_Y;
  p->Lz = HEIGHT_Z;
  p->Mx = GRIDPTS_X;
  p->My = GRIDPTS_Y;
  p->Mz = GRIDPTS_Z;
  p->dx = 2.0 * HALFWIDTH_X / (GRIDPTS_X - 1);
  p->dy = 2.0 * HALFWIDTH_Y / (GRIDPTS_Y - 1);
  p->dz = HEIGHT_Z / (GRIDPTS_Z - 1);
  p->Mbz = 1;
  p->Lbz = p->Mbz * p->dz;
  p->year = 0.0;

/*
  PetscErrorCode ierr;
  ierr = PetscBagRegisterString(b, &p.history, HISTORY_STRING_LENGTH, "\n", "history",
                                "History of commands used to generate this file."); CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.Lx, HALFWIDTH_X, "Lx",
                                "Half width of the ice model grid in x-direction (m).");  CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.Ly, HALFWIDTH_Y, "Ly",
                                "Halfwidth of the ice model grid in y-direction (m).");  CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.Lz, HEIGHT_Z, "Lz",
                                "Extent of the ice model grid in z-direction (m).");  CHKERRQ(ierr);
  const PetscScalar init_dz = HEIGHT_Z / (GRIDPTS_Z - 1);
  ierr = PetscBagRegisterScalar(b, &p.Lbz, GRIDPTS_Z * init_dz, "Lbz",
           "Extent of the bedrock model grid in z-direction (m); DO NOT MODIFY; MOD OF Mbz OK.");  CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(b, &p.Mx, GRIDPTS_X, "Mx",
                             "Number of grid points in x-direction."); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(b, &p.My, GRIDPTS_Y, "My",
                             "Number of grid points in y-direction."); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(b, &p.Mz, GRIDPTS_Z, "Mz",
                             "Number of ice grid points in z-direction."); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(b, &p.Mbz, 1, "Mbz",
                             "Number of bedrock grid points in z-direction."); CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.dx, 2.0 * HALFWIDTH_X / (GRIDPTS_X - 1), "dx",
                                "Grid spacing in x-direction (m); DO NOT MODIFY; MOD OF Lx,Mx OK."); CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.dy, 2.0 * HALFWIDTH_Y / (GRIDPTS_Y - 1), "dy",
                                "Grid spacing in y-direction (m); DO NOT MODIFY; MOD OF Ly,My OK."); CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.dz, init_dz, "dz",
                                "Grid spacing in z-direction (m); DO NOT MODIFY; MOD OF Lz,Mz OK."); CHKERRQ(ierr);
  ierr = PetscBagRegisterScalar(b, &p.year, 0, "year",
                                "Current time in years."); CHKERRQ(ierr);
*/
  return 0;
}


IceGrid::IceGrid(MPI_Comm c,
                 PetscMPIInt r,
                 PetscMPIInt s):
  com(c), rank(r), size(s), createDA_done(PETSC_FALSE) { 
  p = new IceParam;
}


IceGrid::IceGrid(MPI_Comm c,
                 PetscMPIInt r,
                 PetscMPIInt s,
                 IceParam *prm):
  com(c), rank(r), size(s), p(prm), createDA_done(PETSC_FALSE) { }


IceGrid::~IceGrid() {
  if (destroyDA() != 0) {
    PetscPrintf(com, "IceGrid destructor: invalid destroyDA() return; ENDING\n");
    PetscEnd();
  }
  delete p;
/*
  if (PetscBagDestroy(bag) != 0) {
    PetscPrintf(com, "IceGrid destructor: invalid PetscBagDestroy() return; ENDING\n");
    PetscEnd();
  }
*/
}


/* order of allocation and settings (for rev 113 and later, sans PetscBag):
	1. initIceParam() which sets defaults for Mx,My,Mz,Mbz,Lx,Ly,Lz,Lbz,dx,dy,dz,year
	[2]. derivedClass:setFromOptions() to get options special to derived class
	3. setFromOptions() to get all options *including* Mx,My,Mz,Mbz
	[4]. initFromFile_netCDF() which reads Mx,My,Mz,Mbz from file and overwrites previous; if 
	   this represents a change the user should be warned
	5. createDA() uses only Mx,MyMz,Mbz
	6. createVecs() uses DA to create/allocate Vecs
	[7]. derivedClass:createVecs() to create/allocate Vecs special to derived class
	8. afterInitHook() which changes Lx,Ly,Lz if set by user

	Note driver programs call only setFromOptions() and initFromOptions() (for IceModel or derived class).

	initIceParam() is called by the constructor for IceModel.

	IceModel::setFromOptions() should be called at the end of derivedClass:setFromOptions().

	Note 4,5,6 are called from initFromFile() in IceModel *or* initFromOptions() in some derived classes
	(e.g. IceCompModel) but not both.

	Note step 4 is skipped when bootstrapping (-bif and bootstrapFromFile_netCDF() or in
	those derived classes which can start with no input files, e.g. IceCompModel).  That is,
	4 is only done when starting from a saved model state.
*/

PetscErrorCode IceGrid::createDA() {
  PetscErrorCode ierr;
  PetscInt    M, N, m, n;
  DALocalInfo info;

  if (createDA_done == PETSC_TRUE) {
    ierr = destroyDA(); CHKERRQ(ierr);
  }

#if (MARGIN_TRICK_TWO)
  ierr = PetscPrintf(com,"  MARGIN_TRICK_TWO createDA(): setting stencil width to 2 ...\n"); CHKERRQ(ierr);
  PetscInt stencilwidth = 2;
  ierr = DACreate2d(com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    p->My, p->Mx, PETSC_DECIDE, PETSC_DECIDE, 1, stencilwidth,
                    PETSC_NULL, PETSC_NULL, &da2); CHKERRQ(ierr);
#else
  ierr = DACreate2d(com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    p->My, p->Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 1,
                    PETSC_NULL, PETSC_NULL, &da2); CHKERRQ(ierr);
#endif
  ierr = DAGetInfo(da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate3d(com, DA_YZPERIODIC, DA_STENCIL_STAR, p->Mz, N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da3); CHKERRQ(ierr);
  ierr = DACreate3d(com, DA_YZPERIODIC, DA_STENCIL_STAR, p->Mbz, N, M, 1, n, m, 1, 1,
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
  ierr = DASetUniformCoordinates(da3, 0, p->Lz,
                                 -p->Ly, p->Ly, -p->Lx, p->Lx); CHKERRQ(ierr);
  if (p->Mbz > 1) {
    ierr = DASetUniformCoordinates(da3b, -p->Lbz, 0.0,
                                   -p->Ly, p->Ly, -p->Lx, p->Lx); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceGrid::destroyDA() {
  PetscErrorCode ierr;

  ierr = DADestroy(da2); CHKERRQ(ierr);
  ierr = DADestroy(da3); CHKERRQ(ierr);
  ierr = DADestroy(da3b); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceGrid::rescale(const PetscScalar lx, const PetscScalar ly, 
                                const PetscScalar lz) {
  PetscErrorCode ierr;
  
  ierr = rescale(lx,ly,lz,PETSC_FALSE); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceGrid::rescale(const PetscScalar lx, const PetscScalar ly, 
                                const PetscScalar lz, const PetscTruth truelyPeriodic) {
  PetscErrorCode ierr;

  if (lx<=0) {
    SETERRQ(1, "rescale: error, lx must be positive\n");
  }
  if (ly<=0) {
    SETERRQ(1, "rescale: error, ly must be positive\n");
  }
  if (lz<=0) {
    SETERRQ(1, "rescale: error, lz must be positive\n");
  }

  p->Lx = lx; p->Ly = ly; p->Lz = lz;
  if (truelyPeriodic == PETSC_TRUE) {
    p->dx = 2.0 * p->Lx / (p->Mx);
    p->dy = 2.0 * p->Ly / (p->My);
  } else {
    p->dx = 2.0 * p->Lx / (p->Mx - 1);
    p->dy = 2.0 * p->Ly / (p->My - 1);
  }

  p->dz = p->Lz / (p->Mz - 1);
  p->Lbz = p->dz * (p->Mbz - 1);

  ierr = setCoordinatesDA(); CHKERRQ(ierr);
  return 0;
}

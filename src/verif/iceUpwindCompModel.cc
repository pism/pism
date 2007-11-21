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


#include <petscvec.h>
#include "iceCompModel.hh"
#include "iceUpwindCompModel.hh"


IceUpwindCompModel::IceUpwindCompModel(IceGrid &g, ThermoGlenArrIce &i, const char mytest)
  : IceCompModel(g, i, mytest) {
}


PetscErrorCode IceUpwindCompModel::initFromOptions() {
  PetscErrorCode ierr;

  ierr = IceCompModel::initFromOptions(); CHKERRQ(ierr);

  ierr = VecSet(vMask,(PetscScalar) MASK_DRAGGING); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,
           "IceUpwindCompModel:  mask set to all dragging so mass continuity is upwind method\n");
           CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceUpwindCompModel::velocity(bool updateSIAVelocityAtDepth) {
  PetscErrorCode ierr;
  
  ierr = IceCompModel::velocity(updateSIAVelocityAtDepth);  CHKERRQ(ierr);

  PetscScalar **mask, **ubar, **vbar;
  PetscScalar locCFLmaxdt2D = maxdt;
  
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (intMask(mask[i][j]) != MASK_SHEET) {
        PetscScalar denom = PetscAbs(ubar[i][j])/grid.p->dx + PetscAbs(vbar[i][j])/grid.p->dy;
        denom += (0.01/secpera)/(grid.p->dx + grid.p->dy);  // make sure it's pos.
        locCFLmaxdt2D = PetscMin(locCFLmaxdt2D,1.0/denom);
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);

  ierr = PetscGlobalMin(&locCFLmaxdt2D, &CFLmaxdt2D, grid.com); CHKERRQ(ierr);
  //ierr = verbPrintf(1,grid.com,"\nIceUpwindCompModel::velocity():  CFLmaxdt2D = %10.4f",
  //                  CFLmaxdt2D/secpera); CHKERRQ(ierr);

  maxdt_temporary = CFLmaxdt2D;
  return 0;
}




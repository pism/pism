// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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
#include "iceEISplModel.hh"


IceEISplModel::IceEISplModel(IceGrid &g, IceType &i)
  : IceEISModel(g,i) {  // do nothing; note derived classes must have constructors
  expername = 'I';
}


PetscErrorCode IceEISplModel::initFromOptions() {
  PetscErrorCode      ierr;

  useSSAVelocity = PETSC_TRUE;
  doSuperpose = PETSC_TRUE;
  pureSuperpose = PETSC_FALSE;
  doPlasticTill = PETSC_TRUE;

  // these are different from EISMINT I conventions
  updateHmelt = PETSC_TRUE;

  ierr = IceEISModel::initFromOptions(); CHKERRQ(ierr);  

#define DEFAULT_TILL_PHI_LAKE     20.0  // no lake by default
#define DEFAULT_TILL_PHI_STRONG   20.0
#define DEFAULT_TILL_PHI_WEAK      5.0
#define DEFAULT_TILL_PHI_OCEAN     0.0
  PetscScalar default_phi_list[4] = { DEFAULT_TILL_PHI_LAKE, DEFAULT_TILL_PHI_STRONG,
                                      DEFAULT_TILL_PHI_WEAK, DEFAULT_TILL_PHI_OCEAN };
  PetscInt    phi_list_length = 4;
  PetscTruth  phi_listSet;
  
  phi_list = (PetscScalar*) default_phi_list;
  ierr = PetscOptionsGetRealArray(PETSC_NULL, "-till_phi", phi_list, &phi_list_length,
            &phi_listSet); CHKERRQ(ierr);
  if (phi_listSet == PETSC_TRUE) {
    if (phi_list_length != 4) {
      SETERRQ(1,"option -till_phi must be followed by comma-separated list (no spaces!) of\n"
              "   of four real values\n");
    }
    ierr = verbPrintf(2,grid.com, "-till_phi option read; "); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com, "default phi "); CHKERRQ(ierr);
  }
  ierr = verbPrintf(2,grid.com, 
           "values are LAKE = %6.2f, STRONG = %6.2f, WEAK = %6.2f, OCEAN = %6.2f\n",
           phi_list[0],phi_list[1],phi_list[2],phi_list[3]); CHKERRQ(ierr);
  
  ierr = setTillProperties(); CHKERRQ(ierr);
  ierr = resetAccum(); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com, 
           "running plastic till SSA modification of EISMINT II experiment I ...\n");
           CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceEISplModel::resetAccum() {
  PetscErrorCode    ierr;
  PetscScalar       **accum;

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  PetscScalar cx = grid.p->Lx, cy = grid.p->Ly;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // r is distance from center of grid
      const PetscScalar r = sqrt( PetscSqr(-cx + grid.p->dx*i) + PetscSqr(-cy + grid.p->dy*j) );
      // set accumulation outside of sheet to very negative to eliminate "shelf"
      if (r > 650.0e3)   accum[i][j] = -10.0 / secpera;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  return 0;
}



PetscErrorCode IceEISplModel::setTillProperties() {
  PetscErrorCode  ierr;
  
  const PetscScalar    L = 750.0e3;  // half-width of computational domain
  const PetscScalar    dx = grid.p->dx, dy = grid.p->dy;
  const PetscScalar    dx61 = (2*L) / 60; // = 25.0e3
  
  // fill in map of phi = friction angle for till
  useConstantTillPhi = PETSC_FALSE;
  PetscScalar  **tillphi;
  ierr = DAVecGetArray(grid.da2, vtillphi, &tillphi); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar nsd = i * dx, ewd = j *dy;  // north-south and east-west distances
      if (    (nsd >= (29 - 1) * dx61) && (nsd <= (33 - 1) * dx61) 
           && (ewd >= (31 - 1) * dx61)                             ) { // in strip
        if (ewd >= (57 - 1) * dx61) {
          tillphi[i][j] = phi_list[3];  // OCEAN
        } else {
          tillphi[i][j] = phi_list[2];  // WEAK
        }
      } else if (    (nsd >= (33 - 1) * dx61) && (nsd <= (37 - 1) * dx61) 
                  && (ewd >= (35 - 1) * dx61) && (ewd <= (39 - 1) * dx61) )     { // in lake
        tillphi[i][j] = phi_list[0];  // LAKE
      } else {
        tillphi[i][j] = phi_list[1];  // STRONG
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vtillphi, &tillphi); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
         "map of phi = (till friction angle) stored by IceEISplModel\n"); CHKERRQ(ierr);

  return 0;
}


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
#include <petsc.h>
#include "iceEISplModel.hh"

const PetscScalar IceEISplModel::DEFAULT_STREAM_WIDTH = 100.0e3; // 100 km

// default case has no lake and has the upstream and downstream part of the ice
//   stream have the same till strength
const PetscInt    IceEISplModel::phi_list_length = 5;
const PetscScalar IceEISplModel::DEFAULT_TILL_PHI_LAKE        = 20.0;
const PetscScalar IceEISplModel::DEFAULT_TILL_PHI_STRONG      = 20.0;
const PetscScalar IceEISplModel::DEFAULT_TILL_PHI_UPSTREAM    =  5.0;
const PetscScalar IceEISplModel::DEFAULT_TILL_PHI_DOWNSTREAM  =  5.0;
const PetscScalar IceEISplModel::DEFAULT_TILL_PHI_OCEAN       =  0.0;


IceEISplModel::IceEISplModel(IceGrid &g, IceType &i)
  : IceEISModel(g,i) {  // do nothing; note derived classes must have constructors

  expername = 'I';
}


//! Set defaults relative to EISMINT II experiment I and read options -stream_width and -till_phi.
PetscErrorCode IceEISplModel::initFromOptions() {
  PetscErrorCode      ierr;

  useSSAVelocity = PETSC_TRUE;
  doSuperpose = PETSC_TRUE;
  pureSuperpose = PETSC_FALSE;
  doPlasticTill = PETSC_TRUE;

  // these are different from EISMINT I conventions
  updateHmelt = PETSC_TRUE;
  includeBMRinContinuity = PETSC_TRUE;

  ierr = IceEISModel::initFromOptions(); CHKERRQ(ierr);  
  
  stream_width = DEFAULT_STREAM_WIDTH;
  PetscTruth  stream_widthSet;  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-stream_width", &stream_width,
            &stream_widthSet); CHKERRQ(ierr);
  if (stream_widthSet == PETSC_TRUE) {
    stream_width = stream_width * 1000.0;  // user enters stream width in km
  }

  PetscScalar default_phi_list[phi_list_length] = { 
             DEFAULT_TILL_PHI_LAKE, 
             DEFAULT_TILL_PHI_STRONG,
             DEFAULT_TILL_PHI_UPSTREAM, 
             DEFAULT_TILL_PHI_DOWNSTREAM, 
             DEFAULT_TILL_PHI_OCEAN };
  phi_list = (PetscScalar*) default_phi_list;
  PetscTruth  phi_listSet;
  PetscInt    input_phi_list_length = phi_list_length;
  ierr = PetscOptionsGetRealArray(PETSC_NULL, "-till_phi", phi_list, &input_phi_list_length,
            &phi_listSet); CHKERRQ(ierr);
  if (phi_listSet == PETSC_TRUE) {
    if (input_phi_list_length != phi_list_length) {
      SETERRQ2(1,"option -till_phi must be followed by comma-separated list (no spaces!) of\n"
              "   of exactly %d real values; %d were given\n", phi_list_length, input_phi_list_length);
    }
    ierr = verbPrintf(2,grid.com, "[-till_phi option read; "); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com, "[default; "); CHKERRQ(ierr);
  }
  ierr = verbPrintf(2,grid.com, 
           "till friction angles are LAKE = %5.2f, STRONG = %5.2f,\n"
           "   UPSTREAM = %5.2f, DOWNSTREAM = %5.2f, OCEAN = %5.2f]\n",
           phi_list[0],phi_list[1],phi_list[2],phi_list[3],phi_list[4]); CHKERRQ(ierr);
  
  ierr = setTillProperties(); CHKERRQ(ierr);
  ierr = resetAccum(); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com, 
           "running plastic till SSA modification of EISMINT II experiment I ...\n");
           CHKERRQ(ierr);
  return 0;
}


//! Set accumulation outside of sheet very negative to stop flow to edge of computational domain.
PetscErrorCode IceEISplModel::resetAccum() {
  PetscErrorCode    ierr;
  PetscScalar       **accum;

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  PetscScalar cx = grid.p->Lx, cy = grid.p->Ly;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // r is distance from center of grid
      const PetscScalar r = sqrt( PetscSqr(-cx + grid.p->dx*i) + PetscSqr(-cy + grid.p->dy*j) );
      if (r > 650.0e3)   accum[i][j] = -10.0 / secpera;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  return 0;
}


//! Set the till friction angle to five values: LAKE, STRONG, UPSTREAM, DOWNSTREAM, OCEAN.
PetscErrorCode IceEISplModel::setTillProperties() {
  PetscErrorCode  ierr;
  
  const PetscScalar    L = 750.0e3;  // half-width of computational domain
  const PetscScalar    dx = grid.p->dx, dy = grid.p->dy;
  const PetscScalar    dx61 = (2.0 * L) / 60.0; // = 25.0e3
  
  // distances are based on 61 x 61 grid because came from EISMINT II experiment I
  const PetscScalar  streamBegin  = (31 - 1) * dx61,
                     streamTop    = (31 + (int) ((stream_width + 1.0)/(2.0 * dx61)) - 1) * dx61,
                     streamBottom = (31 - (int) ((stream_width + 1.0)/(2.0 * dx61)) - 1) * dx61,
                     lakeWidth    = 4 * dx61,
                     lakeLeft     = (35 - 1) * dx61,
                     lakeRight    = lakeLeft + lakeWidth,
                     upDownDivide = (49 - 1) * dx61,
                     groundLine   = (57 - 1) * dx61;
  ierr = verbPrintf(2,grid.com,
         "[streamBegin = %9.2f, streamTop = %9.2f, streamBottom = %9.2f,\n"
         "    lakeWidth = %9.2f, lakeLeft = %9.2f,  lakeRight = %9.2f,\n"
         "    upDownDivide = %9.2f,  groundLine = %9.2f]\n",
         streamBegin,streamTop,streamBottom,lakeWidth,lakeLeft,lakeRight,
         upDownDivide,groundLine); CHKERRQ(ierr);
                       
  // fill in map of phi = friction angle for till
  useConstantTillPhi = PETSC_FALSE;
  PetscScalar  **tillphi;
  ierr = DAVecGetArray(grid.da2, vtillphi, &tillphi); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar nsd = i * dx, ewd = j * dy;  // north-south and east-west distances
      if (    (nsd >= streamBottom) && (nsd <= streamTop) && (ewd >= streamBegin) ) { // in strip
        if (ewd >= groundLine) {
          tillphi[i][j] = phi_list[4];  // OCEAN
        } else if (ewd >= upDownDivide) {
          tillphi[i][j] = phi_list[3];  // DOWNSTREAM
        } else {
          tillphi[i][j] = phi_list[2];  // UPSTREAM
        }
      } else if (    (nsd >= streamTop) && (nsd <= streamTop + lakeWidth) 
                  && (ewd >= lakeLeft) && (ewd <= lakeRight) )     { // in lake
        tillphi[i][j] = phi_list[0];    // LAKE
      } else {
        tillphi[i][j] = phi_list[1];    // STRONG
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vtillphi, &tillphi); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
         "[map of phi = (till friction angle) stored by IceEISplModel]\n"); CHKERRQ(ierr);

  return 0;
}



PetscErrorCode IceEISplModel::summaryPrintLine(
    const PetscTruth printPrototype, const PetscTruth tempAndAge,
    const PetscScalar year, const PetscScalar dt, 
    const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
    const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0) {

  PetscErrorCode ierr;

  PetscScalar     **H, **ubar, **vbar, **ub, **vb;
  // these are in MKS; sans "g" are local to the processore; with "g" are global
  PetscScalar     maxcbar = 0.0, maxcb = 0.0, avcb = 0.0, Nhaveice = 0.0;
  PetscScalar     gmaxcbar, gmaxcb, gavcb, gNhaveice;
  
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        Nhaveice += 1.0;
        const PetscScalar cbar = sqrt( PetscSqr(ubar[i][j]) + PetscSqr(vbar[i][j]) ),
                          cb   = sqrt( PetscSqr(ub[i][j])   + PetscSqr(vb[i][j])   );
        if (cbar > maxcbar)  maxcbar = cbar;
        if (cb > maxcb)      maxcb = cb;
        avcb += cb;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxcbar, &gmaxcbar, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxcb,   &gmaxcb,   grid.com); CHKERRQ(ierr);
  
  ierr = PetscGlobalSum(&Nhaveice,&gNhaveice,grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avcb,    &gavcb,    grid.com); CHKERRQ(ierr);
  if (gNhaveice > 0) { // area_kmsquare must already have been globalized
    gavcb = gavcb / gNhaveice;
  } else {  // some kind of degenerate case
    gavcb = 0.0;
  }

  if (printPrototype == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,
      "P       YEAR (+     STEP )      VOL    AREA   MELTF   MAXCBAR     MAXCB     AVCB\n");
    ierr = verbPrintf(2,grid.com,
      "U      years       years  10^6_km^3 10^6_km^2 (none)      m/a       m/a      m/a\n");
  } else {
    ierr = verbPrintf(2,grid.com, "S %10.3f (+ %8.5f ) %8.5f %7.4f %7.4f %9.2f %9.2f %8.2f\n",
                      year, dt/secpera, 
                      volume_kmcube/1.0e6, area_kmsquare/1.0e6, meltfrac,
                      gmaxcbar*secpera,gmaxcb*secpera,gavcb*secpera); CHKERRQ(ierr);
  }
  return 0;
}

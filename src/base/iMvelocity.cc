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

#include <cmath>
#include <petscda.h>
#include "iceModel.hh"

//! Manage the computation of velocity and do the necessary communication.
/*!
This procedure calls the important routines velocitySIAStaggered(), velocitySSA(),
and vertVelocityFromIncompressibility() according to various flags. 
 */
PetscErrorCode IceModel::velocity(bool updateVelocityAtDepth) {
  PetscErrorCode ierr;
  static PetscTruth firstTime = PETSC_TRUE;

  // do SIA
  if (computeSIAVelocities == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com, " SIA "); CHKERRQ(ierr);
    ierr = surfaceGradientSIA(); CHKERRQ(ierr); // comm may happen here ...
    // surface gradient temporarily stored in vWork2d[0 1 2 3] 
    ierr = velocitySIAStaggered(); CHKERRQ(ierr);  // fills vWork3d[0 1 2 3]

    // communicate vuvbar[01] for boundary conditions for SSA and vertAveragedVelocityToRegular()
    // and velocities2DSIAToRegular()
    ierr = DALocalToLocalBegin(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);

    // compute (and initialize values in) ub,vb and Rb; zero everything where floating
    ierr = basalSIA(); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vub, INSERT_VALUES, vub); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vub, INSERT_VALUES, vub); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vvb, INSERT_VALUES, vvb); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vvb, INSERT_VALUES, vvb); CHKERRQ(ierr); 
    // no need to communicate vRb since not differenced in horizontal

    // now put staggered grid value of vertically-averaged horizontal velocity on regular grid
    // (this makes ubar and vbar be just the result of SIA; note that for SSA
    // we used *saved* ubar,vbar as first guess for iteration at SSA points)
    // also put staggered value of basal velocity onto regular grid
    ierr = velocities2DSIAToRegular(); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr); 
  
    if (updateVelocityAtDepth) {  
      // communicate h_x[o], h_y[o] on staggered for horizontalVelocitySIARegular()
      for (PetscInt k=0; k<4; ++k) { 
        ierr = DALocalToLocalBegin(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
        ierr = DALocalToLocalEnd(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
      }
      // communicate I on staggered for horizontalVelocitySIARegular()
      //   and also communicate Sigma on staggered for SigmaSIAToRegular()
/*
      for (PetscInt k=0; k<4; ++k) {
        ierr = DALocalToLocalBegin(grid.da3, vWork3d[k], INSERT_VALUES, vWork3d[k]); CHKERRQ(ierr);
        ierr = DALocalToLocalEnd(grid.da3, vWork3d[k], INSERT_VALUES, vWork3d[k]); CHKERRQ(ierr);
      }
*/
      ierr = Sigmastag3[0].beginGhostComm(); CHKERRQ(ierr);
      ierr = Sigmastag3[1].beginGhostComm(); CHKERRQ(ierr);
      ierr = Istag3[0].beginGhostComm(); CHKERRQ(ierr);
      ierr = Istag3[1].beginGhostComm(); CHKERRQ(ierr);
      ierr = Sigmastag3[0].endGhostComm(); CHKERRQ(ierr);
      ierr = Sigmastag3[1].endGhostComm(); CHKERRQ(ierr);
      ierr = Istag3[0].endGhostComm(); CHKERRQ(ierr);
      ierr = Istag3[1].endGhostComm(); CHKERRQ(ierr);

      ierr = SigmaSIAToRegular(); CHKERRQ(ierr);
      ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);
// done with vWork2d[0 1 2 3] and vWork3d[0 1 2 3]
    }
  } else { // if computeSIAVelocities
    ierr = verbPrintf(2,grid.com, "     "); CHKERRQ(ierr);
  }
  
  // do SSA
  if ((firstTime == PETSC_TRUE) && (useSSAVelocity == PETSC_TRUE)) {
    ierr = initSSA(); CHKERRQ(ierr);
  }
  static PetscScalar lastSSAUpdateYear = grid.p->year;  // only set when first called
  if (useSSAVelocity) { // communication happens within SSA
    if ((firstTime == PETSC_TRUE) || (grid.p->year - lastSSAUpdateYear >= ssaIntervalYears)) {
      PetscInt numSSAiter;
      ierr = setupGeometryForSSA(DEFAULT_MINH_SSA); CHKERRQ(ierr);
      ierr = velocitySSA(&numSSAiter); CHKERRQ(ierr); // comm here ...
      ierr = cleanupGeometryAfterSSA(DEFAULT_MINH_SSA); CHKERRQ(ierr);
      lastSSAUpdateYear = grid.p->year;
      ierr = verbPrintf(2,grid.com, "SSA%3d ", numSSAiter); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "       "); CHKERRQ(ierr);
    }
    // even if velocitySSA() did not run, we still need to use stored SSA velocities 
    // to get 3D velocity field, basal velocities, basal frictional heating, 
    // and strain dissipation heating
    ierr = broadcastSSAVelocity(updateVelocityAtDepth); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com, "       "); CHKERRQ(ierr);
  }

  // now communicate modified velocity fields
  ierr = DALocalToLocalBegin(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vub, INSERT_VALUES, vub); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vub, INSERT_VALUES, vub); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vvb, INSERT_VALUES, vvb); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vvb, INSERT_VALUES, vvb); CHKERRQ(ierr);

  if (updateVelocityAtDepth || useSSAVelocity) {  // in latter case u,v are modified by broadcastSSAVelocity()
    ierr = u3.beginGhostComm(); CHKERRQ(ierr);
    ierr = v3.beginGhostComm(); CHKERRQ(ierr);
    ierr = u3.endGhostComm(); CHKERRQ(ierr);
    ierr = v3.endGhostComm(); CHKERRQ(ierr);
  }

  if (useSSAVelocity) {
    ierr = correctSigma(); CHKERRQ(ierr);  // note correctSigma() differences ub,vb in horizontal
    ierr = correctBasalFrictionalHeating(); CHKERRQ(ierr);
  }

  // finally update w
  if (updateVelocityAtDepth) {
    ierr = vertVelocityFromIncompressibility(); CHKERRQ(ierr);
    // no communication needed for w, which is only differenced in the column
  }
  
/* REMOVED TO AVOID STENCIL_BOX COMMUNICATION FOR 3D Vecs: 
  // smoothing Sigma is NOT RECOMMENDED, but do it here
  if (noSpokesLevel > 0) {
    ierr = smoothSigma(); CHKERRQ(ierr); // comm here
  }
*/

  // communication here for global max; sets CFLmaxdt2D
  ierr = computeMax2DSlidingSpeed(); CHKERRQ(ierr);   

  if ((useSSAVelocity == PETSC_TRUE) || (updateVelocityAtDepth == PETSC_TRUE)) {
    // communication here for global max; sets CFLmaxdt
    ierr = computeMax3DVelocities(); CHKERRQ(ierr); 
  }
  
  firstTime = PETSC_FALSE;
  return 0;
}


//! Compute vertical velocity using basal conditions and incompressibility of the ice.
/*! 
The original statement of incompressibility is
   \f[ \nabla\cdot\mathbf{U} + \frac{\partial w}{\partial z} = 0. \f]
This is immediately equivalent to the integral
   \f[ w(x,y,z,t) = - \int_{b(x,y,t)}^{z} \nabla\cdot\mathbf{U}\,d\zeta
                           + w_b(x,y,t). \f]

The basal kinematic equation is
   \f[ w_b = \frac{\partial b}{\partial t} + \mathbf{U}_b \cdot \nabla b - S \f]
where \f$S\f$ is the basal melt rate.  (The inclusion of the basal melt rate is optional.)  
This equation determines the vertical velocity of the ice at the base (\f$w_b\f$), 
which is needed in the incompressibility integral.

@cond VERTCHANGE
Note we do a change of vertical coordinate \f$\tilde z = z - b(x,y,t)\f$.  Thus all 
derivatives, with respect to any variable \f$x,y,z,t\f$, change accordingly 
(see elsewhere).  The revised form of the incompressibility integral is
   \f[ w(x,y,\tilde z,t) = - \int_0^{\tilde z} \left(
          \nabla\cdot\mathbf{U} 
          - \nabla b \cdot \left(\frac{\partial u}{\partial \tilde z},
                                 \frac{\partial v}{\partial \tilde z}\right)
                                                     \right)\,d\zeta
                           + w_b(x,y,t). \f]
@endcond

The vertical integral is computed by the trapezoid rule.
 */
PetscErrorCode IceModel::vertVelocityFromIncompressibility() {
  
  const PetscScalar   dx = grid.p->dx, dy = grid.p->dy, dz = grid.p->dz;
  const PetscInt      Mz = grid.p->Mz;
  PetscErrorCode  ierr;
  PetscScalar **ub, **vb, **basalMeltRate, **b, **dbdt;

  PetscScalar *izz, *u, *v, *w;
  izz = new PetscScalar[grid.p->Mz];
  for (PetscInt k=0; k < grid.p->Mz; k++)   izz[k] = ((PetscScalar) k) * grid.p->dz;
  u = new PetscScalar[grid.p->Mz];
  v = new PetscScalar[grid.p->Mz];
  w = new PetscScalar[grid.p->Mz];

  ierr = u3.needAccessToVals(); CHKERRQ(ierr);
  ierr = v3.needAccessToVals(); CHKERRQ(ierr);
  ierr = w3.needAccessToVals(); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuplift, &dbdt); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
    
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // basal w from basal kinematical equation
      const PetscScalar dbdx = (b[i+1][j] - b[i-1][j]) / (2.0*dx),
                        dbdy = (b[i][j+1] - b[i][j-1]) / (2.0*dy);
//      w[0] = dbdt[i][j] + ub[i][j] * dbdx + vb[i][j] * dbdy;
      w[0] = ub[i][j] * dbdx + vb[i][j] * dbdy;  // DEBUG: remove dbdt
      if (includeBMRinContinuity == PETSC_TRUE) {
        w[0] -= capBasalMeltRate(basalMeltRate[i][j]);
      }

      ierr = u3.getValColumn(i,j,grid.p->Mz,izz,u); CHKERRQ(ierr);
      ierr = v3.getValColumn(i,j,grid.p->Mz,izz,v); CHKERRQ(ierr);

      // compute w above base by trapezoid rule 
      planeStar uss, vss;
      ierr = u3.getPlaneStarZ(i,j,0.0,&uss);
      ierr = v3.getPlaneStarZ(i,j,0.0,&vss);
      PetscScalar OLDintegrand = (uss.ip1 - uss.im1) / (2.0*dx) + (vss.jp1 - vss.jm1) / (2.0*dy);
      // at bottom, difference up:
      OLDintegrand -= dbdx * (u[1] - u[0]) / dz + dbdy * (v[1] - v[0]) / dz;

      for (PetscInt k = 1; k < Mz; ++k) {
        const PetscScalar zk = k * dz;
        ierr = u3.getPlaneStarZ(i,j,zk,&uss);
        ierr = v3.getPlaneStarZ(i,j,zk,&vss);
        PetscScalar NEWintegrand = (uss.ip1 - uss.im1) / (2.0*dx) + (vss.jp1 - vss.jm1) / (2.0*dy);
        if (k == Mz-1) { // at top, difference down:
          NEWintegrand -= dbdx * (u[k] - u[k-1]) / dz + dbdy * (v[k] - v[k-1]) / dz;
        } else { // usual case; central difference
          NEWintegrand -= dbdx * (u[k+1] - u[k-1]) / (2.0*dz) + dbdy * (v[k+1] - v[k-1]) / (2.0*dz);
        }
        w[k] = w[k-1] - 0.5 * (NEWintegrand + OLDintegrand) * dz;
        OLDintegrand = NEWintegrand;
      }
      
      ierr = w3.setValColumn(i,j,grid.p->Mz,izz,w); CHKERRQ(ierr);      
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuplift, &dbdt); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
  
  ierr = u3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = v3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = w3.doneAccessToVals(); CHKERRQ(ierr);

  delete [] izz;  delete [] u;  delete [] v;  delete [] w;
  return 0;
}


PetscScalar IceModel::capBasalMeltRate(const PetscScalar bMR) {
  const PetscScalar MAX_BASALMELTRATE = 0.5/secpera;
  if (bMR > MAX_BASALMELTRATE) {
    return MAX_BASALMELTRATE;
  } else if (bMR < -MAX_BASALMELTRATE) {
    return -MAX_BASALMELTRATE;
  } else {
    return bMR;
  }
}


/* 
THIS IS THE ONLY PLACE WHERE STENCIL_*BOX* COMMUNICATION IS REQUIRED FOR A 3D Vec IS REQUIRED
REMOVING IT THEREFORE HAS A BENEFIT FOR COMMUNICATION

PetscErrorCode IceModel::smoothSigma() {
  // does iterated smoothing of Sigma on regular grid; uses box stencil of width one
  // not exactly recommended; it is rather nonphysical; compare Bueler, Brown, and Lingle 2007
  PetscErrorCode  ierr;
  PetscScalar ***Snew, **H;
  // following is [0.35 0.30 0.35] outer product w itself; ref conversation with Orion L.:
  const PetscScalar c[3][3] = {{0.1225, 0.105,  0.1225},
                               {0.105,  0.09,   0.105},
                               {0.1225, 0.105,  0.1225}};

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);

  for (int count=0; count < noSpokesLevel; ++count) {
    ierr = Sigma3.needAccessToVals(); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vWork3d[0], &Snew); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        for (PetscInt k=0; k<ks; ++k) {
          // note Sigma[neighbor][neighbor][k] will be zero if outside of ice
          planeBox Sbb;
          ierr = Sigma3.getPlaneBoxZ(i,j,k * grid.p->dz,&Sbb);
          const PetscScalar SS =  
              c[0][0]*Sbb.im1jm1 + c[0][1]*Sbb.im1 + c[0][2]*Sbb.im1jp1
            + c[1][0]*Sbb.jm1    + c[1][1]*Sbb.ij  + c[1][2]*Sbb.jp1
            + c[2][0]*Sbb.ip1jm1 + c[2][1]*Sbb.ip1 + c[2][2]*Sbb.ip1jp1;
          const PetscScalar myz = k * grid.p->dz;
          PetscScalar active = 0.0;  // build sum of coeffs for neighbors with ice at or above curr depth
          if (H[i-1][j-1] >= myz) { active += c[0][0]; } 
          if (H[i-1][j] >= myz) { active += c[0][1]; } 
          if (H[i-1][j+1] >= myz) { active += c[0][2]; } 
          if (H[i][j-1] >= myz) { active += c[1][0]; } 
          if (H[i][j+1] >= myz) { active += c[1][2]; } 
          if (H[i+1][j-1] >= myz) { active += c[2][0]; } 
          if (H[i+1][j] >= myz) { active += c[2][1]; } 
          if (H[i+1][j+1] >= myz) { active += c[2][2]; } 
          Snew[i][j][k] = SS / (active + c[1][1]);  // ensures all heating goes into ice not outside it
        }
        for (PetscInt k=ks+1; k<grid.p->Mz; ++k) {
          Snew[i][j][k] = 0;
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da3, vWork3d[0], &Snew); CHKERRQ(ierr);
    ierr = Sigma3.doneAccessToVals(); CHKERRQ(ierr);

    // communicate ghosted values *and* transfer Snew to vSigma
    ierr = DALocalToLocalBegin(grid.da3, vWork3d[0], INSERT_VALUES, Sigma3.v); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vWork3d[0], INSERT_VALUES, Sigma3.v); CHKERRQ(ierr);
  } // for
  
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  return 0;
}
*/


//! Compute the maximum velocities for time-stepping and reporting to user.
PetscErrorCode IceModel::computeMax3DVelocities() {
  // computes max velocities in 3D grid and also sets CFLmaxdt by CFL condition
  PetscErrorCode ierr;
  PetscScalar **H;
  PetscScalar locCFLmaxdt = maxdt;

  PetscScalar *izz, *u, *v, *w;
  izz = new PetscScalar[grid.p->Mz];
  for (PetscInt k=0; k < grid.p->Mz; k++)   izz[k] = ((PetscScalar) k) * grid.p->dz;
  u = new PetscScalar[grid.p->Mz];
  v = new PetscScalar[grid.p->Mz];
  w = new PetscScalar[grid.p->Mz];

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = u3.needAccessToVals(); CHKERRQ(ierr);
  ierr = v3.needAccessToVals(); CHKERRQ(ierr);
  ierr = w3.needAccessToVals(); CHKERRQ(ierr);

  // update global max of abs of velocities for CFL; only velocities under surface
  PetscReal   maxu=0.0, maxv=0.0, maxw=0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt      ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
      const bool isMarginal = checkThinNeigh(H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                             H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);
      ierr = u3.getValColumn(i,j,grid.p->Mz,izz,u); CHKERRQ(ierr);
      ierr = v3.getValColumn(i,j,grid.p->Mz,izz,v); CHKERRQ(ierr);
      ierr = w3.getValColumn(i,j,grid.p->Mz,izz,w); CHKERRQ(ierr);
      for (PetscInt k=0; k<ks; ++k) {
        const PetscScalar au = PetscAbs(u[k]);
        const PetscScalar av = PetscAbs(v[k]);
        maxu = PetscMax(maxu,au);
        maxv = PetscMax(maxv,av);
        PetscScalar tempdenom = (0.001/secpera)/(grid.p->dx + grid.p->dy);  // make sure it's pos.
        tempdenom += PetscAbs(au/grid.p->dx) + PetscAbs(av/grid.p->dy);
        if (!isMarginal) {
          const PetscScalar aw = PetscAbs(w[k]);
          maxw = PetscMax(maxw,aw);
          tempdenom += PetscAbs(aw/grid.p->dz);
        }
        locCFLmaxdt = PetscMin(locCFLmaxdt,1.0 / tempdenom); 
      }
    }
  }

  ierr = u3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = v3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = w3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  delete [] izz;  delete [] u;  delete [] v;  delete [] w;

  ierr = PetscGlobalMax(&maxu, &gmaxu, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxv, &gmaxv, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxw, &gmaxw, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMin(&locCFLmaxdt, &CFLmaxdt, grid.com); CHKERRQ(ierr);
  return 0;
}


//! Because the map-plane mass continuity is advective in sliding case, compute CFL.
/*!
This procedure computes the maximum horizontal speed in the SSA areas so that
the CFL condition for the upwinding (in massBalExplicitStep() and only for 
basal component of mass flux) can be computed.
 */
PetscErrorCode IceModel::computeMax2DSlidingSpeed() {
  PetscErrorCode ierr;
  PetscScalar **ub, **vb;
  PetscScalar locCFLmaxdt2D = maxdt;
  
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar denom = PetscAbs(ub[i][j])/grid.p->dx + PetscAbs(vb[i][j])/grid.p->dy;
      denom += (0.01/secpera)/(grid.p->dx + grid.p->dy);  // make sure it's pos.
      locCFLmaxdt2D = PetscMin(locCFLmaxdt2D,1.0/denom);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);

  ierr = PetscGlobalMin(&locCFLmaxdt2D, &CFLmaxdt2D, grid.com); CHKERRQ(ierr);
  return 0;
}


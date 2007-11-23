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
    ierr = velocitySIAStaggered(); CHKERRQ(ierr);

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
      // communicate I[o],J[o] on staggered for horizontalVelocitySIARegular()
      //   and also communicate Sigma[o] on staggered for SigmaSIAToRegular()
      for (PetscInt k=0; k<6; ++k) {
        ierr = DALocalToLocalBegin(grid.da3, vWork3d[k], INSERT_VALUES, vWork3d[k]); CHKERRQ(ierr);
        ierr = DALocalToLocalEnd(grid.da3, vWork3d[k], INSERT_VALUES, vWork3d[k]); CHKERRQ(ierr);
      }
      ierr = SigmaSIAToRegular(); CHKERRQ(ierr);
      ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);
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
  if (updateVelocityAtDepth || useSSAVelocity) {  // in latter case u,v are modified by broacastSSAVelocity()
    ierr = DALocalToLocalBegin(grid.da3, vu, INSERT_VALUES, vu); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vu, INSERT_VALUES, vu); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da3, vv, INSERT_VALUES, vv); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vv, INSERT_VALUES, vv); CHKERRQ(ierr);
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
  
  // smoothing Sigma is NOT RECOMMENDED, but do it here
  if (noSpokesLevel > 0) {
    ierr = smoothSigma(); CHKERRQ(ierr); // comm here
  }

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
@cond CONTINUUM
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
@endcond

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

@cond NUMERIC
The vertical integral is computed by the trapezoid rule.
@endcond
 */
PetscErrorCode IceModel::vertVelocityFromIncompressibility() {
  
  const PetscScalar   dx = grid.p->dx, dy = grid.p->dy, dz = grid.p->dz;
  const PetscInt      Mz = grid.p->Mz;
  PetscErrorCode  ierr;
  PetscScalar **ub, **vb, **basalMeltRate, **b, **dbdt;
  PetscScalar ***u, ***v, ***w;

  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
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
//      w[i][j][0] = dbdt[i][j] + ub[i][j] * dbdx + vb[i][j] * dbdy;
      w[i][j][0] = ub[i][j] * dbdx + vb[i][j] * dbdy;  // DEBUG: remove dbdt
      if (includeBMRinContinuity == PETSC_TRUE) {
        w[i][j][0] -= capBasalMeltRate(basalMeltRate[i][j]);
      }
      // compute w above base by trapezoid rule 
      // (note vertical variable is  z = zorig - b(x,y,t))
      PetscScalar OLDintegrand =   (u[i+1][j][0] - u[i-1][j][0]) / (2.0*dx)
                                 + (v[i][j+1][0] - v[i][j-1][0]) / (2.0*dy);
      // at bottom, difference up:
      OLDintegrand -= dbdx * (u[i][j][1] - u[i][j][0]) / dz
                      + dbdy * (v[i][j][1] - v[i][j][0]) / dz;
      for (PetscInt k = 1; k < Mz; ++k) {
        PetscScalar NEWintegrand =   (u[i+1][j][k] - u[i-1][j][k]) / (2.0*dx)
                                   + (v[i][j+1][k] - v[i][j-1][k]) / (2.0*dy);
        if (k == Mz-1) { // at top, difference down:
          NEWintegrand -= dbdx * (u[i][j][k] - u[i][j][k-1]) / dz
                          + dbdy * (v[i][j][k] - v[i][j][k-1]) / dz;
        } else { // usual case; central difference
          NEWintegrand -= dbdx * (u[i][j][k+1] - u[i][j][k-1]) / (2.0*dz)
                          + dbdy * (v[i][j][k+1] - v[i][j][k-1]) / (2.0*dz);
        }
        w[i][j][k] = w[i][j][k-1] - 0.5 * (NEWintegrand + OLDintegrand) * dz;
        OLDintegrand = NEWintegrand;
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuplift, &dbdt); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
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


PetscErrorCode IceModel::smoothSigma() {
  // does iterated smoothing of Sigma on regular grid; uses box stencil of width one
  // not exactly recommended; it is rather nonphysical; compare Bueler, Brown, and Lingle 2007
  PetscErrorCode  ierr;
  PetscScalar ***S, ***Snew, **H;
  // following is [0.35 0.30 0.35] outer product w itself; ref conversation with Orion L.:
  const PetscScalar c[3][3] = {{0.1225, 0.105,  0.1225},
                               {0.105,  0.09,   0.105},
                               {0.1225, 0.105,  0.1225}};

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);

  for (int count=0; count < noSpokesLevel; ++count) {
    ierr = DAVecGetArray(grid.da3, vSigma, &S); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vWork3d[4], &Snew); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        for (PetscInt k=0; k<ks; ++k) {
          // note S[neighbor_i][neighbor_j][k] will be zero if outside of ice
          const PetscScalar SS =  c[0][0]*S[i-1][j-1][k] + c[0][1]*S[i-1][j][k] + c[0][2]*S[i-1][j+1][k]
            + c[1][0]*S[i][j-1][k]   + c[1][1]*S[i][j][k]   + c[1][2]*S[i][j+1][k]
            + c[2][0]*S[i+1][j-1][k] + c[2][1]*S[i+1][j][k] + c[2][2]*S[i+1][j+1][k];
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
    ierr = DAVecRestoreArray(grid.da3, vSigma, &S); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vWork3d[4], &Snew); CHKERRQ(ierr);
      
    // communicate ghosted values *and* transfer Snew to vSigma
    ierr = DALocalToLocalBegin(grid.da3, vWork3d[4], INSERT_VALUES, vSigma); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vWork3d[4], INSERT_VALUES, vSigma); CHKERRQ(ierr);
  } // for
  
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  return 0;
}


//! Compute the maximum velocities for time-stepping and reporting to user.
PetscErrorCode IceModel::computeMax3DVelocities() {
  // computes max velocities in 3D grid and also sets CFLmaxdt by CFL condition
  PetscErrorCode ierr;
  PetscScalar ***u, ***v, ***w, **H;
  PetscScalar locCFLmaxdt = maxdt;

  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);

  // update global max of abs of velocities for CFL; only velocities under surface
  PetscReal   maxu=0.0, maxv=0.0, maxw=0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt      ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
      const bool isMarginal = checkThinNeigh(H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                             H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);
      for (PetscInt k=0; k<ks; ++k) {
        const PetscScalar au = PetscAbs(u[i][j][k]);
        const PetscScalar av = PetscAbs(v[i][j][k]);
        maxu = PetscMax(maxu,au);
        maxv = PetscMax(maxv,av);
        PetscScalar tempdenom = (0.001/secpera)/(grid.p->dx + grid.p->dy);  // make sure it's pos.
        tempdenom += PetscAbs(au/grid.p->dx) + PetscAbs(av/grid.p->dy);
        if (!isMarginal) {
          const PetscScalar aw = PetscAbs(w[i][j][k]);
          maxw = PetscMax(maxw,aw);
          tempdenom += PetscAbs(aw/grid.p->dz);
        }
        locCFLmaxdt = PetscMin(locCFLmaxdt,1.0 / tempdenom); 
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

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


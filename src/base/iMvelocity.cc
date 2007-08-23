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

PetscErrorCode IceModel::velocity(bool updateVelocityAtDepth) {
  PetscErrorCode ierr;
  static PetscTruth firstTime = PETSC_TRUE;

  if (computeSIAVelocities == PETSC_TRUE) {
    ierr = velocitySIAStaggered(!updateVelocityAtDepth); CHKERRQ(ierr);

    ierr = DALocalToLocalBegin(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
    for (PetscInt k=4; k<8; ++k) { // communicate ub[o], vb[o] on staggered
      ierr = DALocalToLocalBegin(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
      ierr = DALocalToLocalEnd(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
    }

    ierr = basalSIAConditionsToRegular(); CHKERRQ(ierr);
  
    if (updateVelocityAtDepth) {
  
      for (PetscInt k=0; k<6; ++k) { // communicate I[o], J[o], Sigma[o] on staggered
        ierr = DALocalToLocalBegin(grid.da3, vWork3d[k], INSERT_VALUES, vWork3d[k]); CHKERRQ(ierr);
        ierr = DALocalToLocalEnd(grid.da3, vWork3d[k], INSERT_VALUES, vWork3d[k]); CHKERRQ(ierr);
      }
      for (PetscInt k=0; k<4; ++k) { // communicate h_x[o], h_y[o] on staggered
        ierr = DALocalToLocalBegin(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
        ierr = DALocalToLocalEnd(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
      }
      for (PetscInt k=8; k<10; ++k) { // communicate Rb[o] on staggered
        ierr = DALocalToLocalBegin(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
        ierr = DALocalToLocalEnd(grid.da2, vWork2d[k], INSERT_VALUES, vWork2d[k]); CHKERRQ(ierr);
      }

      ierr = SigmaSIAToRegular(); CHKERRQ(ierr);
      ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);
      ierr = verticalVelocitySIARegular(); CHKERRQ(ierr);

      if (noSpokesLevel > 0) {
        ierr = smoothSigma(); CHKERRQ(ierr); // communication occurs in here!
      }
    }
  }  // if computeSIAVelocities

  if (useSSAVelocity) { // communication happens within SSA
    ierr = setupForSSA(DEFAULT_MINH_SSA); CHKERRQ(ierr);
    if (firstTime) {
      ierr = vertAveragedVelocityToRegular(); CHKERRQ(ierr);  // comm here
    }
    ierr = velocitySSA(); CHKERRQ(ierr); // comm here ...
    ierr = cleanupAfterSSA(DEFAULT_MINH_SSA); CHKERRQ(ierr);
    ierr = broadcastSSAVelocity(); CHKERRQ(ierr); // sets CFLmaxdt2D
    ierr = correctSigma(); CHKERRQ(ierr);
    ierr = correctBasalFrictionalHeating(); CHKERRQ(ierr);
  } else { // Note vertically averaged vels on regular grid (Vecs vubar and 
     // vvbar) are set by SSA procedures above, including (in a reasonable 
     // manner) within the MASK_SHEET regions.  Also vubar and vvbar are used 
     // to initialize the next step of SSA.  SO we want to NOT overwrite 
     // the just-computed regular grid velocities (as they will be needed to
     // initialize another step of SSA; note first time exception).
     // BUT if we don't call SSA in velocity() then
     // we need to get staggered onto regular for viewer and for summary().
    ierr = vertAveragedVelocityToRegular(); CHKERRQ(ierr); // comm here
  }

  if ((useSSAVelocity) || (updateVelocityAtDepth)) {
    // communication here for global max; sets CFLmaxdt
    ierr = computeMax3DVelocities(); CHKERRQ(ierr); 
  }
  
  firstTime = PETSC_FALSE;
  return 0;
}


PetscErrorCode IceModel::smoothSigma() {
  // does iterated smoothing of Sigma on regular grid; uses box stencil of width one
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


PetscErrorCode IceModel::vertAveragedVelocityToRegular() {
  // only 2D regular grid velocities are updated here
  
  PetscErrorCode ierr;
  PetscScalar **u, **v, **uvbar[2];

  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      u[i][j] = 0.5*(uvbar[0][i-1][j] + uvbar[0][i][j]);
      v[i][j] = 0.5*(uvbar[1][i][j-1] + uvbar[1][i][j]);
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vubar, INSERT_VALUES, vubar); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vvbar, INSERT_VALUES, vvbar); CHKERRQ(ierr);

  return 0;
}


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

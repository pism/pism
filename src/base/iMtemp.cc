// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
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

#include <petscda.h>
#include <petscksp.h>
#include "iceModel.hh"


//! Manages the time-stepping and parallel communication for the temperature and age equations.
/*! 
Note that both the temperature equation and the age equation involve advection and have a CFL
condition (Morton & Mayers 1994).  By being slightly conservative we use the same CFL condition
for both (Bueler and others, 2007. "Exact solutions ... thermomechanically-coupled ...," 
J. Glaciol.).  We also report any CFL violations; these can \em only occur when using the 
<tt>-tempskip</tt> option.
 */
PetscErrorCode IceModel::temperatureAgeStep() {
  // update temp and age fields
  PetscErrorCode  ierr;

  PetscScalar  myCFLviolcount = 0.0;  // it is a count but it is type "PetscScalar"
                                      // because that type works with PetscGlobalSum()
  // do CFL and vertical grid blow-out checking only in ageStep()
  ierr = ageStep(&myCFLviolcount); CHKERRQ(ierr);  // puts vtaunew in vWork3d[1]
    
  ierr = temperatureStep(); CHKERRQ(ierr);  // puts vTnew in vWork3d[0]

  // no communication done in ageStep(), temperatureStep(); all done here:
/*
  ierr = DALocalToLocalBegin(grid.da3, Tnew3.v, INSERT_VALUES, T3.v); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, Tnew3.v, INSERT_VALUES, T3.v); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, taunew3.v, INSERT_VALUES, tau3.v); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, taunew3.v, INSERT_VALUES, tau3.v); CHKERRQ(ierr);
*/
  ierr = T3.beginGhostCommTransfer(Tnew3); CHKERRQ(ierr);
  ierr = tau3.beginGhostCommTransfer(taunew3); CHKERRQ(ierr);
  ierr = T3.endGhostCommTransfer(Tnew3); CHKERRQ(ierr);
  ierr = tau3.endGhostCommTransfer(taunew3); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&myCFLviolcount, &CFLviolcount, grid.com); CHKERRQ(ierr);

  return 0;
}


// documentation for temperatureStep() is in pism/src/base/comments.hh
PetscErrorCode IceModel::temperatureStep() {
  // update temp fields Tnew, Tb
  
  // method uses equally-spaced calculation but the methods getValColumn(), setValColumn()
  //   interpolate back and forth between (usually) non-equally space storage grid
  if (grid.dzEQ != grid.dzbEQ) {
    SETERRQ2(9, "IceModel::temperatureStep() method requires grid.dzEQ=%14.11f exactly equal to\n"
                "  grid.dzbEQ=%14.11f\n",
                grid.dzEQ, grid.dzbEQ);
  }

  PetscErrorCode  ierr;

  const PetscInt      Mbz = grid.p->Mbz, 
                      Mz = grid.p->Mz,
                      k0 = Mbz - 1;

  const PetscScalar   dx = grid.p->dx, 
                      dy = grid.p->dy;

  const PetscScalar   rho_c_I = ice.rho * ice.c_p;
  const PetscScalar   rho_c_br = bedrock.rho * bedrock.c_p;
  const PetscScalar   rho_c_av = (rho_c_I + rho_c_br) / 2.0;// only applies for equal-spaced!
  const PetscScalar   iceK = ice.k / rho_c_I;
  const PetscScalar   iceR = iceK * dtTempAge / PetscSqr(grid.dzEQ);
  const PetscScalar   brK = bedrock.k / rho_c_br;
  const PetscScalar   brR = brK * dtTempAge / PetscSqr(grid.dzEQ);

  PetscScalar *Tb, *Tbnew;
  PetscScalar **Ts, **H, **b, **Ghf, **mask, **Hmelt, **Rb, **basalMeltRate;

  PetscScalar *u, *v, *w, *Sigma, *T, *Tnew;
  u = new PetscScalar[grid.p->Mz];
  v = new PetscScalar[grid.p->Mz];
  w = new PetscScalar[grid.p->Mz];
  Sigma = new PetscScalar[grid.p->Mz];
  T = new PetscScalar[grid.p->Mz];
  Tnew = new PetscScalar[grid.p->Mz];

  Tbnew = new PetscScalar[grid.p->Mbz];

  PetscScalar *Lp, *L, *D, *U, *x, *rhs, *work;  
  Lp = new PetscScalar[Mz+k0-1]; L = Lp-1; // ptr arith.; note L[0]=Lp[-1] not alloc
  D = new PetscScalar[Mz+k0];
  U = new PetscScalar[Mz+k0-1];
  x = new PetscScalar[Mz+k0];
  rhs = new PetscScalar[Mz+k0];
  work = new PetscScalar[Mz+k0];

  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vGhf, &Ghf); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vRb, &Rb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
  
  ierr = u3.needAccessToVals(); CHKERRQ(ierr);
  ierr = v3.needAccessToVals(); CHKERRQ(ierr);
  ierr = w3.needAccessToVals(); CHKERRQ(ierr);
  ierr = Sigma3.needAccessToVals(); CHKERRQ(ierr);
  ierr = T3.needAccessToVals(); CHKERRQ(ierr);
  ierr = Tnew3.needAccessToVals(); CHKERRQ(ierr);

  ierr = Tb3.needAccessToVals(); CHKERRQ(ierr);

  PetscInt        myLowTempCount = 0;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt  ks = grid.kBelowHeightEQ(H[i][j]);

      // if isMarginal then only do vertical conduction for ice (i.e. ignor advection
      // and strain heating if isMarginal)
      const bool isMarginal = checkThinNeigh(H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                             H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);
      
      ierr = Tb3.getInternalColumn(i,j,&Tb); CHKERRQ(ierr);

      if (Mbz > 1) { // bedrock present: build k=0:Mbz-1 eqns
        /*
        THIS OLD SCHEME SEEMS ONLY TO GIVE   O(\Delta t,\Delta z^1)  convergence:
        // basal condition (at bottom of bedrock when temp is modelled in bedrock)
        // {from FV (finite volume) across bedrock base z_0}
        D[0] = (1 + brR);  U[0] = -brR;  // note L[0] not an allocated location
        rhs[0] = Tb[i][j][0] + dtTempAge * Ghf[i][j] / (rho_c_br * dz);
        */
        // gives O(\Delta t,\Delta z^2) convergence in Test K
        D[0] = (1.0 + 2.0 * brR);  U[0] = - 2.0 * brR;  // note L[0] not an allocated location
        rhs[0] = Tb[0] + 2.0 * dtTempAge * Ghf[i][j] / (rho_c_br * grid.dzbEQ);
      
        // bedrock only: pure vertical conduction problem
        // {from generic bedrock FV}
        for (PetscInt k=1; k < k0; k++) {
          L[k] = -brR; D[k] = 1+2*brR; U[k] = -brR;
          rhs[k] = Tb[k];
        }
      }
      
      ierr = u3.getValColumn(i,j,Mz,grid.zlevelsEQ,u); CHKERRQ(ierr);
      ierr = v3.getValColumn(i,j,Mz,grid.zlevelsEQ,v); CHKERRQ(ierr);
      ierr = w3.getValColumn(i,j,Mz,grid.zlevelsEQ,w); CHKERRQ(ierr);
      ierr = Sigma3.getValColumn(i,j,Mz,grid.zlevelsEQ,Sigma); CHKERRQ(ierr);
      ierr = T3.getValColumn(i,j,Mz,grid.zlevelsEQ,T); CHKERRQ(ierr);
      
      // bottom part of ice (and top of bedrock in some cases): k=Mbz eqn
      if (ks == 0) { // no ice; set T[0] to surface temp if grounded
        if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
        D[k0] = 1.0; U[k0] = 0.0;
        // if floating and no ice then worry only about bedrock temps;
        // top of bedrock sees ocean
        const PetscScalar floating_base = - (ice.rho/ocean.rho) * H[i][j];
        if (b[i][j] < floating_base - 1.0) {
          rhs[k0] = ice.meltingTemp;
        } else { // top of bedrock sees atmosphere
          rhs[k0] = Ts[i][j];
        }
      } else { // ks > 0; there is ice
        const PetscScalar dz = grid.dzEQ;
        planeStar ss;
        ierr = T3.getPlaneStarZ(i,j,0.0,&ss);
        const PetscScalar UpTu = (u[0] < 0) ? u[0] * (ss.ip1 -  ss.ij) / dx :
                                              u[0] * (ss.ij  - ss.im1) / dx;
        const PetscScalar UpTv = (v[0] < 0) ? v[0] * (ss.jp1 -  ss.ij) / dy :
                                              v[0] * (ss.ij  - ss.jm1) / dy;
        // for w, always upwind *up* from base
        const PetscScalar UpTw = w[0] * (T[1] - T[0]) / dz;
        if (modMask(mask[i][j]) == MASK_FLOATING) {
          // at base of ice shelf, set T = Tpmp but also determine dHmelt/dt
          // by ocean flux; note volume for which energy is being computed is 
          // *half* a segment
          if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
          D[k0] = 1.0 + 2.0 * iceR; U[k0] = - 2.0 * iceR;
          rhs[k0] = T[0] + 2.0 * dtTempAge * DEFAULT_OCEAN_HEAT_FLUX / (rho_c_I * dz);
          if (!isMarginal) {
            rhs[k0] += dtTempAge * (Sigma[0] - UpTu - UpTv - UpTw) / 2;
          }
        } else { // there is *grounded* ice; ice/bedrock interface; from FV across interface
          const PetscScalar rho_c_ratio = rho_c_I / rho_c_av;
          rhs[k0] = T[0] + dtTempAge * (Rb[i][j] / (rho_c_av * dz));
          if (!isMarginal) {
            rhs[k0] += dtTempAge * rho_c_ratio * 0.5 * Sigma[0];
            rhs[k0] -= dtTempAge * rho_c_ratio * (0.5 * (UpTu + UpTv + UpTw) + T[0] * w[0] / dz);
          }
          const PetscScalar iceReff = ice.k * dtTempAge / (rho_c_av * dz * dz);
          const PetscScalar brReff = bedrock.k * dtTempAge / (rho_c_av * dz * dz);
          if (Mbz > 1) { // there is bedrock; apply centered difference with 
                         // jump in diffusivity coefficient
            L[k0] = - brReff; D[k0] = 1 + iceReff + brReff; U[k0] = - iceReff;
          } else { // no bedrock; apply geothermal flux here
            // L[k0] = 0.0;  (note this is not an allocated location!) 
            D[k0] = 1.0 + 2.0 * iceR; U[k0] = - 2.0 * iceR;
            rhs[k0] += 2.0 * dtTempAge * Ghf[i][j] / (rho_c_I * dz);
          }
        }
      }

      // generic ice segment: build k0+1:k0+ks-1 eqns
      for (PetscInt k = 1; k < ks; k++) {
        const PetscScalar dz = grid.dzEQ;
        planeStar ss;
        ierr = T3.getPlaneStarZ(i,j,k * dz,&ss);
        const PetscScalar UpTu = (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx :
                                              u[k] * (ss.ij  - ss.im1) / dx;
        const PetscScalar UpTv = (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy :
                                              v[k] * (ss.ij  - ss.jm1) / dy;
        const PetscScalar UpTw = (w[k] < 0) ? w[k] * (T[k+1] -   T[k]) / dz :
                                              w[k] * (T[k]   - T[k-1]) / dz;
        L[k0+k] = -iceR; D[k0+k] = 1+2*iceR; U[k0+k] = -iceR;
        rhs[k0+k] = T[k];
        if (!isMarginal) {
          rhs[k0+k] += dtTempAge * (Sigma[k] - UpTu - UpTv - UpTw);
        }
      }
      
      // surface b.c.
      if (k0+ks>0) {
        L[k0+ks] = 0;   D[k0+ks] = 1.0;   // ignor U[k0+ks]
        rhs[k0+ks] = Ts[i][j];
        //  HAD NO k0+ks eqn before, and:
        //        rhs[k0+ks-1] += iceR * Ts[i][j];
        // U[k0+ks-1] = 0.0, but never actually eval'ed by tridiag solve
      }

      // solve system; melting not addressed yet
      if (k0+ks>0) {
        ierr = solveTridiagonalSystem(L, D, U, x, rhs, work, k0+ks+1);
        // OLD:       ierr = solveTridiagonalSystem(L, D, U, x, rhs, work, k0+ks);
        if (ierr != 0) {
          SETERRQ3(1, "Tridiagonal solve failed at (%d,%d) with zero pivot in position %d.",
               i, j, ierr);
        }
      }

      // insert bedrock solution; check for too low below
      for (PetscInt k=0; k < k0; k++) {
        Tbnew[k] = x[k];
      }

      // prepare for melting/refreezing
      PetscScalar Hmeltnew = Hmelt[i][j];
      
      // insert solution for generic ice segments
      for (PetscInt k=1; k <= ks; k++) {
//      for (PetscInt k=1; k < ks; k++) {
        if (allowAboveMelting == PETSC_TRUE) {
          Tnew[k] = x[k0 + k];
        } else {
          const PetscScalar dz = grid.dzEQ;
          const PetscScalar depth = H[i][j] - k * dz;
          const PetscScalar Tpmp = ice.meltingTemp - ice.beta_CC_grad * depth;
          if (x[k0 + k] > Tpmp) {
            Tnew[k] = Tpmp;
            PetscScalar Texcess = x[k0 + k] - Tpmp; // always positive
            excessToFromBasalMeltLayer(rho_c_I, k * dz, &Texcess, &Hmeltnew);
            // Texcess  will always come back zero here; ignor it
          } else {
            Tnew[k] = x[k0 + k];
          }
        }
        if (Tnew[k] < globalMinAllowedTemp) {
           ierr = PetscPrintf(PETSC_COMM_SELF,
              "  [[too low (<200) generic segment temp T = %f at %d,%d,%d; proc %d; mask=%f; w=%f]]\n",
              Tnew[k],i,j,k,grid.rank,mask[i][j],w[k]*secpera); CHKERRQ(ierr);
           myLowTempCount++;
        }
      }
      
      // insert solution for ice/rock interface (or base of ice shelf) segment
      if (ks > 0) {
        if (allowAboveMelting == PETSC_TRUE) {
          Tnew[0] = x[k0];
        } else {  // compute diff between x[k0] and Tpmp; melt or refreeze as appropriate
          const PetscScalar Tpmp = ice.meltingTemp - ice.beta_CC_grad * H[i][j];
          PetscScalar Texcess = x[k0] - Tpmp; // positive or negative
          if (modMask(mask[i][j]) == MASK_FLOATING) {
             // when floating, only half a segment has had its temperature raised
             // above Tpmp
             excessToFromBasalMeltLayer(rho_c_I/2, 0.0, &Texcess, &Hmeltnew);
          } else {
             excessToFromBasalMeltLayer(rho_c_av, 0.0, &Texcess, &Hmeltnew);
          }
          Tnew[0] = Tpmp + Texcess;
          if (Tnew[0] > (Tpmp + 0.00001)) {
            SETERRQ(1,"updated temperature came out above Tpmp");
          }
        }
        if (Tnew[0] < globalMinAllowedTemp) {
           ierr = PetscPrintf(PETSC_COMM_SELF,
              "  [[too low (<200) ice/rock segment temp T = %f at %d,%d; proc %d; mask=%f; w=%f]]\n",
              Tnew[0],i,j,grid.rank,mask[i][j],w[0]*secpera); CHKERRQ(ierr);
           myLowTempCount++;
        }
      } else {
        Hmeltnew = 0.0;
      }
      
      // we must agree on redundant values T(z=0) at top of bedrock and at bottom of ice
      if (ks > 0) {
        Tbnew[k0] = Tnew[0];
      } else {
        // if floating then top of bedrock sees ocean
        const PetscScalar floating_base = - (ice.rho/ocean.rho) * H[i][j];
        if (b[i][j] < floating_base - 1.0) {
          Tbnew[k0] = ice.meltingTemp;
        } else { // top of bedrock sees atmosphere
          Tbnew[k0] = Ts[i][j];
        }
      }
      // check bedrock solution        
      for (PetscInt k=0; k <= k0; k++) {
        if (Tbnew[k] < globalMinAllowedTemp) {
           ierr = PetscPrintf(PETSC_COMM_SELF,
              "  [[too low (<200) bedrock temp T = %f at %d,%d,%d; proc %d; mask=%f]]\n",
              Tbnew[k],i,j,k,grid.rank,mask[i][j]); CHKERRQ(ierr);
           myLowTempCount++;
        }
      }

      // transfer column into Tb3; neighboring columns will not reference!
      ierr = Tb3.setInternalColumn(i,j,Tbnew); CHKERRQ(ierr);

      // set to air temp above ice
      for (PetscInt k=ks; k<Mz; k++) {
        Tnew[k] = Ts[i][j];
      }

      // transfer column into Tnew3; communication later
      ierr = Tnew3.setValColumn(i,j,Mz,grid.zlevelsEQ,Tnew); CHKERRQ(ierr);

      // basaMeltRate[][] is rate of change of Hmelt[][]; thus it can be negative
      basalMeltRate[i][j] = (Hmeltnew - Hmelt[i][j]) / dtTempAge;

      // limit Hmelt by default max
      Hmeltnew = PetscMin(DEFAULT_MAX_HMELT,Hmeltnew);

      // eliminate basal water if floating
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        Hmelt[i][j] = 0.0;
      } else {
        Hmelt[i][j] = Hmeltnew;
      }

    } 
  }
  
  if (myLowTempCount > maxLowTempCount) { SETERRQ(1,"too many low temps"); }

  // note that in above code 4 scalar fields were modified: vHmelt, vbasalMeltRate, vTb, and vT
  // but vHmelt, vbasalMeltRate and vTb will never need to communicate ghosted values
  // (i.e. horizontal stencil neighbors);  vT is communicated by temperatureAgeStep()

  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vGhf, &Ghf); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vRb, &Rb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);

  ierr = Tb3.needAccessToVals(); CHKERRQ(ierr);

  ierr = u3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = v3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = w3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = Sigma3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = T3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = Tnew3.doneAccessToVals(); CHKERRQ(ierr);
  
  delete [] Lp; delete [] D; delete [] U; delete [] x; delete [] rhs; delete [] work;
  delete [] u;  delete [] v;  delete [] w;  delete [] Sigma;  delete [] T;
  delete [] Tbnew;  delete [] Tnew;

  return 0;
}


//! Compute the melt water which should go to the base if \f$T\f$ is above pressure-melting.
PetscErrorCode IceModel::excessToFromBasalMeltLayer(
                PetscScalar rho_c, PetscScalar z,
                PetscScalar *Texcess, PetscScalar *Hmelt) {

  const PetscScalar darea = grid.p->dx * grid.p->dy;
  const PetscScalar dvol = darea * grid.dzEQ;
  const PetscScalar dE = rho_c * (*Texcess) * dvol;
  const PetscScalar massmelted = dE / ice.latentHeat;

  if (allowAboveMelting == PETSC_TRUE) {
    SETERRQ(1,"excessToBasalMeltLayer() called but allowAboveMelting==TRUE");
  }
  if (*Texcess >= 0.0) {
    if (updateHmelt == PETSC_TRUE) {
      // T is at or above pressure-melting temp, so temp needs to be set to 
      // pressure-melting, and a fraction of excess energy
      // needs to be turned into melt water at base
      // note massmelted is POSITIVE!
//      const PetscScalar FRACTION_TO_BASE
//                           = (z < 100.0) ? 0.4 * (100.0 - z) / 100.0 : 0.0;
      const PetscScalar FRACTION_TO_BASE
                           = (z < 100.0) ? 0.2 * (100.0 - z) / 100.0 : 0.0;
      *Hmelt += (FRACTION_TO_BASE * massmelted) / (ice.rho * darea);  // note: ice-equiv thickness
    }
    *Texcess = 0.0;
  } else if (updateHmelt == PETSC_TRUE) {  // neither Texcess nor Hmelt need to change 
                                           // if Texcess < 0.0
    // Texcess negative; only refreeze (i.e. reduce Hmelt) if at base and Hmelt > 0.0
    // note ONLY CALLED IF AT BASE!   note massmelted is NEGATIVE!
    if (z > 0.00001) {
      SETERRQ(1, "excessToBasalMeltLayer() called with z not at base and negative Texcess");
    }
    if (*Hmelt > 0.0) {
      const PetscScalar thicknessToFreezeOn = - massmelted / (ice.rho * darea);
      if (thicknessToFreezeOn <= *Hmelt) { // the water *is* available to freeze on
        *Hmelt -= thicknessToFreezeOn;
        *Texcess = 0.0;
      } else { // only refreeze Hmelt thickness of water; update Texcess
        *Hmelt = 0.0;
        const PetscScalar dTemp = ice.latentHeat * ice.rho * (*Hmelt) / (rho_c * grid.dzEQ);
        *Texcess += dTemp;
      }
    } 
    // note: if *Hmelt == 0 and Texcess < 0.0 then Texcess unmolested; temp will go down
  }
  return 0;
}                           


//! Take an explicit time-step for the age equation.  Also check the CFL for advection.
/*! 
@cond CONTINUUM
The age equation is\f$d\tau/dt = 1\f$, that is,
    \f[ \frac{\partial \tau}{\partial t} + u \frac{\partial \tau}{\partial x} + v \frac{\partial \tau}{\partial y} + w \frac{\partial \tau}{\partial z} = 1\f]
where \f$\tau(t,x,y,z)\f$ is the age of the ice and \f$(u,v,w)\f$  is the three dimensional
velocity field.  This equation is hyperbolic (purely advective).  
The boundary condition is that when the ice fell as snow it had age zero.  
That is, \f$\tau(t,x,y,h(t,x,y)) = 0\f$ in accumulation areas, while there is no 
boundary condition elsewhere (as the characteristics go outward elsewhere).  At this point 
the refreeze case, either grounded basal ice or marine basal ice, is not handled correctly.

By default, when computing the grain size for the Goldsby-Kohlstedt flow law, the age \f$\tau\f$ is not used.
Instead a pseudo age is computed by updateGrainSizeNow().  If you want the age computed by this routine to be
used for the grain size estimation, from the Vostok core relation as in grainSizeVostok(), add 
option <tt>-real_age_grainsize</tt>.
@endcond

@cond NUMERIC
The numerical method is first-order upwind.
@endcond
 */
PetscErrorCode IceModel::ageStep(PetscScalar* CFLviol) {
  // update age field taunew3

  // method uses equally-spaced calculation but the methods getValColumn(), setValColumn()
  //   interpolate back and forth between (usually) non-equally space storage grid

  PetscErrorCode  ierr;

  const PetscInt      Mz = grid.p->Mz;
  const PetscScalar   dx = grid.p->dx, 
                      dy = grid.p->dy, 
                      dz = grid.dzEQ;
  const PetscScalar   cflx = dx/dtTempAge, cfly = dy/dtTempAge, cflz = dz/dtTempAge;

  PetscScalar **H;

  PetscScalar *tau, *u, *v, *w, *taunew;
  tau = new PetscScalar[grid.p->Mz];
  u = new PetscScalar[grid.p->Mz];
  v = new PetscScalar[grid.p->Mz];
  w = new PetscScalar[grid.p->Mz];
  taunew = new PetscScalar[grid.p->Mz];
  
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = tau3.needAccessToVals(); CHKERRQ(ierr);
  ierr = u3.needAccessToVals(); CHKERRQ(ierr);
  ierr = v3.needAccessToVals(); CHKERRQ(ierr);
  ierr = w3.needAccessToVals(); CHKERRQ(ierr);
  ierr = taunew3.needAccessToVals(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt  ks = grid.kBelowHeightEQ(H[i][j]);
/*
      if (ks >= Mz) {
        ierr = PetscPrintf(grid.com,
                 "[[error LOCATION: i, j, ks, H = %5d %5d %5d %10.2f]]\n",i, j, ks, H[i][j]); 
        SETERRQ(1,"Vertical grid exceeded");
      }
*/
    
      // only effects of this is whether vertical velocities are used in advection
      const bool isMarginal = checkThinNeigh(H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                             H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);

      ierr = tau3.getValColumn(i,j,grid.p->Mz,grid.zlevelsEQ,tau); CHKERRQ(ierr);
      ierr = u3.getValColumn(i,j,grid.p->Mz,grid.zlevelsEQ,u); CHKERRQ(ierr);
      ierr = v3.getValColumn(i,j,grid.p->Mz,grid.zlevelsEQ,v); CHKERRQ(ierr);
      ierr = w3.getValColumn(i,j,grid.p->Mz,grid.zlevelsEQ,w); CHKERRQ(ierr);
      for (PetscInt k=0; k<ks; k++) {
        // age evolution is pure advection (so provides check on temp calculation)
        // check CFL conditions at each point, then upwind for age
        if (PetscAbs(u[k]) > cflx)  *CFLviol += 1.0;
        if (PetscAbs(v[k]) > cfly)  *CFLviol += 1.0;
        
        // note ss.ij = tau[k]
        planeStar ss;
        const PetscScalar zk = grid.zlevelsEQ[k];
        ierr = tau3.getPlaneStarZ(i,j,zk,&ss);

        PetscScalar     rtau;
        rtau =  (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx
                           : u[k] * (ss.ij  - ss.im1) / dx;
        rtau += (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy
                           : v[k] * (ss.ij  - ss.jm1) / dy;
        if (!isMarginal) {
          if (PetscAbs(w[k]) > cflz)  *CFLviol += 1.0;
          // if w upward at k=0 then ignor contribution to age
          if ((k > 0) || (w[k] < 0)) {
            // note: this is lowest-order upwinding
            rtau += (w[k] < 0) ? w[k] * (tau[k+1] - tau[k]) / dz
                               : w[k] * (tau[k] - tau[k-1]) / dz;
          }
        }
        taunew[k] = tau[k] + dtTempAge * (1.0 - rtau);
      }      
      for (PetscInt k=ks; k<Mz; k++) {
        taunew[k] = 0.0;  // age of ice above (and at) surface is zero years
      }
      
      ierr = taunew3.setValColumn(i,j,grid.p->Mz,grid.zlevelsEQ,taunew); CHKERRQ(ierr);
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = tau3.doneAccessToVals();  CHKERRQ(ierr);
  ierr = u3.doneAccessToVals();  CHKERRQ(ierr);
  ierr = v3.doneAccessToVals();  CHKERRQ(ierr);
  ierr = w3.doneAccessToVals();  CHKERRQ(ierr);
  ierr = taunew3.doneAccessToVals();  CHKERRQ(ierr);

  delete [] tau;  delete [] u;  delete [] v;  delete [] w;  delete [] taunew;
  
  return 0;
}


bool IceModel::checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                              PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE) {
  const PetscScalar THIN = 100.0;  // thin = (at most 100m thick)
  return (   (E < THIN) || (NE < THIN) || (N < THIN) || (NW < THIN)
          || (W < THIN) || (SW < THIN) || (S < THIN) || (SE < THIN) );
}


PetscErrorCode IceModel::solveTridiagonalSystem(
         const PetscScalar* L, const PetscScalar* D, const PetscScalar* U,
         PetscScalar* x, const PetscScalar* r, PetscScalar* a, const int n) const {
  // modified slightly from Numerical Recipes version

  PetscScalar b;
  b = D[0];
  if (b == 0.0) { return 1; }
  x[0] = r[0]/b;
  for (int i=1; i<n; ++i) {
    a[i] = U[i-1]/b;
    b = D[i] - L[i]*a[i];
    if (b == 0) { return i+1; }
    x[i] = (r[i] - L[i]*x[i-1]) / b;
  }
  for (int i=n-2; i>=0; --i) {
    x[i] -= a[i+1] * x[i+1];
  }

  return 0;
}


// test IceModelVec: assuming this is called from pisms, try
//   pisms -eisII A -y 1 -Mz 11    # no errors when grid coincides; significant otherwise
//   pisms -eisII A -y 1 -Mz 101   # no errors when grid coincides; small otherwise
//   pisms -eisII A -y 1 -Mz 501   # no errors
//   pisms -eisII A -y 1 -Mz 500   # small errors (appropriate; from linear interpolation)
PetscErrorCode IceModel::testIceModelVec()    {
  PetscErrorCode ierr;
  IceModelVec3 test3;

  ierr = test3.create(grid,"testMe",false); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,"\n\ntesting IceModelVec3; setting to constant %f",
                    60402.70804); CHKERRQ(ierr);
  ierr = test3.setToConstant(60402.70804); CHKERRQ(ierr);

  ierr = test3.needAccessToVals(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,"\n\nIceModelVec3::getValZ() says value is %f",
                    test3.getValZ(grid.xs,grid.ys,0.0) ); CHKERRQ(ierr);
  ierr = test3.doneAccessToVals(); CHKERRQ(ierr);

  ierr = test3.beginGhostComm(); CHKERRQ(ierr);
  ierr = test3.endGhostComm(); CHKERRQ(ierr);

  ierr = test3.needAccessToVals(); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,"\n\ntesting IceModelVec3::setValColumn() and getValColumn()\n");
    CHKERRQ(ierr);
  PetscScalar levels[10] = {0.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, -1.0};
  levels[9] = grid.p->Lz;
  PetscScalar valsIN[10], valsOUT[10];
  for (PetscInt k=0; k < 10; k++) {
    valsIN[k] = sin(levels[k]/1000.0);
    //valsIN[k] = levels[k];
  }
  ierr = test3.setValColumn(grid.xs, grid.ys, 10, levels, valsIN); CHKERRQ(ierr);

  /*
  PetscInt myMz;
  PetscScalar *mylev, *myval;
  ierr = test3.getInternalColumn(grid.xs, grid.ys, &myMz, &mylev, &myval); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,"from getInternalColumn(), myMz=%d; returned levels and values are\n",
                    myMz); CHKERRQ(ierr);
  for (PetscInt k=0; k < myMz; k++) {
    ierr = verbPrintf(1,grid.com,"   k=%d:   z[k]=%7.2f   test3[i][j][k]=%7.2f\n",
                      k,mylev[k],myval[k]); CHKERRQ(ierr);
  }
  */
  
  ierr = test3.getValColumn(grid.xs, grid.ys, 10, levels, valsOUT); CHKERRQ(ierr);
  for (PetscInt k=0; k < 10; k++) {
    ierr = verbPrintf(1,grid.com,"   k=%d:   level=%7.2f   valsIN=%7.4f   valsOUT=%7.4f   |diff|=%5.4e\n",
                      k,levels[k],valsIN[k],valsOUT[k],PetscAbs(valsIN[k]-valsOUT[k]) ); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com,"done testing IceModelVec3::setValColumn() and getValColumn()\n\n\n");
    CHKERRQ(ierr);
  ierr = test3.doneAccessToVals(); CHKERRQ(ierr);

  ierr = test3.destroy(); CHKERRQ(ierr);
  return 0;
}


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
  ierr = DALocalToLocalBegin(grid.da3, vWork3d[0], INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vWork3d[0], INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vWork3d[1], INSERT_VALUES, vtau); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vWork3d[1], INSERT_VALUES, vtau); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&myCFLviolcount, &CFLviolcount, grid.com); CHKERRQ(ierr);

  return 0;
}


// documentation for temperatureStep() is in pism/src/base/comments.hh
PetscErrorCode IceModel::temperatureStep() {
  // update temp fields vTnew, vTb
  
  PetscErrorCode  ierr;

  const PetscInt      Mbz = grid.p->Mbz, Mz = grid.p->Mz;
  const PetscScalar   dx = grid.p->dx, dy = grid.p->dy, dz = grid.p->dz;
  const PetscScalar   rho_c_I = ice.rho * ice.c_p;
  const PetscScalar   rho_c_br = bedrock.rho * bedrock.c_p;
  const PetscScalar   rho_c_av = (rho_c_I + rho_c_br) / 2.0;// only applies for equal-spaced!
  const PetscScalar   iceK = ice.k / rho_c_I;
  const PetscScalar   iceR = iceK * dtTempAge / PetscSqr(dz); // only equal-spaced
  const PetscScalar   brK = bedrock.k / rho_c_br;
  const PetscScalar   brR = brK * dtTempAge / PetscSqr(dz); // only equal-spaced

  Vec     vTnew = vWork3d[0];  // will be communicated by temperatureAgeStep()
  PetscScalar **Ts, **H, **b, **Ghf, **mask, **Hmelt, **Rb, **basalMeltRate;
  PetscScalar ***T, ***Tb, ***Tnew, ***u, ***v, ***w, ***Sigma;

  const PetscInt  k0 = Mbz - 1;
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
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vTnew, &Tnew); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);

  PetscInt        myLowTempCount = 0;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/dz));
      // go ahead and assume ks < Mz; see ageStep() for check

      // if isMarginal then only do vertical conduction for ice (i.e. ignor advection
      // and strain heating if isMarginal)
      const bool isMarginal = checkThinNeigh(H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                             H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);
      
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
        rhs[0] = Tb[i][j][0] + 2.0 * dtTempAge * Ghf[i][j] / (rho_c_br * dz);
      
        // bedrock only: pure vertical conduction problem
        // {from generic bedrock FV}
        for (PetscInt k=1; k < k0; k++) {
          L[k] = -brR; D[k] = 1+2*brR; U[k] = -brR;
          rhs[k] = Tb[i][j][k];
        }
      }
      
      // bottom part of ice (and top of bedrock in some cases): k=Mbz eqn
      if (ks == 0) { // no ice; set T[i][j][0] to surface temp if grounded
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
        if (modMask(mask[i][j]) == MASK_FLOATING) {
          // at base of ice shelf, set T = Tpmp but also determine dHmelt/dt
          // by ocean flux; note volume for which energy is being computed is 
          // *half* a segment
          if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
          D[k0] = 1.0 + 2.0 * iceR; U[k0] = - 2.0 * iceR;
          rhs[k0] = T[i][j][0] + 2.0 * dtTempAge * DEFAULT_OCEAN_HEAT_FLUX / (rho_c_I * dz);
          if (!isMarginal) {
            const PetscScalar UpTu = (u[i][j][0] < 0) ?
                u[i][j][0] * (T[i+1][j][0] - T[i][j][0]) / dx :
                u[i][j][0] * (T[i][j][0] - T[i-1][j][0]) / dx;
            const PetscScalar UpTv = (v[i][j][0] < 0) ?
                v[i][j][0] * (T[i][j+1][0] - T[i][j][0]) / dy :
                v[i][j][0] * (T[i][j][0] - T[i][j-1][0]) / dy;
            // for w, always upwind *up* from base
            const PetscScalar UpTw = w[i][j][0] * (T[i][j][1] - T[i][j][0]) / dz;
            rhs[k0] += dtTempAge * (Sigma[i][j][0] - UpTu - UpTv - UpTw) / 2;
          }
        } else { // there is *grounded* ice; ice/bedrock interface
          // {from FV across interface}
          const PetscScalar rho_c_ratio = rho_c_I / rho_c_av;
          rhs[k0] = T[i][j][0] + dtTempAge * (Rb[i][j] / (rho_c_av * dz));
          if (!isMarginal) {
            // for w, always upwind *up* from interface
            const PetscScalar UpTu = (u[i][j][0] < 0) ?
                u[i][j][0] * (T[i+1][j][0] - T[i][j][0]) / dx :
                u[i][j][0] * (T[i][j][0] - T[i-1][j][0]) / dx;
            const PetscScalar UpTv = (v[i][j][0] < 0) ?
                v[i][j][0] * (T[i][j+1][0] - T[i][j][0]) / dy :
                v[i][j][0] * (T[i][j][0] - T[i][j-1][0]) / dy;
            const PetscScalar UpTw = w[i][j][0] * (T[i][j][1] - T[i][j][0]) / dz;
            rhs[k0] += dtTempAge * rho_c_ratio * 0.5 * Sigma[i][j][0];
            rhs[k0] -= dtTempAge * rho_c_ratio
                            * (0.5 * (UpTu + UpTv + UpTw) + T[i][j][0] * w[i][j][0] / dz);
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
        L[k0+k] = -iceR; D[k0+k] = 1+2*iceR; U[k0+k] = -iceR;
        rhs[k0+k] = T[i][j][k];
        if (!isMarginal) {
          const PetscScalar UpTu = (u[i][j][k] < 0) ?
              u[i][j][k] * (T[i+1][j][k] - T[i][j][k]) / dx :
              u[i][j][k] * (T[i][j][k] - T[i-1][j][k]) / dx;
          const PetscScalar UpTv = (v[i][j][k] < 0) ?
              v[i][j][k] * (T[i][j+1][k] - T[i][j][k]) / dy :
              v[i][j][k] * (T[i][j][k] - T[i][j-1][k]) / dy;
          const PetscScalar UpTw = (w[i][j][k] < 0) ?
              w[i][j][k] * (T[i][j][k+1] - T[i][j][k]) / dz :
              w[i][j][k] * (T[i][j][k] - T[i][j][k-1]) / dz;
          rhs[k0+k] += dtTempAge * (Sigma[i][j][k] - UpTu - UpTv - UpTw);
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
        Tb[i][j][k] = x[k];
      }

      // prepare for melting/refreezing
      PetscScalar Hmeltnew = Hmelt[i][j];
      
      // insert solution for generic ice segments
      for (PetscInt k=1; k <= ks; k++) {
//      for (PetscInt k=1; k < ks; k++) {
        if (allowAboveMelting == PETSC_TRUE) {
          Tnew[i][j][k] = x[k0 + k];
        } else {
          const PetscScalar depth = H[i][j] - k * dz;
          const PetscScalar Tpmp = ice.meltingTemp - ice.beta_CC_grad * depth;
          if (x[k0 + k] > Tpmp) {
            Tnew[i][j][k] = Tpmp;
            PetscScalar Texcess = x[k0 + k] - Tpmp; // always positive
            excessToFromBasalMeltLayer(rho_c_I, k * dz, &Texcess, &Hmeltnew);
            // Texcess  will always come back zero here; ignor it
          } else {
            Tnew[i][j][k] = x[k0 + k];
          }
        }
        if (Tnew[i][j][k] < globalMinAllowedTemp) {
           ierr = PetscPrintf(PETSC_COMM_SELF,
              "  [[too low (<200) generic segment temp T = %f at %d,%d,%d; proc %d; mask=%f; w=%f]]\n",
              Tnew[i][j][k],i,j,k,grid.rank,mask[i][j],w[i][j][k]*secpera); CHKERRQ(ierr);
           myLowTempCount++;
        }
      }
      
      // insert solution for ice/rock interface (or base of ice shelf) segment
      if (ks > 0) {
        if (allowAboveMelting == PETSC_TRUE) {
          Tnew[i][j][0] = x[k0];
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
          Tnew[i][j][0] = Tpmp + Texcess;
          if (Tnew[i][j][0] > (Tpmp + 0.00001)) {
            SETERRQ(1,"updated temperature came out above Tpmp");
          }
        }
        if (Tnew[i][j][0] < globalMinAllowedTemp) {
           ierr = PetscPrintf(PETSC_COMM_SELF,
              "  [[too low (<200) ice/rock segment temp T = %f at %d,%d; proc %d; mask=%f; w=%f]]\n",
              Tnew[i][j][0],i,j,grid.rank,mask[i][j],w[i][j][0]*secpera); CHKERRQ(ierr);
           myLowTempCount++;
        }
      } else {
        Hmeltnew = 0.0;
      }
      
      // we must agree on redundant values T(z=0) at top of bedrock and at bottom of ice
      if (ks > 0) {
        Tb[i][j][k0] = Tnew[i][j][0];
      } else {
        // if floating then top of bedrock sees ocean
        const PetscScalar floating_base = - (ice.rho/ocean.rho) * H[i][j];
        if (b[i][j] < floating_base - 1.0) {
          Tb[i][j][k0] = ice.meltingTemp;
        } else { // top of bedrock sees atmosphere
          Tb[i][j][k0] = Ts[i][j];
        }
      }
      // check bedrock solution        
      for (PetscInt k=0; k <= k0; k++) {
        if (Tb[i][j][k] < globalMinAllowedTemp) {
           ierr = PetscPrintf(PETSC_COMM_SELF,
              "  [[too low (<200) bedrock temp T = %f at %d,%d,%d; proc %d; mask=%f]]\n",
              Tb[i][j][k],i,j,k,grid.rank,mask[i][j]); CHKERRQ(ierr);
           myLowTempCount++;
        }
      }
      
      // set to air temp above ice
      for (PetscInt k=ks; k<Mz; k++) {
        Tnew[i][j][k] = Ts[i][j];
      }

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
  // but (11/16/06) vHmelt, vbasalMeltRate and vTb will never need to communicate ghosted values
  // (i.e. horizontal stencil neighbors);  vT is communicated by temperatureAgeStep()

  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vGhf, &Ghf); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vRb, &Rb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vTnew, &Tnew); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  delete [] Lp; delete [] D; delete [] U; delete [] x; delete [] rhs; delete [] work;
  return 0;
}


//! Compute the melt water which should go to the base if \f$T\f$ is above pressure-melting.
PetscErrorCode IceModel::excessToFromBasalMeltLayer(
                PetscScalar rho_c, PetscScalar z,
                PetscScalar *Texcess, PetscScalar *Hmelt) {

  const PetscScalar darea = grid.p->dx * grid.p->dy;
  const PetscScalar dvol = darea * grid.p->dz;  // only for equally-spaced!
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
        const PetscScalar dTemp = ice.latentHeat * ice.rho * (*Hmelt) / (rho_c * grid.p->dz);
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
The age equation is\f$ \frac{d\tau}{dt} = 1\f$, that is,
    \f[ \frac{\partial \tau}{\partial t} + u \frac{\partial \tau}{\partial x} + v \frac{\partial \tau}{\partial y} + w \frac{\partial \tau}{\partial z} = 1\f]
where \f$\tau(t,x,y,z)\f$ is the age of the ice and \f$(u,v,w)\f$  is the three dimensional
velocity field.  This equation is hyperbolic (purely advective).  
The boundary condition is that when the ice fell as snow it had age zero.  
That is, \f$\tau(t,x,y,h(t,x,y)) = 0\f$ in accumulation areas, while there is no 
boundary condition elsewhere (as the characteristics go outward elsewhere).  At this point 
the refreeze case, either grounded basal ice or marine basal ice, is not handled correctly.
@endcond

@cond NUMERIC
The numerical method is first-order upwind.
@endcond
 */
PetscErrorCode IceModel::ageStep(PetscScalar* CFLviol) {
  // update age field vtaunew
  PetscErrorCode  ierr;

  const PetscInt      Mz = grid.p->Mz;
  const PetscScalar   dx = grid.p->dx, dy = grid.p->dy, dz = grid.p->dz;
  const PetscScalar   cflx = dx/dtTempAge, cfly = dy/dtTempAge, cflz = dz/dtTempAge;

  Vec     vtaunew = vWork3d[1];  // will be communicated by temperatureAgeStep()
  PetscScalar ***u, ***v, ***w, ***tau, ***taunew, **H;
  
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vtau, &tau); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vtaunew, &taunew); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/dz));
      if (ks >= Mz) {
        ierr = PetscPrintf(grid.com,
                 "[[error LOCATION: i, j, ks, H = %5d %5d %5d %10.2f]]\n",i, j, ks, H[i][j]); 
        SETERRQ(1,"Vertical grid exceeded");
      }
      
      // only effects of this is whether vertical velocities are used in advection
      const bool isMarginal = checkThinNeigh(H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                             H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);

      for (PetscInt k=0; k<ks; k++) {
        // age evolution is pure advection (so provides check on temp calculation)
        // check CFL conditions at each point, then upwind for age
        if (PetscAbs(u[i][j][k]) > cflx)  *CFLviol += 1.0;
        if (PetscAbs(v[i][j][k]) > cfly)  *CFLviol += 1.0;
        PetscScalar     rtau;
        rtau = (u[i][j][k] < 0) ?
          u[i][j][k] * (tau[i+1][j][k] - tau[i][j][k]) / dx :
          u[i][j][k] * (tau[i][j][k] - tau[i-1][j][k]) / dx;
        rtau += (v[i][j][k] < 0) ?
          v[i][j][k] * (tau[i][j+1][k] - tau[i][j][k]) / dy :
          v[i][j][k] * (tau[i][j][k] - tau[i][j-1][k]) / dy;
        if (!isMarginal) {
          if (PetscAbs(w[i][j][k]) > cflz)  *CFLviol += 1.0;
          // if w upward at k=0 then ignor contribution to age
          if ((k > 0) || (w[i][j][k] < 0)) {
            rtau += (w[i][j][k] < 0) ?    // note: this is lowest-order upwinding
              w[i][j][k] * (tau[i][j][k+1] - tau[i][j][k]) / dz :
              w[i][j][k] * (tau[i][j][k] - tau[i][j][k-1]) / dz;
          }
        }
        taunew[i][j][k] = tau[i][j][k] + dtTempAge * (1.0 - rtau);
      }      
      for (PetscInt k=ks; k<Mz; k++) {
        taunew[i][j][k] = 0.0;  // age of ice above (and at) surface is zero years
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vtau, &tau); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vtaunew, &taunew); CHKERRQ(ierr);
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

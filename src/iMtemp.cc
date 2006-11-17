// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
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

#include <petscda.h>
#include <petscksp.h>
#include "iceModel.hh"


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


PetscErrorCode IceModel::temperatureStep() {
  // update temp fields vTnew, vTb

  // This procedure involves many choices.  See the "On the temperature problem in
  //    a column of flowing, sliding ice over bedrock" section in the model notes 
  //    eqns3D.tex.
  // Note that we work from the bottom of the column upward in building the system 
  //    to solve for the semi-implicit time-stepping scheme.
  
  PetscErrorCode  ierr;

  const PetscInt      Mbz = grid.p->Mbz, Mz = grid.p->Mz;
  const PetscScalar   dx = grid.p->dx, dy = grid.p->dy, dz = grid.p->dz;
  const PetscScalar   dvol = dx * dy * dz;  // only for equally-spaced!
  const PetscScalar   iceK = ice.k / (ice.rho * ice.c_p);
  const PetscScalar   iceR = iceK * dtTempAge / PetscSqr(dz);
  const PetscScalar   brK = bedrock.k / (bedrock.rho * bedrock.c_p);
  const PetscScalar   brR = brK * dtTempAge / PetscSqr(dz);
  const PetscScalar   rho_c_I = ice.rho * ice.c_p;
  const PetscScalar   rho_c_av = (rho_c_I + bedrock.rho * bedrock.c_p) / 2.0;

  Vec     vTnew = vWork3d[0];  // will be communicated by temperatureAgeStep()
  PetscScalar **Ts, **H, **h, **bed, **Ghf, **mask, **Hmelt, **Rb, **basalMeltRate;
  PetscScalar ***T, ***Tb, ***Tnew, ***u, ***v, ***w, ***Sigma;

  PetscScalar *Lp, *L, *D, *U, *x, *rhs, *work;  
  Lp = new PetscScalar[Mz+Mbz-1]; L = Lp-1; // pointer arithmetic
  D = new PetscScalar[Mz+Mbz];
  U = new PetscScalar[Mz+Mbz-1];
  x = new PetscScalar[Mz+Mbz];
  rhs = new PetscScalar[Mz+Mbz];
  work = new PetscScalar[Mz+Mbz];

  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
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
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vTnew, &Tnew); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/dz));
      // go ahead and assume ks < Mz; see ageStep() for check
      
      if (Mbz > 0) { // bedrock present: do k=0:Mbz-1
        // basal condition (at bottom of bedrock when temp is modelled in bedrock)
        // {from FV (finite volume) across bedrock base z_0}
        D[0] = (1 + brR);  U[0] = -brR;  // note L[0] not an allocated location
        rhs[0] = Tb[i][j][0] + (dtTempAge * brK * Ghf[i][j] / dz);
      
        // bedrock only: pure vertical conduction problem
        // {from generic bedrock FV}
        for (PetscInt k=1; k < Mbz; k++) {
          L[k] = -brR; D[k] = 1+2*brR; U[k] = -brR;
          rhs[k] = Tb[i][j][k];
        }
      }
      
      // bottom part of ice (and top of bedrock in some cases): k=Mbz
      if (ks == 0) {
        if (Mbz > 0) { L[Mbz] = 0.0; } // note L[0] not allocated 
        D[Mbz] = 1.0; U[Mbz] = 0.0;
        if (modMask(mask[i][j]) == MASK_FLOATING) {
          rhs[Mbz] = ocean.homol_temp;
        } else {
          rhs[Mbz] = Ts[i][j]; 
        } 
      } else { // ks > 0; there is ice
        if (modMask(mask[i][j]) == MASK_FLOATING) {
          // base of ice shelf: set T[Mbz] to 2 deg above T_pmp
          if (Mbz > 0) { L[Mbz] = 0.0; }
          D[Mbz] = 1.0; U[Mbz] = 0.0;
          rhs[Mbz] = ocean.homol_temp - ice.beta_CC_grad * H[i][j];
        } else { // there is *grounded* ice
          // ice/bedrock interface
          // {from FV across interface}
          const PetscScalar rho_c_ratio = rho_c_I / rho_c_av;
          const PetscScalar UpTu = (u[i][j][0] < 0) ?
              u[i][j][0] * (T[i+1][j][0] - T[i][j][0]) / dx :
              u[i][j][0] * (T[i][j][0] - T[i-1][j][0]) / dx;
          const PetscScalar UpTv = (v[i][j][0] < 0) ?
              v[i][j][0] * (T[i][j+1][0] - T[i][j][0]) / dy :
              v[i][j][0] * (T[i][j][0] - T[i][j-1][0]) / dy;
          // for w, always upwind *up* from interface
          const PetscScalar UpTw = w[i][j][0] * (T[i][j][1] - T[i][j][0]) / dz;
          const PetscScalar iceReff = ice.k * dtTempAge / (rho_c_av * dz * dz);
          const PetscScalar brReff = bedrock.k * dtTempAge / (rho_c_av * dz * dz);
          rhs[Mbz] = T[i][j][0] 
                     + dtTempAge * ( rho_c_ratio * 0.5 * Sigma[i][j][0] 
                                     + Rb[i][j] / (rho_c_av * dz) 
                                     - rho_c_ratio * 0.5 * (UpTu + UpTv + UpTw)
                                     - rho_c_ratio * T[i][j][0] * w[i][j][0] / dz );
          if (Mbz > 0) { // there is bedrock
            L[Mbz] = - brReff; D[Mbz] = 1 + iceReff + brReff; U[Mbz] = - iceReff;
          } else { // no bedrock
            // L[Mbz] = 0.0;  (note this is not an allocated location!) 
            D[Mbz] = 1 + iceReff; U[Mbz] = - iceReff;
            rhs[Mbz] += dtTempAge * Ghf[i][j] / (rho_c_av * dz);
          }
        }
      }

      // generic ice segment
      for (PetscInt k=1; k<ks; k++) {
        const PetscScalar UpTu = (u[i][j][k] < 0) ?
            u[i][j][k] * (T[i+1][j][k] - T[i][j][k]) / dx :
            u[i][j][k] * (T[i][j][k] - T[i-1][j][k]) / dx;
        const PetscScalar UpTv = (v[i][j][k] < 0) ?
            v[i][j][k] * (T[i][j+1][k] - T[i][j][k]) / dy :
            v[i][j][k] * (T[i][j][k] - T[i][j-1][k]) / dy;
        const PetscScalar UpTw = (w[i][j][k] < 0) ?
            w[i][j][k] * (T[i][j][k+1] - T[i][j][k]) / dz :
            w[i][j][k] * (T[i][j][k] - T[i][j][k-1]) / dz;
        L[Mbz+k] = -iceR; D[Mbz+k] = 1+2*iceR; U[Mbz+k] = -iceR;
        rhs[Mbz+k] = T[i][j][k] + dtTempAge * (Sigma[i][j][k] - UpTu - UpTv - UpTw);
      }
      
      // surface b.c.
      if (Mbz+ks>0) {
        rhs[Mbz+ks-1] += iceR * Ts[i][j]; // L[Mbz+ks-1], set above, will be ignored in this case
      }

      // solve system; melting not addressed yet
      if (Mbz+ks>0) {
        ierr = solveTridiagonalSystem(L, D, U, x, rhs, work, Mbz+ks);
        if (ierr != 0) {
          SETERRQ3(1, "Tridiagonal solve failed at (%d,%d) with zero pivot in position %d.",
               i, j, ierr);
        }
      }

      // insert bedrock solution        
      for (PetscInt k=0; k < Mbz; k++) {
        Tb[i][j][k] = x[k];
      }

      // prepare for melting/refreezing
      PetscScalar Hmeltnew = Hmelt[i][j];
      
      // insert ice solution for generic ice segments
      for (PetscInt k=1; k < ks; k++) {
        if (allowAboveMelting == PETSC_TRUE) {
          Tnew[i][j][k] = x[Mbz+k];
        } else {
          const PetscScalar depth = H[i][j] - k * dz;
          const PetscScalar Tpmp = ice.meltingTemp - ice.beta_CC_grad * depth;
          if (x[Mbz+k] <= Tpmp) {
            Tnew[i][j][k] = x[Mbz + k];
          } else { // T is above pressure-melting temp, so excess temperature
                   // needs to be turned into melt water at base
            const PetscScalar dE = ice.rho * ice.c_p * (x[Mbz+k] - Tpmp) * dvol;
            const PetscScalar massmelted = dE / ice.latentHeat;
            Hmeltnew += massmelted / (ice.rho * dvol);
            Tnew[i][j][k] = Tpmp;
          }
        }
      }
      
      // insert solution for ice/rock interface (or base of ice shelf) segment
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        Tnew[i][j][0] = x[Mbz]; // see above for equation which determines this
      } else {
        const PetscScalar Tpmp = ice.meltingTemp - ice.beta_CC_grad * H[i][j];
        if (Hmeltnew > 0.0) { // have melt water; conservation of energy determines Hmelt; T=Tpmp
          const PetscScalar dEtemp = rho_c_av * (x[Mbz] - Tpmp) * dvol; // could be negative!
          const PetscScalar massmelted = dEtemp / ice.latentHeat; // could be negative
          Hmeltnew += massmelted / (ice.rho * dvol);
          if (Hmeltnew <= 0.0) {
            Hmeltnew = 0.0;
            Tnew[i][j][Mbz] = Tpmp - 1e-3;
          } else {
            Tnew[i][j][Mbz] = Tpmp;
          }
        } else {  // Hmelt=0.0 from previous step, and no melt from above
          if (allowAboveMelting == PETSC_TRUE) {
            Tnew[i][j][0] = x[Mbz];
          } else {
            if (x[Mbz] <= Tpmp) {
              Tnew[i][j][0] = x[Mbz];
            } else { // T is above pressure-melting temp, so excess temperature
                     // needs to be turned into melt water at base
              const PetscScalar dE = rho_c_av * (x[Mbz] - Tpmp) * dvol;
              const PetscScalar massmelted = dE / ice.latentHeat;
              Hmeltnew += massmelted / (ice.rho * dvol);
              Tnew[i][j][Mbz] = Tpmp;
            }
          }
        }
      }
       
      // set to air temp above ice
      for (PetscInt k=ks; k<Mz; k++) {
        Tnew[i][j][k] = Ts[i][j];
      }

      // basaMeltRate[][] is rate of change of Hmelt[][]; thus it can be negative
      basalMeltRate[i][j] = (Hmeltnew - Hmelt[i][j]) / dtTempAge;
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        Hmelt[i][j] = 0.0;
      } else {
        Hmelt[i][j] = PetscMin(1.0,Hmeltnew);  // Hmelt at most 1.0 m
      }
    } 
  }
  
  // note that in above code 4 scalar fields were modified: vHmelt, vbasalMeltRate, vTb, and vT
  // but (11/16/06) vHmelt, vbasalMeltRate and vTb will never need to communicate ghosted values
  // (i.e. horizontal stencil neighbors);  vT is communicated by temperatureAgeStep()

  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
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
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vTnew, &Tnew); CHKERRQ(ierr);
  delete [] Lp; delete [] D; delete [] U; delete [] x; delete [] rhs; delete [] work;
  return 0;
}


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
      for (PetscInt k=0; k<ks; k++) {
        // age evolution is pure advection (so provides check on temp calculation)
        // check CFL conditions at each point, then upwind for age
        if (PetscAbs(u[i][j][k]) > cflx)  *CFLviol += 1.0;
        if (PetscAbs(v[i][j][k]) > cfly)  *CFLviol += 1.0;
        if (PetscAbs(w[i][j][k]) > cflz)  *CFLviol += 1.0;
        PetscScalar     rtau;
        rtau = (u[i][j][k] < 0) ?
          u[i][j][k] * (tau[i+1][j][k] - tau[i][j][k]) / dx :
          u[i][j][k] * (tau[i][j][k] - tau[i-1][j][k]) / dx;
        rtau += (v[i][j][k] < 0) ?
          v[i][j][k] * (tau[i][j+1][k] - tau[i][j][k]) / dy :
          v[i][j][k] * (tau[i][j][k] - tau[i][j-1][k]) / dy;
        // if w upward at k=0 then ignor contribution to age
        if ((k > 0) || (w[i][j][k] < 0)) {
          rtau += (w[i][j][k] < 0) ?    // note: this is lowest-order upwinding
            w[i][j][k] * (tau[i][j][k+1] - tau[i][j][k]) / dz :
            w[i][j][k] * (tau[i][j][k] - tau[i][j][k-1]) / dz;
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


PetscErrorCode IceModel::solveTridiagonalSystem(
         const PetscScalar* L, const PetscScalar* D, const PetscScalar* U,
         PetscScalar* x, const PetscScalar* r, PetscScalar* a, const int n) const {

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

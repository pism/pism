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

  PetscScalar  myCFLviolcount = 0.0;  // it is a count but it is "PetscScalar"
                                      // because that type works with PetscGlobalSum()
  // do CFL and vertical grid blow-out checking in ageStep() but not temperatureStep()
  ierr = ageStep(&myCFLviolcount); CHKERRQ(ierr);
  ierr = temperatureStep(); CHKERRQ(ierr);

  // no communication done in ageStep(), temperatureStep(); all done here:
  ierr = DALocalToLocalBegin(grid.da3, vWork3d[0], INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vWork3d[0], INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vWork3d[1], INSERT_VALUES, vtau); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vWork3d[1], INSERT_VALUES, vtau); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&myCFLviolcount, &CFLviolcount, grid.com); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::temperatureStep() {
  // update temp field
  PetscErrorCode  ierr;

  const PetscInt      Mbz = grid.p->Mbz, Mz = grid.p->Mz;
  const PetscScalar   dx = grid.p->dx, dy = grid.p->dy, dz = grid.p->dz;
  const PetscScalar   iceK = ice.k / (ice.rho * ice.c_p);
  const PetscScalar   iceR = iceK * dtTempAge / PetscSqr(dz);
  const PetscScalar   brK = bedrock.k / (bedrock.rho * bedrock.c_p);
  const PetscScalar   brR = brK * dtTempAge / PetscSqr(dz);

  Vec     vTnew = vWork3d[0];  // will be communicated by temperatureAgeStep()
  PetscScalar **Ts, **H, **h, **bed, **Ghf, **mask, **basalMeltRate;
  PetscScalar ***T, ***Tb, ***Tnew, ***u, ***v, ***w, ***Sigma;

  PetscScalar *Lp, *L, *D, *U, *x, *rhs, *work;  
  Lp = new PetscScalar[Mz+Mbz-1]; L = Lp-1;
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
      if (ks >= Mz) {
        ierr = PetscPrintf(grid.com,
                 "LOCATION: i, j, ks, H = %5d %5d %5d %10.2f\n",i, j, ks, H[i][j]); 
        SETERRQ(1,"Vertical grid exceeded");
      }

      // Basal condition (at bottom of bedrock when temp is modelled in bedrock)
      D[0] = (Mbz > 0) ? bedrock.k/dz : ice.k/dz;
      U[0] = -D[0];
      rhs[0] = Ghf[i][j];  // note this value does not affect floating ice.
                           // (Temps in bedrock *are* modelled when bed is
                           // covered by ocean not ice, but they are not coupled to ice.)
      
      if (Mbz > 0) {
        // Bedrock only: pure vertical conduction problem
        for (PetscInt k=1; k < Mbz; k++) {
          L[k] = -brR; D[k] = 1+2*brR; U[k] = -brR;
          rhs[k] = Tb[i][j][k];
        }

        // Bedrock:Ice interface
        if (ks > 0) { // i.e. if there is ice
          L[Mbz] = -brR; D[Mbz] = 1+brR+iceR; U[Mbz] = -iceR;
          PetscScalar r = (u[i][j][0] < 0) ?
            u[i][j][0] * (T[i+1][j][0] - T[i][j][0]) / dx :
            u[i][j][0] * (T[i][j][0] - T[i-1][j][0]) / dx;
          r += (v[i][j][0] < 0) ?
            v[i][j][0] * (T[i][j+1][0] - T[i][j][0]) / dy :
            v[i][j][0] * (T[i][j][0] - T[i][j-1][0]) / dy;
          r += w[i][j][0] * (T[i][j][0+1] - T[i][j][0]) / dz;  // always difference up
          rhs[Mbz] = T[i][j][0] + dtTempAge * (Sigma[i][j][0] - r);
          if (modMask(mask[i][j]) != MASK_FLOATING) {
            const PetscScalar ub = u[i][j][0];
            const PetscScalar vb = v[i][j][0];
            PetscScalar basal_stress_x, basal_stress_y;
            if ((modMask(mask[i][j]) == MASK_DRAGGING) && (useMacayealVelocity)) {
              basal_stress_x = basalDrag(ub, vb) * ub;
              basal_stress_y = basalDrag(ub, vb) * vb;
            } else { // usual SIA assumption: basal driving shear stress
              basal_stress_x = -ice.rho * ice.grav * H[i][j]
                                        * (h[i+1][j] - h[i-1][j]) / (2 * dx);
              basal_stress_y = -ice.rho * ice.grav * H[i][j]
                                        * (h[i][j+1] - h[i][j-1]) / (2 * dy);
            }
            const PetscScalar basal_heating = basal_stress_x * ub + basal_stress_y * vb;
//            if (basal_heating < 0.0) {
//              SETERRQ(1,"basal heating negative in IceModel::temperatureStep()");
//            }
            const PetscScalar avrhocp = (ice.rho*ice.c_p + bedrock.rho*bedrock.c_p) / 2.0;
            rhs[Mbz] += dtTempAge * (basal_heating / (avrhocp * dz));
          } else {  // If floating, don't add basal strain heating but do
            // set temp at top of bedrock.  (It is in contact with ocean.)
            // Set to *pressure-melting temp at bottom of ice shelf*.
            // Note that in this case Sigma[i][j][0] will have ocean heating added, too.
            D[Mbz-1]=1.0; L[Mbz-1]=0.0; U[Mbz-1]=0.0;
            rhs[Mbz - 1] = 273.15 - ice.beta_CC_grad * H[i][j];  
          }
        } else {  // if no ice or essentially no ice set bottom of ice to surface temp
          L[Mbz] = 0.0; D[Mbz] = 1.0; U[Mbz] = 0.0;
          rhs[Mbz] = Ts[i][j];  
          if (modMask(mask[i][j]) == MASK_FLOATING) { // if no ice, or essentially no ice, 
                                                      // and if floating, set top of bedrock
                                                      // to triple point
            D[Mbz-1]=1.0; L[Mbz-1]=0.0; U[Mbz-1]=0.0;
            rhs[Mbz - 1] = ice.meltingTemp;  
          }
        }
      }

      if ((ks < 2) && (Mbz == 0)) { // FIXME: if no bedrock, this is not right when ice is floating!
        // if zero or one point above bed, assume T_zz=0, T(z=H)=Ts, T_z(z=0)=-G/k
        const PetscScalar GoK = Ghf[i][j]/ice.k;
        x[0] = Ts[i][j] + GoK * H[i][j];
        if (ks > 0)   x[1] = T[i][j][0] - GoK * dz; 
      } else { // ... and solve tridiagonal system for temps
        // Ice only: diffusion plus advection and strain heating
        for (PetscInt k=1; k<ks; k++) {
          L[Mbz+k] = -iceR; D[Mbz+k] = 1+2*iceR; U[Mbz+k] = -iceR;
          // upwind for temperature:
          PetscScalar r = (u[i][j][k] < 0) ?
            u[i][j][k] * (T[i+1][j][k] - T[i][j][k]) / dx :
            u[i][j][k] * (T[i][j][k] - T[i-1][j][k]) / dx;
          r += (v[i][j][k] < 0) ?
            v[i][j][k] * (T[i][j+1][k] - T[i][j][k]) / dy :
            v[i][j][k] * (T[i][j][k] - T[i][j-1][k]) / dy;
          r += (w[i][j][k] < 0) ?         // note: this is only standard upwinding
            w[i][j][k] * (T[i][j][k+1] - T[i][j][k]) / dz :
            w[i][j][k] * (T[i][j][k] - T[i][j][k-1]) / dz;
          rhs[Mbz+k] = T[i][j][k] + dtTempAge * (Sigma[i][j][k] - r);
        }
        // Surface condition and solve
        if (Mbz+ks>0) {
          rhs[Mbz+ks-1] += iceR * Ts[i][j];
          ierr = solveTridiagonalSystem(L, D, U, x, rhs, work, Mbz+ks);
          if (ierr != 0) {
            SETERRQ3(1, "Tridiagonal solve failed at (%d,%d) with zero pivot in position %d.",
                   i, j, ierr);
          }
        }
      }

      for (PetscInt k=0; k < Mbz; k++) {
        Tb[i][j][k] = PetscMax(x[k], globalMinTemp);
      }

      PetscScalar meltingEnergyFlux = 0; // J / m^2
      if (allowAboveMelting == PETSC_TRUE) {
        for (PetscInt k=0; k < ks; k++) {
          Tnew[i][j][k] = PetscMax(x[Mbz+k], globalMinTemp);
        }
      } else {
        for (PetscInt k=0; k < ks; k++) {
          // tchange is temperature difference between pressure melting at depth z
          //     and ice.meltingTemp
          const PetscScalar tchange = ice.beta_CC_grad * (H[i][j] - k*dz);
          PetscScalar excess = x[Mbz + k] - (ice.meltingTemp - tchange);
          if (excess > 0) {
            Tnew[i][j][k] = x[Mbz + k] - excess;
            // Use the trapezoid rule
            meltingEnergyFlux += (excess * ((k > 0) ? 1 : 0.5)
                                  * dz * ice.rho * ice.c_p); // J / m^2
          } else {
            Tnew[i][j][k] = x[Mbz + k];
          }
          // Tnew[i][j][k] = (x[Mbz+k] + tchange < ice.meltingTemp)
          //   ? x[Mbz+k] : ice.meltingTemp - tchange;
          Tnew[i][j][k] = PetscMax(Tnew[i][j][k], globalMinTemp);
        }
      }
      basalMeltRate[i][j] = meltingEnergyFlux / (ice.rho * ice.latentHeat * dtTempAge); // m / s

      for (PetscInt k=ks; k<Mz; k++) {
        Tnew[i][j][k] = PetscMax(Ts[i][j], globalMinTemp);
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vGhf, &Ghf); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vTnew, &Tnew); CHKERRQ(ierr);

  // note that in above code 3 scalar fields were modified: vbasalMeltRate, vTb, and vT
  // but (8/03/06) vbasalMeltRate and vTb will never need to communicate ghosted values
  // (i.e. horizontal stencil neighbors);  vT is communicated by temperatureAgeStep()

  delete [] Lp; delete [] D; delete [] U; delete [] x; delete [] rhs; delete [] work;

  return 0;
}


PetscErrorCode IceModel::ageStep(PetscScalar* CFLviol) {
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

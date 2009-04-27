// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petsc.h>
#include "pism_const.hh"   // e.g. MASK_FLOATING and PismModMask()
#include "iceModelVec.hh"
#include "columnSystem.hh"


columnSystemCtx::columnSystemCtx(int my_nmax) : nmax(my_nmax) {
  if (nmax < 1) {
    PetscPrintf(PETSC_COMM_WORLD,
      "columnSystem ERROR: nmax of system too small\n");
    PetscEnd();
  }
  if (nmax > 1000000) {
    PetscPrintf(PETSC_COMM_WORLD,
      "columnSystem ERROR: nmax of system unreasonable (> 10^6)\n");
    PetscEnd();
  }
  Lp   = new PetscScalar[nmax-1];
  L    = Lp-1; // ptr arith.; note L[0]=Lp[-1] not allocated
  D    = new PetscScalar[nmax];
  U    = new PetscScalar[nmax-1];
  rhs  = new PetscScalar[nmax];
  work = new PetscScalar[nmax];
}


columnSystemCtx::~columnSystemCtx() {
  delete [] Lp;
  delete [] D;
  delete [] U;
  delete [] rhs;
  delete [] work;
}


//! The actual code for solving a tridiagonal system.  Return code has diagnostic importance.
/*!
This is modified slightly from a Numerical Recipes version.

Input size n is size of instance.  Requires n <= columnSystemCtx::nmax.

Solution of system in x.

Success is return code zero.  Positive return code gives location of zero pivot.
Negative return code indicates a software problem.
 */
PetscErrorCode columnSystemCtx::solveTridiagonalSystem(
                  PetscInt n, PetscScalar **x) {
  if (x == NULL) { SETERRQ(-999,"x is NULL in columnSystemCtx"); }
  if (*x == NULL) { SETERRQ(-998,"*x is NULL in columnSystemCtx"); }
  if (n < 1) { SETERRQ(-997,"instance size n < 1 in columnSystemCtx"); }
  if (n > nmax) { SETERRQ(-996,"instance size n too large in columnSystemCtx"); }
  PetscScalar b;
  b = D[0];
  if (b == 0.0) { return 1; }
  (*x)[0] = rhs[0]/b;
  for (int i=1; i<n; ++i) {
    work[i] = U[i-1]/b;
    b = D[i] - L[i] * work[i];
    if (b == 0.0) { return i+1; }
    (*x)[i] = (rhs[i] - L[i] * (*x)[i-1]) / b;
  }
  for (int i=n-2; i>=0; --i) {
    (*x)[i] -= work[i+1] * (*x)[i+1];
  }
  return 0;
}


ageSystemCtx::ageSystemCtx(int my_Mz)
      : columnSystemCtx(my_Mz) { // size of system is Mz
  callcount = 0;
  dx = -1.0;
  dy = -1.0;
  dtTempAge = -1.0;
  dzEQ = -1.0;
  nuEQ = -1.0;
}


PetscErrorCode ageSystemCtx::ageSetConstants(
    PetscScalar my_dx, PetscScalar my_dy, PetscScalar my_dtTempAge,
    PetscScalar my_dzEQ) {

  if (callcount > 1) {
    SETERRQ(1,"ageSetConstants() should only be called once");
  }
  callcount++;

  // is all this checking necessary?
  if (my_dx <= 0.0) { SETERRQ(2,"invalid dx in ageSetConstants()"); }
  dx = my_dx;
  if (my_dy <= 0.0) { SETERRQ(3,"invalid dy in ageSetConstants()"); }
  dy = my_dy;
  if (my_dtTempAge <= 0.0) { SETERRQ(4,"invalid dtTempAge in ageSetConstants()"); }
  dtTempAge = my_dtTempAge;
  if (my_dzEQ <= 0.0) { SETERRQ(6,"invalid dzEQ in ageSetConstants()"); }
  dzEQ = my_dzEQ;
  
  nuEQ = dtTempAge / dzEQ;
  
  return 0;
}


PetscErrorCode ageSystemCtx::ageColumnSetUpAndSolve(
    PetscInt i, PetscInt j, PetscInt ks,
    PetscScalar *u, PetscScalar *v, PetscScalar *w,
    IceModelVec3 &tau3,
    PetscScalar **x) {

  PetscErrorCode ierr;
  if (callcount == 0) {
    SETERRQ(1,"ageColumnSetUpAndSolve() says ageSetConstants() has not been called");
  }

  // set up system: 0 <= k < ks
  for (PetscInt k=0; k<ks; k++) {
    planeStar ss;  // note ss.ij = tau[k]
    ierr = tau3.getPlaneStarZ(i,j,k * dzEQ,&ss); CHKERRQ(ierr);
    // do lowest-order upwinding, explicitly for horizontal
    rhs[k] =  (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx
                         : u[k] * (ss.ij  - ss.im1) / dx;
    rhs[k] += (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy
                         : v[k] * (ss.ij  - ss.jm1) / dy;
    // note it is the age eqn: dage/dt = 1.0 and we have moved the hor.
    //   advection terms over to right:
    rhs[k] = ss.ij + dtTempAge * (1.0 - rhs[k]);

    // do lowest-order upwinding, *implicitly* for vertical
    PetscScalar AA = nuEQ * w[k];
    if (k > 0) {
      if (AA >= 0) { // upward velocity
        L[k] = - AA;
        D[k] = 1.0 + AA;
        U[k] = 0.0;
      } else { // downward velocity; note  -AA >= 0
        L[k] = 0.0;
        D[k] = 1.0 - AA;
        U[k] = + AA;
      }
    } else { // k == 0 case
      // note L[0] not an allocated location
      if (AA > 0) { // if strictly upward velocity apply boundary condition:
                    // age = 0 because ice is being added to base
        D[0] = 1.0;
        U[0] = 0.0;
        rhs[0] = 0.0;
      } else { // downward velocity; note  -AA >= 0
        D[0] = 1.0 - AA;
        U[0] = + AA;
        // keep rhs[0] as is
      }
    }
  }  // done "set up system: 0 <= k < ks"
      
  // surface b.c. at ks
  if (ks>0) {
    L[ks] = 0;
    D[ks] = 1.0;   // ignore U[ks]
    rhs[ks] = 0.0;  // age zero at surface
  }

  // solve it
  return solveTridiagonalSystem(ks+1,x);
}


tempSystemCtx::tempSystemCtx(int my_Mz, int my_Mbz)
      : columnSystemCtx(my_Mz + my_Mbz - 1) {

  Mz = my_Mz;
  Mbz = my_Mbz;
  k0 = Mbz - 1; // max size nmax of system is Mz + k0 = Mz + Mbz - 1

  callcount = 0;

  dx = -1;
  dy = -1;
  dtTempAge = -1;
  dzEQ = -1;
  dzbEQ = -1;
  ice_rho = -1;
  ice_c_p = -1;
  ice_k = -1;
  bed_thermal_rho = -1;
  bed_thermal_c_p = -1;
  bed_thermal_k = -1;
}


PetscErrorCode tempSystemCtx::tempSetConstants(
    PetscScalar my_dx, PetscScalar my_dy, PetscScalar my_dtTempAge,
    PetscScalar my_dzEQ, PetscScalar my_dzbEQ,
    PetscScalar my_ice_rho, PetscScalar my_ice_c_p, PetscScalar my_ice_k,
    PetscScalar my_bed_thermal_rho, PetscScalar my_bed_thermal_c_p,
       PetscScalar my_bed_thermal_k) {

  if (callcount > 1) {
    SETERRQ(1,"ageSetConstants() should only be called once");
  }
  callcount++;

  dx = my_dx;
  dy = my_dy;
  dtTempAge = my_dtTempAge;
  dzEQ = my_dzEQ;
  dzbEQ = my_dzbEQ;
  ice_rho = my_ice_rho;
  ice_c_p = my_ice_c_p;
  ice_k = my_ice_k;
  bed_thermal_rho = my_bed_thermal_rho;
  bed_thermal_c_p = my_bed_thermal_c_p;
  bed_thermal_k = my_bed_thermal_k;
  
  nuEQ = dtTempAge / dzEQ;
  rho_c_I = ice_rho * ice_c_p;
  rho_c_br = bed_thermal_rho * bed_thermal_c_p;
  rho_c_av = (dzEQ * rho_c_I + dzbEQ * rho_c_br) / (dzEQ + dzbEQ);
  iceK = ice_k / rho_c_I;
  iceR = iceK * dtTempAge / PetscSqr(dzEQ);
  brK = bed_thermal_k / rho_c_br;
  brR = brK * dtTempAge / PetscSqr(dzbEQ);
  rho_c_ratio = rho_c_I / rho_c_av;
  dzav = 0.5 * (dzEQ + dzbEQ);
  iceReff = ice_k * dtTempAge / (rho_c_av * dzEQ * dzEQ);
  brReff = bed_thermal_k * dtTempAge / (rho_c_av * dzbEQ * dzbEQ);

  return 0;
}


PetscErrorCode tempSystemCtx::tempColumnSetUpAndSolve(
    PetscInt i, PetscInt j,
    PetscInt ks, bool isMarginal, PetscScalar lambda,
    PetscScalar *T, PetscScalar *Tb,
    PetscScalar *u, PetscScalar *v, PetscScalar *w, PetscScalar *Sigma,
    PetscScalar Ghf_ij, PetscScalar Ts_ij, PetscScalar mask_ij,
    PetscScalar Tshelfbase_ij, PetscScalar Rb_ij,
    IceModelVec3 &T3,
    PetscScalar **x) {

  PetscErrorCode ierr;
  if (callcount == 0) {
    SETERRQ(1,"tempColumnSetUpAndSolve() says tempSetConstants() has not been called");
  }

  if (Mbz > 1) { // bedrock present: build k=0:Mbz-2 eqns
    // gives O(\Delta t,\Delta z^2) convergence in Test K for equal spaced grid;
    // note L[0] not an allocated location:
    D[0] = (1.0 + 2.0 * brR);
    U[0] = - 2.0 * brR;  
    rhs[0] = Tb[0] + 2.0 * dtTempAge * Ghf_ij / (rho_c_br * dzbEQ);

    // bedrock only; pure vertical conduction problem
    for (PetscInt k=1; k < k0; k++) {
      L[k] = -brR;
      D[k] = 1.0 + 2.0 * brR;
      U[k] = -brR;
      rhs[k] = Tb[k];
    }
  }

  // bottom part of ice (and top of bedrock in some cases): k=k0=Mbz-1 eqn
  if (ks == 0) { // no ice; set T[0] to surface temp if grounded
    if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
    D[k0] = 1.0;
    U[k0] = 0.0;
    // if floating and no ice then worry only about bedrock temps
    if (PismModMask(mask_ij) == MASK_FLOATING) {
      // essentially no ice but floating ... ask PISMOceanCoupler
      rhs[k0] = Tshelfbase_ij;
      // FIXME: split k0 into two grid points?
    } else { // top of bedrock sees atmosphere
      rhs[k0] = Ts_ij; 
    }
  } else { // ks > 0; there is ice
    planeStar ss;
    ierr = T3.getPlaneStarZ(i,j,0.0,&ss);
    const PetscScalar UpTu = (u[0] < 0) ? u[0] * (ss.ip1 -  ss.ij) / dx :
                                          u[0] * (ss.ij  - ss.im1) / dx;
    const PetscScalar UpTv = (v[0] < 0) ? v[0] * (ss.jp1 -  ss.ij) / dy :
                                          v[0] * (ss.ij  - ss.jm1) / dy;
    // for w, always difference *up* from base, but make it implicit
    if (PismModMask(mask_ij) == MASK_FLOATING) {
      // just apply Dirichlet condition to base of column of ice in an ice shelf
      if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
      D[k0] = 1.0;
      U[k0] = 0.0;
      rhs[k0] = Tshelfbase_ij; // set by PISMOceanCoupler
    } else { 
      // there is *grounded* ice; ice/bedrock interface; from FV across interface
      rhs[k0] = T[0] + dtTempAge * (Rb_ij / (rho_c_av * dzav));
      if (!isMarginal) {
        rhs[k0] += dtTempAge * rho_c_ratio * 0.5 * (Sigma[0] / rho_c_I);
        // WARNING: subtle consequences of finite volume argument across interface
        rhs[k0] -= dtTempAge * rho_c_ratio * (0.5 * (UpTu + UpTv));
      }
      const PetscScalar AA = dtTempAge * rho_c_ratio * w[0] / (2.0 * dzEQ);
      if (Mbz > 1) { // there is bedrock; apply upwinding if w[0]<0,
                     // otherwise ignore advection; note 
                     // jump in diffusivity coefficient
        L[k0] = - brReff;
        if (w[0] >= 0.0) {  // velocity upward
          D[k0] = 1.0 + iceReff + brReff;
          U[k0] = - iceReff;
        } else { // velocity downward
          D[k0] = 1.0 + iceReff + brReff - AA;
          U[k0] = - iceReff + AA;
        }
      } else { // no bedrock; apply geothermal flux here
        // L[k0] = 0.0;  (note this is not an allocated location!) 
        if (w[0] >= 0.0) {  // velocity upward
          D[k0] = 1.0 + 2.0 * iceR;
          U[k0] = - 2.0 * iceR;
        } else { // velocity downward
          D[k0] = 1.0 + 2.0 * iceR - AA;
          U[k0] = - 2.0 * iceR + AA;
        }
        rhs[k0] += 2.0 * dtTempAge * Ghf_ij / (rho_c_I * dzEQ);
      }
    }
  }

  // generic ice segment: build k0+1:k0+ks-1 eqns
  for (PetscInt k = 1; k < ks; k++) {
    planeStar ss;
    ierr = T3.getPlaneStarZ(i,j,k * dzEQ,&ss);
    const PetscScalar UpTu = (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx :
                                          u[k] * (ss.ij  - ss.im1) / dx;
    const PetscScalar UpTv = (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy :
                                          v[k] * (ss.ij  - ss.jm1) / dy;
    const PetscScalar AA = nuEQ * w[k];      
    if (w[k] >= 0.0) {  // velocity upward
      L[k0+k] = - iceR - AA * (1.0 - lambda/2.0);
      D[k0+k] = 1.0 + 2.0 * iceR + AA * (1.0 - lambda);
      U[k0+k] = - iceR + AA * (lambda/2.0);
    } else {  // velocity downward
      L[k0+k] = - iceR - AA * (lambda/2.0);
      D[k0+k] = 1.0 + 2.0 * iceR - AA * (1.0 - lambda);
      U[k0+k] = - iceR + AA * (1.0 - lambda/2.0);
    }
    rhs[k0+k] = T[k];
    if (!isMarginal) {
      rhs[k0+k] += dtTempAge * (Sigma[k] / rho_c_I - UpTu - UpTv);
    }
  }
      
  // surface b.c.
  if (ks>0) {
    L[k0+ks] = 0.0;
    D[k0+ks] = 1.0;
    // ignore U[k0+ks]
    rhs[k0+ks] = Ts_ij;
  }

  // solve it; note melting not addressed yet
  return solveTridiagonalSystem(k0+ks+1,x);
}


// Copyright (C) 2004-2011, 2013 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "pism_const.hh"
#include "iceModelVec.hh"
#include "tempSystem.hh"
#include "Mask.hh"

tempSystemCtx::tempSystemCtx(PetscInt my_Mz, string my_prefix)
      : columnSystemCtx(my_Mz, my_prefix) {
  Mz = my_Mz;
  // set flags to indicate nothing yet set
  initAllDone = false;
  schemeParamsValid = false;
  surfBCsValid = false;
  basalBCsValid = false;
  // set values so we can check if init was called on all
  dx = -1;
  dy = -1;
  dtTemp = -1;
  dzEQ = -1;
  ice_rho = -1;
  ice_c_p = -1;
  ice_k = -1;
  T = NULL;
  u = NULL;
  v = NULL;
  w = NULL;
  strain_heating = NULL;
  T3 = NULL;
}


PetscErrorCode tempSystemCtx::initAllColumns() {
  // check whether each parameter & pointer got set
  if (dx <= 0.0) { SETERRQ(PETSC_COMM_SELF, 2,"un-initialized dx in tempSystemCtx"); }
  if (dy <= 0.0) { SETERRQ(PETSC_COMM_SELF, 3,"un-initialized dy in tempSystemCtx"); }
  if (dtTemp <= 0.0) { SETERRQ(PETSC_COMM_SELF, 4,"un-initialized dtTemp in tempSystemCtx"); }
  if (dzEQ <= 0.0) { SETERRQ(PETSC_COMM_SELF, 5,"un-initialized dzEQ in tempSystemCtx"); }
  if (ice_rho <= 0.0) { SETERRQ(PETSC_COMM_SELF, 7,"un-initialized ice_rho in tempSystemCtx"); }
  if (ice_c_p <= 0.0) { SETERRQ(PETSC_COMM_SELF, 8,"un-initialized ice_c_p in tempSystemCtx"); }
  if (ice_k <= 0.0) { SETERRQ(PETSC_COMM_SELF, 9,"un-initialized ice_k in tempSystemCtx"); }
  if (T == NULL) { SETERRQ(PETSC_COMM_SELF, 13,"un-initialized pointer T in tempSystemCtx"); }
  if (u == NULL) { SETERRQ(PETSC_COMM_SELF, 15,"un-initialized pointer u in tempSystemCtx"); }
  if (v == NULL) { SETERRQ(PETSC_COMM_SELF, 16,"un-initialized pointer v in tempSystemCtx"); }
  if (w == NULL) { SETERRQ(PETSC_COMM_SELF, 17,"un-initialized pointer w in tempSystemCtx"); }
  if (strain_heating == NULL) { SETERRQ(PETSC_COMM_SELF, 18,"un-initialized pointer strain_heating in tempSystemCtx"); }
  if (T3 == NULL) { SETERRQ(PETSC_COMM_SELF, 19,"un-initialized pointer T3 in tempSystemCtx"); }
  // set derived constants
  nuEQ = dtTemp / dzEQ;
  rho_c_I = ice_rho * ice_c_p;
  iceK = ice_k / rho_c_I;
  iceR = iceK * dtTemp / PetscSqr(dzEQ);
  // done
  initAllDone = true;
  return 0;
}


PetscErrorCode tempSystemCtx::setSchemeParamsThisColumn(
                     PismMask my_mask, bool my_isMarginal, PetscScalar my_lambda) {
  if (!initAllDone) {  SETERRQ(PETSC_COMM_SELF, 2,
     "setSchemeParamsThisColumn() should only be called after initAllColumns() in tempSystemCtx"); }
  if (schemeParamsValid) {  SETERRQ(PETSC_COMM_SELF, 3,
     "setSchemeParamsThisColumn() called twice (?) in tempSystemCtx"); }
  mask = my_mask;
  isMarginal = my_isMarginal;
  lambda = my_lambda;
  schemeParamsValid = true;
  return 0;
}


PetscErrorCode tempSystemCtx::setSurfaceBoundaryValuesThisColumn(PetscScalar my_Ts) {
  if (!initAllDone) {  SETERRQ(PETSC_COMM_SELF, 2,
     "setSurfaceBoundaryValuesThisColumn() should only be called after initAllColumns() in tempSystemCtx"); }
  if (surfBCsValid) {  SETERRQ(PETSC_COMM_SELF, 3,
     "setSurfaceBoundaryValuesThisColumn() called twice (?) in tempSystemCtx"); }
  Ts = my_Ts;
  surfBCsValid = true;
  return 0;
}


PetscErrorCode tempSystemCtx::setBasalBoundaryValuesThisColumn(
                     PetscScalar my_G0, PetscScalar my_Tshelfbase, PetscScalar my_Rb) {
  if (!initAllDone) {  SETERRQ(PETSC_COMM_SELF, 2,
     "setIndicesThisColumn() should only be called after initAllColumns() in tempSystemCtx"); }
  if (basalBCsValid) {  SETERRQ(PETSC_COMM_SELF, 3,
     "setBasalBoundaryValuesThisColumn() called twice (?) in tempSystemCtx"); }
  G0 = my_G0;
  Tshelfbase = my_Tshelfbase;
  Rb = my_Rb;
  basalBCsValid = true;
  return 0;
}


PetscErrorCode tempSystemCtx::solveThisColumn(PetscScalar **x, PetscErrorCode &pivoterrorindex) {

  if (!initAllDone) {  SETERRQ(PETSC_COMM_SELF, 2,
     "solveThisColumn() should only be called after initAllColumns() in tempSystemCtx"); }
  if (!schemeParamsValid) {  SETERRQ(PETSC_COMM_SELF, 3,
     "solveThisColumn() should only be called after setSchemeParamsThisColumn() in tempSystemCtx"); }
  if (!surfBCsValid) {  SETERRQ(PETSC_COMM_SELF, 3,
     "solveThisColumn() should only be called after setSurfaceBoundaryValuesThisColumn() in tempSystemCtx"); }
  if (!basalBCsValid) {  SETERRQ(PETSC_COMM_SELF, 3,
     "solveThisColumn() should only be called after setBasalBoundaryValuesThisColumn() in tempSystemCtx"); }

  Mask M;

  // bottom of ice; k=0 eqn
  if (ks == 0) { // no ice; set T[0] to surface temp if grounded
    // note L[0] not allocated 
    D[0] = 1.0;
    U[0] = 0.0;
    // if floating and no ice then worry only about bedrock temps
    if (M.ocean(mask)) {
      // essentially no ice but floating ... ask PISMOceanCoupler
      rhs[0] = Tshelfbase;
    } else { // top of bedrock sees atmosphere
      rhs[0] = Ts; 
    }
  } else { // ks > 0; there is ice
    // for w, always difference *up* from base, but make it implicit
    if (M.ocean(mask)) {
      // just apply Dirichlet condition to base of column of ice in an ice shelf
      // note L[0] not allocated 
      D[0] = 1.0;
      U[0] = 0.0;
      rhs[0] = Tshelfbase; // set by PISMOceanCoupler
    } else { 
      // there is *grounded* ice; from FV across interface
      rhs[0] = T[0] + dtTemp * (Rb / (rho_c_I * dzEQ));
      if (!isMarginal) {
        rhs[0] += dtTemp * 0.5 * strain_heating[0]/ rho_c_I;
        planeStar<PetscScalar> ss;
        T3->getPlaneStar(i,j,0,&ss);
        const PetscScalar UpTu = (u[0] < 0) ? u[0] * (ss.e -  ss.ij) / dx :
                                              u[0] * (ss.ij  - ss.w) / dx;
        const PetscScalar UpTv = (v[0] < 0) ? v[0] * (ss.n -  ss.ij) / dy :
                                              v[0] * (ss.ij  - ss.s) / dy;
        rhs[0] -= dtTemp  * (0.5 * (UpTu + UpTv));
      }
      // vertical upwinding
      // L[0] = 0.0;  (is not an allocated location!) 
      D[0] = 1.0 + 2.0 * iceR;
      U[0] = - 2.0 * iceR;
      if (w[0] < 0.0) { // velocity downward: add velocity contribution
        const PetscScalar AA = dtTemp * w[0] / (2.0 * dzEQ);
        D[0] -= AA;
        U[0] += AA;
      }
      // apply geothermal flux G0 here
      rhs[0] += 2.0 * dtTemp * G0 / (rho_c_I * dzEQ);
    }
  }

  // generic ice segment; build 1:ks-1 eqns
  for (PetscInt k = 1; k < ks; k++) {
    planeStar<PetscScalar> ss;
    T3->getPlaneStar_fine(i,j,k,&ss);
    const PetscScalar UpTu = (u[k] < 0) ? u[k] * (ss.e -  ss.ij) / dx :
                                          u[k] * (ss.ij  - ss.w) / dx;
    const PetscScalar UpTv = (v[k] < 0) ? v[k] * (ss.n -  ss.ij) / dy :
                                          v[k] * (ss.ij  - ss.s) / dy;
    const PetscScalar AA = nuEQ * w[k];      
    if (w[k] >= 0.0) {  // velocity upward
      L[k] = - iceR - AA * (1.0 - lambda/2.0);
      D[k] = 1.0 + 2.0 * iceR + AA * (1.0 - lambda);
      U[k] = - iceR + AA * (lambda/2.0);
    } else {  // velocity downward
      L[k] = - iceR - AA * (lambda/2.0);
      D[k] = 1.0 + 2.0 * iceR - AA * (1.0 - lambda);
      U[k] = - iceR + AA * (1.0 - lambda/2.0);
    }
    rhs[k] = T[k];
    if (!isMarginal) {
      rhs[k] += dtTemp * (strain_heating[k] / rho_c_I - UpTu - UpTv);
    }
  }
      
  // surface b.c.
  if (ks>0) {
    L[ks] = 0.0;
    D[ks] = 1.0;
    // ignore U[ks]
    rhs[ks] = Ts;
  }

  // mark column as done
  schemeParamsValid = false;
  surfBCsValid = false;
  basalBCsValid = false;

  // solve it; note melting not addressed yet
  pivoterrorindex = solveTridiagonalSystem(ks+1,x);
  return 0;
}


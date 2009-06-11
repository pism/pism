// Copyright (C) 2009 Andreas Aschwanden and Ed Bueler
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

#include "iceEnthalpyModel.hh"
#include "columnSystem.hh"

#define DEBUGVERB 2

//! Tridiagonal linear system for vertical column of enthalpy-based conservation of energy.
class enthSystemCtx : public columnSystemCtx {

public:
  enthSystemCtx(int my_Mz, int my_Mbz);
  PetscErrorCode initAllColumns();
  PetscErrorCode setIndicesThisColumn(PetscInt i, PetscInt j, PetscInt ks);  
  PetscErrorCode setSchemeParamsThisColumn(
                     PetscScalar my_mask, bool my_isMarginal, PetscScalar my_lambda);  
  PetscErrorCode setSurfaceBoundaryValuesThisColumn(PetscScalar my_Ts);
  PetscErrorCode setBasalBoundaryValuesThisColumn(
                     PetscScalar my_Ghf, PetscScalar my_Tshelfbase, PetscScalar my_Rb);
  PetscErrorCode solveThisColumn(PetscScalar **x);  

public:
  // constants which should be set before calling initForAllColumns()
  PetscScalar  dx,
               dy,
               dtTemp,
               dzEQ,
               dzbEQ,
               ice_rho,
               ice_c,
               ice_k,
               bed_thermal_rho,
               bed_thermal_c,
               bed_thermal_k;
  // pointers which should be set before calling initForAllColumns()
  PetscScalar  *Enth,
               *Tb,
               *u,
               *v,
               *w,
               *Sigma;
  IceModelVec3 *Enth3;

protected: // used internally
  PetscInt    Mz, Mbz, k0;
  PetscInt    i, j, ks;
  PetscScalar mask, lambda, Ts, Ghf, Tshelfbase, Rb;
  bool        isMarginal;
  PetscScalar nuEQ,
              rho_c_I,
              rho_c_br,
              rho_c_av,
              iceK,
              iceR,
              brK,
              brR,
              rho_c_ratio,
              dzav,
              iceReff,
              brReff;
  bool        initAllDone,
              indicesValid,
              schemeParamsValid,
              surfBCsValid,
              basalBCsValid;
};


enthSystemCtx::enthSystemCtx(int my_Mz, int my_Mbz)
      : columnSystemCtx(my_Mz + my_Mbz - 1) {
  Mz = my_Mz;
  Mbz = my_Mbz;
  k0 = Mbz - 1; // max size nmax of system is Mz + k0 = Mz + Mbz - 1
  // set flags to indicate nothing yet set
  initAllDone = false;
  indicesValid = false;
  schemeParamsValid = false;
  surfBCsValid = false;
  basalBCsValid = false;
  // set values so we can check if init was called on all
  dx = -1;
  dy = -1;
  dtTemp = -1;
  dzEQ = -1;
  dzbEQ = -1;
  ice_rho = -1;
  ice_c   = -1;
  ice_k   = -1;
  bed_thermal_rho = -1;
  bed_thermal_c   = -1;
  bed_thermal_k   = -1;
  Enth = NULL;
  Tb = NULL;
  u = NULL;
  v = NULL;
  w = NULL;
  Sigma = NULL;
  Enth3 = NULL;
}


PetscErrorCode enthSystemCtx::initAllColumns() {
  // check whether each parameter & pointer got set
  if (dx <= 0.0) { SETERRQ(2,"un-initialized dx in enthSystemCtx"); }
  if (dy <= 0.0) { SETERRQ(3,"un-initialized dy in enthSystemCtx"); }
  if (dtTemp <= 0.0) { SETERRQ(4,"un-initialized dtTemp in enthSystemCtx"); }
  if (dzEQ <= 0.0) { SETERRQ(5,"un-initialized dzEQ in enthSystemCtx"); }
  if (dzbEQ <= 0.0) { SETERRQ(6,"un-initialized dzbEQ in enthSystemCtx"); }
  if (ice_rho <= 0.0) { SETERRQ(7,"un-initialized ice_rho in enthSystemCtx"); }
  if (ice_c <= 0.0) { SETERRQ(8,"un-initialized ice_c_p in enthSystemCtx"); }
  if (ice_k <= 0.0) { SETERRQ(9,"un-initialized ice_k in enthSystemCtx"); }
  if (bed_thermal_rho <= 0.0) { SETERRQ(10,"un-initialized bed_thermal_rho in enthSystemCtx"); }
  if (bed_thermal_c <= 0.0) { SETERRQ(11,"un-initialized bed_thermal_c_p in enthSystemCtx"); }
  if (bed_thermal_k <= 0.0) { SETERRQ(12,"un-initialized bed_thermal_k in enthSystemCtx"); }
  if (Enth == NULL) { SETERRQ(13,"un-initialized pointer T in enthSystemCtx"); }
  if (Tb == NULL) { SETERRQ(14,"un-initialized pointer Tb in enthSystemCtx"); }
  if (u == NULL) { SETERRQ(15,"un-initialized pointer u in enthSystemCtx"); }
  if (v == NULL) { SETERRQ(16,"un-initialized pointer v in enthSystemCtx"); }
  if (w == NULL) { SETERRQ(17,"un-initialized pointer w in enthSystemCtx"); }
  if (Sigma == NULL) { SETERRQ(18,"un-initialized pointer Sigma in enthSystemCtx"); }
  if (Enth3 == NULL) { SETERRQ(19,"un-initialized pointer T3 in enthSystemCtx"); }
  // set derived constants
  nuEQ = dtTemp / dzEQ;
  rho_c_I = ice_rho * ice_c;
  rho_c_br = bed_thermal_rho * bed_thermal_c;
  rho_c_av = (dzEQ * rho_c_I + dzbEQ * rho_c_br) / (dzEQ + dzbEQ);
  iceK = ice_k / rho_c_I;
  iceR = iceK * dtTemp / PetscSqr(dzEQ);
  brK = bed_thermal_k / rho_c_br;
  brR = brK * dtTemp / PetscSqr(dzbEQ);
  rho_c_ratio = rho_c_I / rho_c_av;
  dzav = 0.5 * (dzEQ + dzbEQ);
  iceReff = ice_k * dtTemp / (rho_c_av * dzEQ * dzEQ);
  brReff = bed_thermal_k * dtTemp / (rho_c_av * dzbEQ * dzbEQ);
  // done
  initAllDone = true;
  return 0;
}


PetscErrorCode enthSystemCtx::setIndicesThisColumn(
                  PetscInt my_i, PetscInt my_j, PetscInt my_ks) {
  if (!initAllDone) {  SETERRQ(2,
     "setIndicesThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (indicesValid) {  SETERRQ(3,
     "setIndicesThisColumn() called twice in same column (?) in enthSystemCtx"); }
  i = my_i;
  j = my_j;
  ks = my_ks;
  indicesValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::setSchemeParamsThisColumn(
                     PetscScalar my_mask, bool my_isMarginal, PetscScalar my_lambda) {
  if (!initAllDone) {  SETERRQ(2,
     "setSchemeParamsThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (schemeParamsValid) {  SETERRQ(3,
     "setSchemeParamsThisColumn() called twice (?) in enthSystemCtx"); }
  mask = my_mask;
  isMarginal = my_isMarginal;
  lambda = my_lambda;
  schemeParamsValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::setSurfaceBoundaryValuesThisColumn(PetscScalar my_Ts) {
  if (!initAllDone) {  SETERRQ(2,
     "setSurfaceBoundaryValuesThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (surfBCsValid) {  SETERRQ(3,
     "setSurfaceBoundaryValuesThisColumn() called twice (?) in enthSystemCtx"); }
  Ts = my_Ts;
  surfBCsValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::setBasalBoundaryValuesThisColumn(
                     PetscScalar my_Ghf, PetscScalar my_Tshelfbase, PetscScalar my_Rb) {
  if (!initAllDone) {  SETERRQ(2,
     "setIndicesThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (basalBCsValid) {  SETERRQ(3,
     "setBasalBoundaryValuesThisColumn() called twice (?) in enthSystemCtx"); }
  Ghf = my_Ghf;
  Tshelfbase = my_Tshelfbase;
  Rb = my_Rb;
  basalBCsValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::solveThisColumn(PetscScalar **x) {
  SETERRQ(1,"enthSystemCtx::solveThisColumn()   NOT IMPLEMENTED");
#if 0
  PetscErrorCode ierr;
  if (!initAllDone) {  SETERRQ(2,
     "solveThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (!indicesValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setIndicesThisColumn() in enthSystemCtx"); }
  if (!schemeParamsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setSchemeParamsThisColumn() in enthSystemCtx"); }
  if (!surfBCsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setSurfaceBoundaryValuesThisColumn() in enthSystemCtx"); }
  if (!basalBCsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setBasalBoundaryValuesThisColumn() in enthSystemCtx"); }

  if (Mbz > 1) { // bedrock present: build k=0:Mbz-2 eqns
    // gives O(\Delta t,\Delta z^2) convergence in Test K for equal spaced grid;
    // note L[0] not an allocated location:
    D[0] = (1.0 + 2.0 * brR);
    U[0] = - 2.0 * brR;  
    rhs[0] = Tb[0] + 2.0 * dtTemp * Ghf / (rho_c_br * dzbEQ);

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
    if (PismModMask(mask) == MASK_FLOATING) {
      // essentially no ice but floating ... ask PISMOceanCoupler
      rhs[k0] = Tshelfbase;
      // FIXME: split k0 into two grid points?
    } else { // top of bedrock sees atmosphere
      rhs[k0] = Ts; 
    }
  } else { // ks > 0; there is ice
    planeStar ss;
    ierr = T3->getPlaneStarZ(i,j,0.0,&ss);
    const PetscScalar UpTu = (u[0] < 0) ? u[0] * (ss.ip1 -  ss.ij) / dx :
                                          u[0] * (ss.ij  - ss.im1) / dx;
    const PetscScalar UpTv = (v[0] < 0) ? v[0] * (ss.jp1 -  ss.ij) / dy :
                                          v[0] * (ss.ij  - ss.jm1) / dy;
    // for w, always difference *up* from base, but make it implicit
    if (PismModMask(mask) == MASK_FLOATING) {
      // just apply Dirichlet condition to base of column of ice in an ice shelf
      if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
      D[k0] = 1.0;
      U[k0] = 0.0;
      rhs[k0] = Tshelfbase; // set by PISMOceanCoupler
    } else { 
      // there is *grounded* ice; ice/bedrock interface; from FV across interface
      rhs[k0] = T[0] + dtTemp * (Rb / (rho_c_av * dzav));
      if (!isMarginal) {
        rhs[k0] += dtTemp * rho_c_ratio * 0.5 * (Sigma[0] / rho_c_I);
        // WARNING: subtle consequences of finite volume argument across interface
        rhs[k0] -= dtTemp * rho_c_ratio * (0.5 * (UpTu + UpTv));
      }
      const PetscScalar AA = dtTemp * rho_c_ratio * w[0] / (2.0 * dzEQ);
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
        rhs[k0] += 2.0 * dtTemp * Ghf / (rho_c_I * dzEQ);
      }
    }
  }

  // generic ice segment: build k0+1:k0+ks-1 eqns
  for (PetscInt k = 1; k < ks; k++) {
    planeStar ss;
    ierr = T3->getPlaneStarZ(i,j,k * dzEQ,&ss);
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
      rhs[k0+k] += dtTemp * (Sigma[k] / rho_c_I - UpTu - UpTv);
    }
  }
      
  // surface b.c.
  if (ks>0) {
    L[k0+ks] = 0.0;
    D[k0+ks] = 1.0;
    // ignore U[k0+ks]
    rhs[k0+ks] = Ts;
  }

  // mark column as done
  indicesValid = false;
  schemeParamsValid = false;
  surfBCsValid = false;
  basalBCsValid = false;

  // solve it; note melting not addressed yet
  return solveTridiagonalSystem(k0+ks+1,x);
#endif
  return 999;
}


PetscErrorCode IceEnthalpyModel::enthalpyStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount) {

  SETERRQ(999,"IceEnthalpyModel::enthalpyStep() NOT IMPLEMENTED");

#if 0
  if (doColdIceMethods) {
    PetscPrintf(grid.com,
      "\n\n    IceEnthalpyModel::enthalpyStep() called but doColdIceMethods==true ... ending\n");
    PetscEnd();
  }
  
  PetscErrorCode  ierr;

  // set up fine grid in ice and bedrock
  PetscInt    fMz, fMbz;
  PetscScalar fdz, *fzlev, fdzb, *fzblev;
  ierr = grid.getFineEqualVertCounts(fMz,fMbz); CHKERRQ(ierr);
  fzlev = new PetscScalar[fMz];
  fzblev = new PetscScalar[fMbz];
  ierr = grid.getFineEqualVertLevs(fMz,fMbz,fdz,fdzb,fzlev,fzblev); CHKERRQ(ierr);

  ierr = verbPrintf(DEBUGVERB,grid.com,
    "\n  [entering enthalpyStep(); fMz = %d, fdz = %5.3f, fMbz = %d, fdzb = %5.3f]",
    fMz, fdz, fMbz, fdzb); CHKERRQ(ierr);

  enthSystemCtx system(fMz,fMbz);
  system.dx              = grid.dx;
  system.dy              = grid.dy;
  system.dtTemp          = dtTempAge; // same time step for temp and age, currently
  system.dzEQ            = fdz;
  system.dzbEQ           = fdzb;
  system.ice_rho         = config.get("ice_density"); // ice->rho;
  system.ice_c           = config.get("ice_specific_heat_capacity"); // ice->c_p;
  system.ice_k           = config.get("ice_thermal_conductivity"); // ice->k;
  system.bed_thermal_rho = config.get("bedrock_thermal_density"); // bed_thermal.rho;
  system.bed_thermal_c   = config.get("bedrock_thermal_specific_heat_capacity"); // bed_thermal.c_p;
  system.bed_thermal_k   = config.get("bedrock_thermal_conductivity"); // bed_thermal.k;

  const PetscInt k0 = fMbz - 1;
  PetscScalar *x;  
  x = new PetscScalar[fMz + k0]; // space for solution of system; length = fMz + fMbz - 1 

  // constants needed after solution of system, in insertion phase
  const PetscScalar rho_c_I = system.ice_rho * system.ice_c,
                    rho_c_br = system.bed_thermal_rho * system.bed_thermal_c,
                    rho_c_av = (fdz * rho_c_I + fdzb * rho_c_br) / (fdz + fdzb);
  // this is bulge limit constant in J kg-1; is max amount by which ice
  //   or bedrock can be lower than surface temperature
  const PetscScalar bulgeMax   = system.ice_c * 15.0;  // enthalpy change equivalent to change in cold ice temp by 15 K

  PetscScalar *Enthnew, *Tbnew;
  // pointers to values in current column
  system.u     = new PetscScalar[fMz];
  system.v     = new PetscScalar[fMz];
  system.w     = new PetscScalar[fMz];
  system.Sigma = new PetscScalar[fMz];
  system.Enth  = new PetscScalar[fMz];
  Enthnew      = new PetscScalar[fMz];

  system.Tb    = new PetscScalar[fMbz];
  Tbnew        = new PetscScalar[fMbz];
  
  // system needs access to Enth3 for planeStar()
  system.Enth3 = &Enth3;

  // checks that all needed constants and pointers got set:
  ierr = system.initAllColumns(); CHKERRQ(ierr);

  // now get map-plane fields, starting with coupler fields
  PetscScalar  **Ts, **Tshelfbase, **H, **Ghf, **mask, **Hmelt, **Rb,
               **basalMeltRate, **bmr_float;

  IceModelVec2 *pccTs, *pccsbt, *pccsbmf;
  if (atmosPCC != PETSC_NULL) {
    // call sets pccTs to point to IceModelVec2 with current surface temps
    ierr = atmosPCC->updateSurfTempAndProvide(
              grid.year, dtTempAge / secpera, &info_coupler, pccTs);
              CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");
  }
  if (oceanPCC != PETSC_NULL) {
    ierr = oceanPCC->updateShelfBaseTempAndProvide(
              grid.year, dt / secpera, &info_coupler, pccsbt);
              CHKERRQ(ierr);
    ierr = oceanPCC->updateShelfBaseMassFluxAndProvide(
              grid.year, dt / secpera, &info_coupler, pccsbmf);
              CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: oceanPCC == PETSC_NULL");
  }
  ierr = pccTs->get_array(Ts);  CHKERRQ(ierr);
  ierr = pccsbt->get_array(Tshelfbase);  CHKERRQ(ierr);
  ierr = pccsbmf->get_array(bmr_float);  CHKERRQ(ierr);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vHmelt.get_array(Hmelt); CHKERRQ(ierr);
  ierr = vbasalMeltRate.get_array(basalMeltRate); CHKERRQ(ierr);
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vRb.get_array(Rb); CHKERRQ(ierr);
  ierr = vGhf.get_array(Ghf); CHKERRQ(ierr);

  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = Sigma3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = EnthNew3.begin_access(); CHKERRQ(ierr);

  ierr = Tb3.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // this should *not* be replaced by call to grid.kBelowHeight():
      const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/fdz));
      
      if (k0+ks>0) { // if there are enough points in bedrock&ice to bother ...
        ierr = system.setIndicesThisColumn(i,j,ks); CHKERRQ(ierr);
        ierr = Tb3.getValColumn(i,j,fMbz,fzblev,system.Tb); CHKERRQ(ierr);

        if (grid.vertical_spacing == EQUAL) {
          ierr = u3.getValColumnPL(i,j,fMz,fzlev,system.u); CHKERRQ(ierr);
          ierr = v3.getValColumnPL(i,j,fMz,fzlev,system.v); CHKERRQ(ierr);
          ierr = w3.getValColumnPL(i,j,fMz,fzlev,system.w); CHKERRQ(ierr);
          ierr = Sigma3.getValColumnPL(i,j,fMz,fzlev,system.Sigma); CHKERRQ(ierr);
          ierr = Enth3.getValColumnPL(i,j,fMz,fzlev,system.Enth); CHKERRQ(ierr);
        } else {
          // slower, but right for not-equal spaced
          ierr = u3.getValColumnQUAD(i,j,fMz,fzlev,system.u); CHKERRQ(ierr);
          ierr = v3.getValColumnQUAD(i,j,fMz,fzlev,system.v); CHKERRQ(ierr);
          ierr = w3.getValColumnQUAD(i,j,fMz,fzlev,system.w); CHKERRQ(ierr);
          ierr = Sigma3.getValColumnQUAD(i,j,fMz,fzlev,system.Sigma); CHKERRQ(ierr);
          ierr = Enth3.getValColumnQUAD(i,j,fMz,fzlev,system.Enth); CHKERRQ(ierr);
        }

        // FIXME: following mechanism only make sense for comparing cold ice diffusion
        //   to advection; temperate ice diffusion using moisture transport diffusion
        //   coefficient might give different result
        // go through column and find appropriate lambda for BOMBPROOF
        PetscScalar lambda = 1.0;  // start with centered implicit for more accuracy
        for (PetscInt k = 1; k < ks; k++) {   
          const PetscScalar denom = (PetscAbs(system.w[k]) + 0.000001/secpera)
                                      * rho_c_I * fdz;  
          lambda = PetscMin(lambda, 2.0 * system.ice_k / denom);
        }
        if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1
        // if isMarginal then only do vertical conduction for ice; ignore advection
        //   and strain heating if isMarginal
        const bool isMarginal = checkThinNeigh(H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                               H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);
        ierr = system.setSchemeParamsThisColumn(mask[i][j], isMarginal, lambda);
                 CHKERRQ(ierr);  

        // set boundary values for tridiagonal system
        ierr = system.setSurfaceBoundaryValuesThisColumn(Ts[i][j]); CHKERRQ(ierr);
        ierr = system.setBasalBoundaryValuesThisColumn(
                 Ghf[i][j],Tshelfbase[i][j],Rb[i][j]); CHKERRQ(ierr);

        // solve the system for this column; melting not addressed yet
        ierr = system.solveThisColumn(&x); // no CHKERRQ(ierr) immediately because:
        if (ierr > 0) {
          SETERRQ3(2,
            "Tridiagonal solve failed at (%d,%d) with zero pivot position %d.\n",
            i, j, ierr);
        } else { CHKERRQ(ierr); }
      }

      // insert bedrock solution; check for too low below
      for (PetscInt k=0; k < k0; k++) {
        Tbnew[k] = x[k];  // FIXME: CAREFUL HERE: is system getting temperature in bedrock into x?
      }

      // prepare for melting/refreezing
      PetscScalar Hmeltnew = Hmelt[i][j];
      
      // insert solution for generic ice segments
      for (PetscInt k=1; k <= ks; k++) {
        if (allowAboveMelting == PETSC_TRUE) {
          Enthnew[k] = x[k0 + k];  // FIXME: CAREFUL HERE: is system getting enthalpy in ice into x?
        } else {
          PetscScalar T_m, H_l, H_s;
          // FIXME: next line only uses H_s, not T_m, H_l; no check for liquid water
          getEnthalpyInterval(config, getPressureFromDepth(config, H[i][j] - fzlev[k]), T_m, H_l, H_s); 
          // OLD:  if (x[k0 + k] > Tpmp) {
          if (x[k0 + k] > H_s) {
FIXME FROM HERE: use drainageToBaseModelEnth(H[i][j],fzlev[k],fdz,L_latent,enthalpyTOUPDATE,HmeltTOUPDATE)
            Tnew[k] = Tpmp;
            PetscScalar Texcess = x[k0 + k] - Tpmp; // always positive
            excessToFromBasalMeltLayer(rho_c_I, fzlev[k], fdz, &Texcess, &Hmeltnew);
            // Texcess  will always come back zero here; ignore it
          } else {
            Tnew[k] = x[k0 + k];
          }
        }
        if (Tnew[k] < Ts[i][j] - bulgeMax) {
          Tnew[k] = Ts[i][j] - bulgeMax;  bulgeCount++;   }
      }
      
      // insert solution for ice/rock interface (or base of ice shelf) segment
      if (ks > 0) {
        if (allowAboveMelting == PETSC_TRUE) {
          Tnew[0] = x[k0];
        } else {  // compute diff between x[k0] and Tpmp; melt or refreeze as appropriate
          const PetscScalar Tpmp = ice->meltingTemp - ice->beta_CC_grad * H[i][j];
          PetscScalar Texcess = x[k0] - Tpmp; // positive or negative
          if (PismModMask(mask[i][j]) == MASK_FLOATING) {
             // when floating, only half a segment has had its temperature raised
             // above Tpmp
             excessToFromBasalMeltLayer(rho_c_I/2, 0.0, fdz, &Texcess, &Hmeltnew);
          } else {
             excessToFromBasalMeltLayer(rho_c_av, 0.0, fdz, &Texcess, &Hmeltnew);
          }
          Tnew[0] = Tpmp + Texcess;
          if (Tnew[0] > (Tpmp + 0.00001)) {
            SETERRQ(1,"updated temperature came out above Tpmp");
          }
        }
        if (Tnew[0] < Ts[i][j] - bulgeMax) {
          Tnew[0] = Ts[i][j] - bulgeMax;   bulgeCount++;   }
      } else {
        Hmeltnew = 0.0;
      }
      
      // we must agree on redundant values T(z=0) at top of bedrock and at bottom of ice
      if (ks > 0) {
        Tbnew[k0] = Tnew[0];
      } else {
        if (PismModMask(mask[i][j]) == MASK_FLOATING) { // top of bedrock sees ocean
          Tbnew[k0] = Tshelfbase[i][j]; // set by PISMOceanCoupler
        } else { // top of bedrock sees atmosphere
          Tbnew[k0] = Ts[i][j];
        }
      }
      // check bedrock solution        
      for (PetscInt k=0; k <= k0; k++) {
        if (Tbnew[k] < Ts[i][j] - bulgeMax) {
          Tbnew[k] = Ts[i][j] - bulgeMax;   bulgeCount++;   }
      }

      // transfer column into Tb3; neighboring columns will not reference!
      ierr = Tb3.setValColumn(i,j,fMbz,fzblev,Tbnew); CHKERRQ(ierr);

      // set to air temp above ice
      for (PetscInt k=ks; k<fMz; k++) {
        Tnew[k] = Ts[i][j];
      }

      // transfer column into Tnew3; communication later
      ierr = Tnew3.setValColumnPL(i,j,fMz,fzlev,Tnew); CHKERRQ(ierr);

      // basalMeltRate[][] is rate of mass loss at bottom of ice everywhere;
      //   note massContExplicitStep() calls PISMOceanCoupler separately
      if (PismModMask(mask[i][j]) == MASK_FLOATING) {
        // rate of mass loss at bottom of ice shelf;  can be negative (marine freeze-on)
        basalMeltRate[i][j] = bmr_float[i][j]; // set by PISMOceanCoupler
      } else {
        // rate of change of Hmelt[][];  can be negative (till water freeze-on)
        basalMeltRate[i][j] = (Hmeltnew - Hmelt[i][j]) / dtTempAge;
      }

      if (PismModMask(mask[i][j]) == MASK_FLOATING) {
        // eliminate basal lubrication water if floating; 
        Hmelt[i][j] = 0.0;
      } else {
        // limit Hmelt by default max and store
        Hmelt[i][j] = PetscMin(Hmelt_max, Hmeltnew);
      }

    } 
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vHmelt.end_access(); CHKERRQ(ierr);
  ierr = vRb.end_access(); CHKERRQ(ierr);
  ierr = vGhf.end_access(); CHKERRQ(ierr);
  ierr = vbasalMeltRate.end_access(); CHKERRQ(ierr);

  ierr = pccTs->end_access(); CHKERRQ(ierr);
  ierr = pccsbt->end_access();  CHKERRQ(ierr);
  ierr = pccsbmf->end_access();  CHKERRQ(ierr);

  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = Sigma3.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = EnthNew3.end_access(); CHKERRQ(ierr);
  
  delete [] x;
  delete [] system.Enth;  delete [] system.Tb;  
  delete [] system.u;  delete [] system.v;  delete [] system.w;
  delete [] system.Sigma;
  
  delete [] Tbnew;  delete [] Enthnew;

  delete [] fzlev;   delete [] fzblev;

#endif
  return 0;
}



//! Move some of the liquid water fraction to the base according to heuristics.
PetscErrorCode IceEnthalpyModel::drainageToBaseModelEnth(
                const PetscScalar thickness, const PetscScalar z, const PetscScalar dz, const PetscScalar L,
                PetscScalar &enthalpy, PetscScalar &Hmelt) {

  SETERRQ(99,"drainageToBaseModelEnth() NOT IMPLEMENTED");

#if 0
  if (allowAboveMelting == PETSC_TRUE) {
    SETERRQ(1,"drainageToBaseModelEnth() called but allowAboveMelting==TRUE");
  }

  // meltlayerequiv is ice-equivalent thickness of water in current layer that is
  //   present in ice matrix in [z,z+dz] part of column
  const PetscScalar p              = getPressureFromDepth(config, thickness - z),
                    meltlayerequiv = getWaterFraction(config, enthalpy, p) * dz;

  if (meltlayerequiv >= 0.0) {
    if (updateHmelt == PETSC_TRUE) {
      // if low enough in the ice column, a fraction of liquid content at depth in column
      //   is transported to base
      const PetscScalar FRACTION_TO_BASE = (z < 100.0) ? 0.2 * (100.0 - z) / 100.0 : 0.0;
      // note: ice-equiv thickness:
      melttobase = FRACTION_TO_BASE * meltlayerequiv;
      Hmelt    += melttobase;  
      enthalpy -= melttobase * L;
    }
  } else if (updateHmelt == PETSC_TRUE) {  // no liquid fraction, but we could freeze some of the
                                           // basal layer
    // Texcess negative; only refreeze (i.e. reduce Hmelt) if at base and Hmelt > 0.0
    // note massmelted is NEGATIVE!
    if (z <= 0.00001) { // only freeze on if at base
      FIXME FROM HERE: in case of cold ice base we freeze on; use \Delta enthalpy = c \Delta T;
      if (Hmelt > 0.0) {
        const PetscScalar thicknessToFreezeOn = - massmelted / (ice->rho * darea);
        if (thicknessToFreezeOn <= *Hmelt) { // the water *is* available to freeze on
          *Hmelt -= thicknessToFreezeOn;
          *Texcess = 0.0;
        } else { // only refreeze Hmelt thickness of water; update Texcess
          *Hmelt = 0.0;
          const PetscScalar dTemp = ice->latentHeat * ice->rho * (*Hmelt) / (rho_c * dz);
          *Texcess += dTemp;
        }
      } 
    }
    // note: if *Hmelt == 0 and Texcess < 0.0 then Texcess unmolested; temp will go down
  }
#endif

  return 0;
}




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

#include "pism_const.hh"        // e.g. for MASK_FLOATING and PismModMask()
#include "enthalpyConverter.hh"
#include "enthColumnSystem.hh"


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
  ice_nu  = -1;
  bed_thermal_rho = -1;
  bed_thermal_c   = -1;
  bed_thermal_k   = -1;
  Enth = NULL;
  Enth_s = NULL;
  Enth_b = NULL;
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
  if (ice_nu < 0.0) { SETERRQ(21,"un-initialized ice_nu in enthSystemCtx"); } // 0.0 is allowed
  if (bed_thermal_rho <= 0.0) { SETERRQ(10,"un-initialized bed_thermal_rho in enthSystemCtx"); }
  if (bed_thermal_c <= 0.0) { SETERRQ(11,"un-initialized bed_thermal_c_p in enthSystemCtx"); }
  if (bed_thermal_k <= 0.0) { SETERRQ(12,"un-initialized bed_thermal_k in enthSystemCtx"); }
  if (Enth == NULL) { SETERRQ(13,"un-initialized pointer Enth in enthSystemCtx"); }
  if (Enth_s == NULL) { SETERRQ(20,"un-initialized pointer Enth_s in enthSystemCtx"); }
  if (Enth_b == NULL) { SETERRQ(14,"un-initialized pointer Enth_b in enthSystemCtx"); }
  if (u == NULL) { SETERRQ(15,"un-initialized pointer u in enthSystemCtx"); }
  if (v == NULL) { SETERRQ(16,"un-initialized pointer v in enthSystemCtx"); }
  if (w == NULL) { SETERRQ(17,"un-initialized pointer w in enthSystemCtx"); }
  if (Sigma == NULL) { SETERRQ(18,"un-initialized pointer Sigma in enthSystemCtx"); }
  if (Enth3 == NULL) { SETERRQ(19,"un-initialized pointer Enth3 in enthSystemCtx"); }
  // set derived constants
  nuEQ = dtTemp / dzEQ;
  rho_c_I = ice_rho * ice_c;
  rho_c_br = bed_thermal_rho * bed_thermal_c;
  iceK = ice_k / rho_c_I;
  iceRcold = iceK * dtTemp / PetscSqr(dzEQ);
  iceRtemp = ice_nu * dtTemp / PetscSqr(dzEQ);
  brK = bed_thermal_k / rho_c_br;
  brR = brK * dtTemp / PetscSqr(dzbEQ);
  dzav = 0.5 * (dzEQ + dzbEQ);
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


PetscErrorCode enthSystemCtx::setSurfaceBoundaryValuesThisColumn(PetscScalar my_Enth_surface) {
  if (!initAllDone) {  SETERRQ(2,
     "setSurfaceBoundaryValuesThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (surfBCsValid) {  SETERRQ(3,
     "setSurfaceBoundaryValuesThisColumn() called twice (?) in enthSystemCtx"); }
  Enth_ks = my_Enth_surface;
  surfBCsValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::setBasalBoundaryValuesThisColumn(
                     PetscScalar my_Ghf, PetscScalar my_Enth_shelfbase, PetscScalar my_Rb) {
  if (!initAllDone) {  SETERRQ(2,
     "setIndicesThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (basalBCsValid) {  SETERRQ(3,
     "setBasalBoundaryValuesThisColumn() called twice (?) in enthSystemCtx"); }
  Ghf = my_Ghf;
  Enth_shelfbase = my_Enth_shelfbase;
  Rb = my_Rb;
  basalBCsValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::viewConstants(PetscViewer viewer) {
  PetscErrorCode ierr;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }

  PetscTruth iascii;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for enthSystemCtx::viewConstants()\n"); }
  
  ierr = PetscViewerASCIIPrintf(viewer,"\n\n<<VIEWING enthSystemCtx:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  initAllDone = %d\n",initAllDone); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  dx,dy,dtTemp,dzEQ,dzbEQ = %8.2f,%8.2f,%10.3e,%8.2f,%8.2f\n",
                     dx,dy,dtTemp,dzEQ,dzbEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  ice_rho,ice_c,ice_k = %10.3e,%10.3e,%10.3e\n",
                     ice_rho,ice_c,ice_k); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  bed_thermal_rho,bed_thermal_c,bed_thermal_k = %10.3e,%10.3e,%10.3e\n",
                     bed_thermal_rho,bed_thermal_c,bed_thermal_k); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  nuEQ,dzav,rho_c_I,rho_c_br, = %10.3e,%10.3e,%10.3e,%10.3e\n",
                     nuEQ,dzav,rho_c_I,rho_c_br); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  iceK,iceRcold,iceRtemp,brK,brR = %10.3e,%10.3e,%10.3e,%10.3e,%10.3e\n",
		     iceK,iceRcold,iceRtemp,brK,brR); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  i,j,ks = %d,%d,%d\n",
                     i,j,ks); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  Enth_ks,Ghf,Enth_shelfbase,Rb = %10.3e,%10.3e,%10.3e,%10.3e\n",
                     Enth_ks,Ghf,Enth_shelfbase,Rb); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,">>\n\n"); CHKERRQ(ierr);
  return 0;
}


//! Solve the tridiagonal system which determines the new values of enthalpy in the column of ice and bedrock.
/*!
Highly nontrivial choices are made here, especially at the ice/bedrock
transition.

See the page documenting BOMBPROOF.  We implement equations (\ref bombtwo), (\ref bedrockeqn),
(\ref geothermalbedeqn), (\ref icebedfinalcold), (\ref icebedfinaltemperate),
(\ref neartopofbedrock), and (\ref icebasenobedrock).
 */
PetscErrorCode enthSystemCtx::solveThisColumn(PetscScalar **x) {
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

/* PRINCIPLES ABOUT THESE MODIFICATIONS: 
1)  coefficients in system are unitless and therefore most D,L,U expressions are not altered
2)  old temperature equation had units of Kelvin on each side; new equation has units
    of enthalpy, namely J kg-1, on each side
3)  item 2) means all rhs[] expressions must be modified, and expressions not proportional
    to a temperature or enthalpy generally are multiplied by c (which has units J kg-1 K-1)
*/

  if (Mbz > 1) { // bedrock present: build k=0:Mbz-2 eqns

    // case where Mbz==2 is different; it is stopped by IceEnthalpy model

    // should gives O(\Delta t,\Delta z^2) convergence
    // note L[0] not an allocated location:
    D[0] = (1.0 + 2.0 * brR);
    U[0] = - 2.0 * brR;  
    rhs[0] = Enth_b[0] + 2.0 * dtTemp * Ghf / (bed_thermal_rho * dzbEQ);

    // bedrock only; pure vertical conduction problem
    // because Mbz >=3, k0 = Mbz-1 >= 2, and the loop goes at least once
    for (PetscInt k = 1; k < k0; k++) {
      L[k] = -brR;
      D[k] = 1.0 + 2.0 * brR;
      if (k == k0-1) {
        U[k] = -brR * (bed_thermal_c / ice_c);
      } else {
        U[k] = -brR;
      }
      rhs[k] = Enth_b[k];
    }
  }

  // bottom part of ice (and top of bedrock in some cases): k=k0=Mbz-1 eqn
  if (ks == 0) { // no ice; set Enth[0] to surface enthalpy if grounded
    if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
    D[k0] = 1.0;
    U[k0] = 0.0;
    // if floating and no ice then worry only about bedrock temps
    if (PismModMask(mask) == MASK_FLOATING) {
      // essentially no ice but floating ... ask PISMOceanCoupler
      rhs[k0] = Enth_shelfbase;
    } else { // top of bedrock sees atmosphere
      rhs[k0] = Enth_ks; 
    }
  } else { // ks > 0; there is ice
    if (PismModMask(mask) == MASK_FLOATING) {
      // just apply Dirichlet condition to base of column of ice in an ice shelf
      if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
      D[k0] = 1.0;
      U[k0] = 0.0;
      rhs[k0] = Enth_shelfbase;
    } else if (Mbz == 1) {
      // grounded but no bedrock layer; apply geothermal flux here
      // WARNING: subtle consequences of finite volume argument for basal segment;
      // see BOMBPROOF docs
      const PetscScalar R = (Enth[k0] > Enth_s[k0]) ? iceRtemp : iceRcold;
      // L[k0] = 0.0;  (note this is not an allocated location!) 
      D[k0] = 1.0 + 2.0 * R;
      U[k0] = - 2.0 * R;
      if (w[0] < 0.0) { // velocity downward: upwind vertical
        const PetscScalar AA = nuEQ * w[0];
        D[k0] -= AA;
        U[k0] += AA;
      }
      rhs[k0] = Enth[k0];
      rhs[k0] += (2.0 * nuEQ / ice_rho) * (Ghf + 0.5 * Rb); // geothermal and half of frictional heat
      if (!isMarginal) {
        planeStar ss;
        ierr = Enth3->getPlaneStarZ(i,j,0.0,&ss);
        const PetscScalar UpEnthu = (u[0] < 0) ? u[0] * (ss.ip1 -  ss.ij) / dx :
                                                 u[0] * (ss.ij  - ss.im1) / dx;
        const PetscScalar UpEnthv = (v[0] < 0) ? v[0] * (ss.jp1 -  ss.ij) / dy :
                                                 v[0] * (ss.ij  - ss.jm1) / dy;
        rhs[k0] -= dtTemp * (UpEnthu + UpEnthv);     // hor. advection
        rhs[k0] += dtTemp * Sigma[0] / ice_rho;      // strain heat
      }
    } else { 
      // there is *grounded* ice AND there is bedrock at interface
      // WARNING: subtle consequences of finite volume argument across interface;
      // in this segment even temperate ice is conductive (same value as cold ice)
      // see BOMBPROOF docs
      // note sure how to proceed here, keep cold value or replace with temperate?
      const PetscScalar rho_ratio  = ice_rho / bed_thermal_rho,
                        c_ratioINV = bed_thermal_c / ice_c;
      L[k0] = - 2.0 * brR;
      if (Enth[k0] > Enth_s[k0]) {
        D[k0] = rho_ratio * (1.0 + 2.0 * iceRcold);
      } else {
        D[k0] = (1.0 + 2.0 * brR) * c_ratioINV + rho_ratio * (1.0 + 2.0 * iceRcold);
      }
      U[k0] = - rho_ratio * 2.0 * iceRcold;
      if (w[0] < 0.0) { // velocity downward: upwind vertical
        const PetscScalar AA = rho_ratio * dtTemp * w[0] / dzav;
        D[k0] -= AA;
        U[k0] += AA;
      }
      if (Enth[k0] > Enth_s[k0]) { // decide cold vs temperate using prev enthalpy
        rhs[k0] = rho_ratio * Enth[0] - 2 * brR * c_ratioINV * Enth_s[k0];
      } else {
        rhs[k0] = (rho_ratio + c_ratioINV) * Enth[0];
      }
      rhs[k0] += 2.0 * dtTemp * Rb / (bed_thermal_rho * dzav); // frictional heat
      if (!isMarginal) {
        planeStar ss;
        ierr = Enth3->getPlaneStarZ(i,j,0.0,&ss);
        const PetscScalar UpEnthu = (u[0] < 0) ? u[0] * (ss.ip1 -  ss.ij) / dx :
                                                 u[0] * (ss.ij  - ss.im1) / dx;
        const PetscScalar UpEnthv = (v[0] < 0) ? v[0] * (ss.jp1 -  ss.ij) / dy :
                                                 v[0] * (ss.ij  - ss.jm1) / dy;
        rhs[k0] -= dtTemp * rho_ratio * (UpEnthu + UpEnthv); // hor. advection
        rhs[k0] += dtTemp * Sigma[0] / bed_thermal_rho;      // strain heat
      }
    }
  }

  // generic ice segment: build k0+1:k0+ks-1 eqns
  for (PetscInt k = 1; k < ks; k++) {
    const PetscScalar AA = nuEQ * w[k],
                      // zero conduction if omega > 0  <==>  E > E_s(p) :
                      R = (Enth[k] > Enth_s[k]) ? iceRtemp : iceRcold;
    if (w[k] >= 0.0) {  // velocity upward
      L[k0+k] = - R - AA * (1.0 - lambda/2.0);
      D[k0+k] = 1.0 + 2.0 * R + AA * (1.0 - lambda);
      U[k0+k] = - R + AA * (lambda/2.0);
    } else {  // velocity downward
      L[k0+k] = - R - AA * (lambda/2.0);
      D[k0+k] = 1.0 + 2.0 * R - AA * (1.0 - lambda);
      U[k0+k] = - R + AA * (1.0 - lambda/2.0);
    }
    rhs[k0+k] = Enth[k];
    if (!isMarginal) {
      planeStar ss;
      ierr = Enth3->getPlaneStarZ(i,j,k * dzEQ,&ss);
      const PetscScalar UpEnthu = (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx :
                                               u[k] * (ss.ij  - ss.im1) / dx;
      const PetscScalar UpEnthv = (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy :
                                               v[k] * (ss.ij  - ss.jm1) / dy;
      rhs[k0+k] += dtTemp * ((Sigma[k] / ice_rho) - UpEnthu - UpEnthv);
    }
  }
      
  // surface b.c.
  if (ks>0) {
    L[k0+ks] = 0.0;
    D[k0+ks] = 1.0;
    // ignore U[k0+ks]
    rhs[k0+ks] = Enth_ks;
  }

  // mark column as done
  indicesValid = false;
  schemeParamsValid = false;
  surfBCsValid = false;
  basalBCsValid = false;

  // solve it; note melting not addressed yet
  return solveTridiagonalSystem(k0+ks+1,x);
}



#if 1

bedrockOnlySystemCtx::bedrockOnlySystemCtx(int my_Mbz)
      : columnSystemCtx(my_Mbz) {
  Mbz = my_Mbz;
  if (Mbz <= 1) {
    PetscErrorPrintf(
       "\n\n\n  Mbz > 1 required in bedrockOnlySystemCtx (constructor error) \n\n");
    endPrintRank();
  }  
  k0 = Mbz - 1;
  // set flags to indicate nothing yet set
  initAllDone = false;
  topBCValid = false;
  basalBCValid = false;
  // set values so we can check if init was called on all
  dtTemp = -1;
  dzbEQ = -1;
  bed_thermal_rho = -1;
  bed_thermal_c   = -1;
  bed_thermal_k   = -1;
  T_b = NULL;
}


PetscErrorCode bedrockOnlySystemCtx::initAllColumns() {
  // check whether each parameter & pointer got set
  if (dtTemp <= 0.0) { SETERRQ(4,"un-initialized dtTemp in bedrockOnlySystemCtx"); }
  if (dzbEQ <= 0.0) { SETERRQ(6,"un-initialized dzbEQ in bedrockOnlySystemCtx"); }
  if (bed_thermal_rho <= 0.0) { SETERRQ(10,"un-initialized bed_thermal_rho in bedrockOnlySystemCtx"); }
  if (bed_thermal_c <= 0.0) { SETERRQ(11,"un-initialized bed_thermal_c_p in bedrockOnlySystemCtx"); }
  if (bed_thermal_k <= 0.0) { SETERRQ(12,"un-initialized bed_thermal_k in bedrockOnlySystemCtx"); }
  if (T_b == NULL) { SETERRQ(14,"un-initialized pointer T_b in bedrockOnlySystemCtx"); }
  // set derived constants
  brK = bed_thermal_k / (bed_thermal_rho * bed_thermal_c);
  brR = brK * dtTemp / PetscSqr(dzbEQ);
  // done
  initAllDone = true;
  return 0;
}


PetscErrorCode bedrockOnlySystemCtx::setTopBoundaryValueThisColumn(PetscScalar my_Ttop) {
  if (!initAllDone) {  SETERRQ(2,
     "setTopBoundaryValueThisColumn() should only be called after initAllColumns() in bedrockOnlySystemCtx"); }
  if (topBCValid) {  SETERRQ(3,
     "setTopBoundaryValueThisColumn() called twice (?) in bedrockOnlySystemCtx"); }
  Ttop = my_Ttop;
  topBCValid = true;
  return 0;
}


PetscErrorCode bedrockOnlySystemCtx::setBasalBoundaryValueThisColumn(PetscScalar my_Ghf) {
  if (!initAllDone) {  SETERRQ(2,
     "setBasalBoundaryValueThisColumn() should only be called after initAllColumns() in bedrockOnlySystemCtx"); }
  if (basalBCValid) {  SETERRQ(3,
     "setBasalBoundaryValueThisColumn() called twice (?) in bedrockOnlySystemCtx"); }
  Ghf = my_Ghf;
  basalBCValid = true;
  return 0;
}


PetscErrorCode bedrockOnlySystemCtx::viewConstants(PetscViewer viewer) {
  PetscErrorCode ierr;
  
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }

  PetscTruth iascii;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for bedrockOnlySystemCtx::viewConstants()\n"); }

  ierr = PetscViewerASCIIPrintf(viewer,"\n\n<<VIEWING bedrockOnlySystemCtx:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  initAllDone = %d\n",initAllDone); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  dtTemp,dzbEQ = %10.3e,%8.2f\n",
                     dtTemp,dzbEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  bed_thermal_rho,bed_thermal_c,bed_thermal_k = %10.3e,%10.3e,%10.3e\n",
                     bed_thermal_rho,bed_thermal_c,bed_thermal_k); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  brK,brR = %10.3e,%10.3e\n",
                     brK,brR); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  Ttop,Ghf = %10.3e,%10.3e\n",
                     Ttop,Ghf); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,">>\n\n"); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode bedrockOnlySystemCtx::solveThisColumn(PetscScalar **x) {
  if (!initAllDone) {  SETERRQ(2,
     "solveThisColumn() should only be called after initAllColumns() in bedrockOnlySystemCtx"); }
  if (!topBCValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setTopBoundaryValueThisColumn() in bedrockOnlySystemCtx"); }
  if (!basalBCValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setBasalBoundaryValueThisColumn() in bedrockOnlySystemCtx"); }

  if (Mbz <= 1) { SETERRQ(1,"Mbz<=1; constructor should have caught; how did I get here?"); }

  // geothermal flux at base; method should gives O(\Delta t,\Delta z^2)
  //   local truncation error; comes from introducing one lower level, etc.
  D[0] = (1.0 + 2.0 * brR);
  U[0] = - 2.0 * brR;  
  // L[0] not allocated anyway
  rhs[0] = T_b[0] + 2.0 * dtTemp * Ghf / (bed_thermal_rho * bed_thermal_c * dzbEQ);

  // generic segments; because Mbz >= 2, k0 = Mbz-1 >= 1
  //   the loop does not go at all in Mbz==2 case
  for (PetscInt k = 1; k < k0; k++) {
    L[k] = -brR;
    D[k] = 1.0 + 2.0 * brR;
    U[k] = -brR;
    rhs[k] = T_b[k];
  }

  // known temperature (Dirichlet condition) at top
  L[k0] = 0.0;
  D[k0] = 1.0;
  // ignore U[k0]; not accessed
  rhs[k0] = Ttop;

  // mark column as done
  topBCValid = false;
  basalBCValid = false;

  // solve it
  return solveTridiagonalSystem(Mbz,x);
}

#endif


// Copyright (C) 2009-2010 Andreas Aschwanden and Ed Bueler
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

#include "pism_const.hh"        // e.g. for MASK_FLOATING
#include "enthalpyConverter.hh"
#include "enthColumnSystem.hh"


enthSystemCtx::enthSystemCtx(int my_Mz, int my_Mbz)
      : columnSystemCtx(my_Mbz + 1 + my_Mz) {  // <- critical: sets size of sys
  Mz = my_Mz;
  Mbz = my_Mbz;

  // set flags to indicate nothing yet set
  initAllDone = false;
  schemeParamsValid = false;
  BCsValid = false;
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
  bed_rho = -1;
  bed_c   = -1;
  bed_k   = -1;
  Enth = NULL;
  Enth_s = NULL;
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
  if (ice_c <= 0.0) { SETERRQ(8,"un-initialized ice_c in enthSystemCtx"); }
  if (ice_k <= 0.0) { SETERRQ(9,"un-initialized ice_k in enthSystemCtx"); }
  if (ice_nu < 0.0) { SETERRQ(21,"un-initialized ice_nu in enthSystemCtx"); } // 0.0 is allowed
  if (bed_rho <= 0.0) { SETERRQ(10,"un-initialized bed_rho in enthSystemCtx"); }
  if (bed_c <= 0.0) { SETERRQ(11,"un-initialized bed_c in enthSystemCtx"); }
  if (bed_k <= 0.0) { SETERRQ(12,"un-initialized bed_k in enthSystemCtx"); }
  if (Enth == NULL) { SETERRQ(13,"un-initialized pointer Enth in enthSystemCtx"); }
  if (Enth_s == NULL) { SETERRQ(20,"un-initialized pointer Enth_s in enthSystemCtx"); }
  if (Tb == NULL) { SETERRQ(14,"un-initialized pointer Tb in enthSystemCtx"); }
  if (u == NULL) { SETERRQ(15,"un-initialized pointer u in enthSystemCtx"); }
  if (v == NULL) { SETERRQ(16,"un-initialized pointer v in enthSystemCtx"); }
  if (w == NULL) { SETERRQ(17,"un-initialized pointer w in enthSystemCtx"); }
  if (Sigma == NULL) { SETERRQ(18,"un-initialized pointer Sigma in enthSystemCtx"); }
  if (Enth3 == NULL) { SETERRQ(19,"un-initialized pointer Enth3 in enthSystemCtx"); }
  // set derived constants
  nuEQ = dtTemp / dzEQ;
  iceK = ice_k / (ice_rho * ice_c);
  iceRcold = iceK * dtTemp / PetscSqr(dzEQ);
  iceRtemp = ice_nu * dtTemp / PetscSqr(dzEQ);
  bedK = bed_k / (bed_rho * bed_c);
  bedR = bedK * dtTemp / PetscSqr(dzbEQ);
  // done
  initAllDone = true;
  return 0;
}


PetscErrorCode enthSystemCtx::setSchemeParamsThisColumn(
                     const bool my_isfloating, bool my_ismarginal,
                     const PetscScalar my_lambda) {
  if (!initAllDone) {  SETERRQ(2,
     "setSchemeParamsThisColumn() should only be called after\n"
     "  initAllColumns() in enthSystemCtx"); }
  if (schemeParamsValid) {  SETERRQ(3,
     "setSchemeParamsThisColumn() called twice (?) in enthSystemCtx"); }
  isfloating = my_isfloating;
  ismarginal = my_ismarginal;
  lambda = my_lambda;
  schemeParamsValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::setBoundaryValuesThisColumn(
                     const PetscScalar my_Enth_surface,
                     const PetscScalar my_Ghf, const PetscScalar my_Rb) {
  if (!initAllDone) {  SETERRQ(2,
     "setBoundaryValuesThisColumn() should only be called after\n"
     "  initAllColumns() in enthSystemCtx"); }
  if (BCsValid) {  SETERRQ(3,
     "setBoundaryValuesThisColumn() called twice (?) in enthSystemCtx"); }
  Enth_ks = my_Enth_surface;
  Ghf = my_Ghf;
  Rb = my_Rb;
  BCsValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::viewConstants(
                     PetscViewer viewer, bool show_col_dependent) {
  PetscErrorCode ierr;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }

  PetscTruth iascii;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for enthSystemCtx::viewConstants()\n"); }
  
  ierr = PetscViewerASCIIPrintf(viewer,"\n<<VIEWING enthSystemCtx:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "initAllDone = %d\n"
                     "for ALL columns:\n",
                     initAllDone); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  dx,dy,dtTemp,dzEQ,dzbEQ = %8.2f,%8.2f,%10.3e,%8.2f,%8.2f\n",
                     dx,dy,dtTemp,dzEQ,dzbEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  ice_rho,ice_c,ice_k = %10.3e,%10.3e,%10.3e\n",
                     ice_rho,ice_c,ice_k); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  bed_rho,bed_c,bed_k = %10.3e,%10.3e,%10.3e\n",
                     bed_rho,bed_c,bed_k); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  nuEQ = %10.3e\n",
                     nuEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  iceK,iceRcold,iceRtemp = %10.3e,%10.3e,%10.3e,\n",
		     iceK,iceRcold,iceRtemp); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  bedK,bedR = %10.3e,%10.3e\n",
		     bedK,bedR); CHKERRQ(ierr);
  if (show_col_dependent) {
    ierr = PetscViewerASCIIPrintf(viewer,
                     "for THIS column:\n"
                     "  schemeParamsValid,BCsValid = (%d,%d)\n",
                     (int)schemeParamsValid,(int)BCsValid); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  i,j,ks = %d,%d,%d\n",
                     i,j,ks); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  isfloating,ismarginal,lambda, = %d,%d,%10.3f\n",
                     (int)isfloating,(int)ismarginal,lambda); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  Enth_ks,Ghf,Rb = %10.3e,%10.3e,%10.3e\n",
                     Enth_ks,Ghf,Rb); CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,">>\n\n"); CHKERRQ(ierr);
  return 0;
}


/*! \brief Solve the tridiagonal system, in a single column, which determines
the new values of the ice enthalpy and the bedrock temperature. */
PetscErrorCode enthSystemCtx::solveThisColumn(PetscScalar **x) {

  if (!initAllDone) {  SETERRQ(2,
     "solveThisColumn() should only be called after\n"
     "  initAllColumns() in enthSystemCtx"); }
  if (!schemeParamsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after\n"
     "  setSchemeParamsThisColumn() in enthSystemCtx"); }
  if (!BCsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after\n"
     "  setBoundaryValuesThisColumn() in enthSystemCtx"); }

  if (Mbz > 1) { // bedrock present: build k=0:Mbz-1 eqns
    // eqn:  - k_b (d Tb / d zb) = G
    // L[0] is not an allocated or used location
    D[0] = (1.0 + 2.0 * bedR);
    U[0] = - 2.0 * bedR;  
    rhs[0] = Tb[0] + 2.0 * dtTemp * Ghf / (bed_rho * bed_c * dzbEQ);

    // k=1:Mbz-2  bedrock only; pure vertical conduction problem
    for (PetscInt k = 1; k < Mbz-1; k++) {
      L[k] = -bedR;
      D[k] = 1.0 + 2.0 * bedR;
      U[k] = -bedR;
      rhs[k] = Tb[k];
    }
    
    // k=Mbz-1 equation says flux delivered to base of ice is determined by
    //         bedrock temperature gradient and heat from sliding
    if (isfloating) {
      L[Mbz-1] = 0.0;
      D[Mbz-1] = 0.0;
      U[Mbz-1] = - 1.0;  // "+X" on right
      rhs[Mbz-1] = - Ghf;  // "+G" on left
    } else {
      L[Mbz-1] = bed_k / dzbEQ;
      D[Mbz-1] = - bed_k / dzbEQ;
      U[Mbz-1] = - 1.0;  // "-X" on left in cold case, "+X" on right in warm case
      rhs[Mbz-1] = - Rb;
    }

  } else { // no bedrock; k=0 equation *is* k=Mbz-1 equation
    // k=0 equation says flux delivered to base of ice is determined by
    //     bedrock temperature gradient and heat from sliding
    if (isfloating) {
      // L[0] = 0.0;  // not allocated
      D[0] = 0.0;
      U[0] = - 1.0;  // "+X" on right
      rhs[0] = - Ghf;  // "+G" on left
    } else {
      // L[0] = 0.0;  // not allocated
      D[0] = 0.0;
      U[0] = - 1.0;  // "-X" on left in cold case, "+X" on right in warm case
      rhs[0] = - Ghf - Rb;  // "+G" on left
    }    
  }

  // k=Mbz equation says top-of-bedrock temperature equals base of ice temperature
  D[Mbz] = 0.0;       // unknown X does not play role in this equation; CAUTION
  if (isfloating) {
    // eqn:  Tb[Mbz-1] + 0 X + 0 Enth[0] = c_i^{-1} Enth_s[0]    (=T_m(p))
    L[Mbz] = 1.0;
    U[Mbz] = 0.0;  
    rhs[Mbz] = (1.0 / ice_c) * Enth_s[0];
  } else {
    if (Enth[0] < Enth_s[0]) {  // cold base
      // eqn:  Tb[Mbz-1] + 0 X - c_i^{-1} Enth[0] = 0
      L[Mbz] = 1.0;
      U[Mbz] = - 1.0 / ice_c;  
      rhs[Mbz] = 0.0;
    } else { // warm base
      // eqn:  Tb[Mbz-1] + 0 X + 0 Enth[0] = c_i^{-1} Enth_s[0]   (=T_m(p))
      L[Mbz] = 1.0;
      U[Mbz] = 0.0;  
      rhs[Mbz] = (1.0 / ice_c) * Enth_s[0];
    }
  }

  // k=Mbz+1 eqn: usually says heat flux into ice is zero (melting/warm base
  //   case) or given by X (see above)
  if (ks == 0) {
    // essentially no ice; set Enth[0] to air value if grounded or to pressure-
    //   melting if sub-ice shelf
    L[Mbz+1] = 0.0; 
    D[Mbz+1] = 1.0;
    U[Mbz+1] = 0.0;
    rhs[Mbz+1] = (isfloating) ? (1.0 / ice_c) * Enth_s[0] : Enth_ks;
  } else { 
    // there is more than one grid space of ice
    if (isfloating) {
      // apply Neumann condition to base of column of ice in an ice shelf
      // eqn:   - E[0] + E[1] = 0
      L[Mbz+1] = 0.0; 
      D[Mbz+1] = - 1.0;
      U[Mbz+1] = 1.0;
      rhs[Mbz+1] = 0.0;
    } else {
      if (Enth[0] < Enth_s[0]) {  // cold base
        // eqn:  X - (k_i/(c_i*dz)) E[0] + (k_i/(c_i*dz)) E[0] = 0
        const PetscScalar A = ice_k / (ice_c * dzEQ);
        L[Mbz+1] = 1.0;
        D[Mbz+1] = - A;
        U[Mbz+1] = A;
        rhs[Mbz+1] =  0.0;
      } else { // warm base
        // apply Neumann condition to base of column of ice because melt eats heat flux
        // eqn:   - E[0] + E[1] = 0
        L[Mbz+1] = 0.0; 
        D[Mbz+1] = - 1.0;
        U[Mbz+1] = 1.0;
        rhs[Mbz+1] = 0.0;
      }
    }
  }

  // generic ice segment in Mbz+1+k location (if any; only runs if ks >= 2)
  for (PetscInt k = 1; k < ks; k++) {
    const PetscScalar AA = nuEQ * w[k],
                      R = (Enth[k] > Enth_s[k]) ? iceRtemp : iceRcold;
    if (w[k] >= 0.0) {  // velocity upward
      L[Mbz+1+k] = - R - AA * (1.0 - lambda/2.0);
      D[Mbz+1+k] = 1.0 + 2.0 * R + AA * (1.0 - lambda);
      U[Mbz+1+k] = - R + AA * (lambda/2.0);
    } else {            // velocity downward
      L[Mbz+1+k] = - R - AA * (lambda/2.0);
      D[Mbz+1+k] = 1.0 + 2.0 * R - AA * (1.0 - lambda);
      U[Mbz+1+k] = - R + AA * (1.0 - lambda/2.0);
    }
    rhs[Mbz+1+k] = Enth[k];
    if (!ismarginal) {
      planeStar ss;
      Enth3->getPlaneStarZ(i,j,k * dzEQ,&ss);
      const PetscScalar UpEnthu = (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx :
                                               u[k] * (ss.ij  - ss.im1) / dx;
      const PetscScalar UpEnthv = (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy :
                                               v[k] * (ss.ij  - ss.jm1) / dy;
      rhs[Mbz+1+k] += dtTemp * ((Sigma[k] / ice_rho) - UpEnthu - UpEnthv);
    }
  }

  // air above
  for (PetscInt k = ks; k < Mz; k++) {
    L[Mbz+1+k] = 0.0;
    D[Mbz+1+k] = 1.0;
    if (k < Mz) U[Mbz+1+k] = 0.0;
    rhs[Mbz+1+k] = Enth_ks;
  }

  // mark column as done
  schemeParamsValid = false;
  BCsValid = false;

  // solve it; note drainage is not addressed yet and post-processing may occur
  return solveTridiagonalSystem(Mbz + 1 + Mz, x);
}


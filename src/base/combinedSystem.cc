// Copyright (C) 2009-2010 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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

#include "enthalpyConverter.hh"
#include "combinedSystem.hh"
#include <gsl/gsl_math.h>


combinedSystemCtx::combinedSystemCtx(
          const NCConfigVariable &config, IceModelVec3 &my_Enth3, 
          int my_Mz, int my_Mbz)
  : columnSystemCtx(my_Mbz + my_Mz - 1) {  // <- critical: sets size of sys

  Mz = my_Mz;
  Mbz = my_Mbz;

  // set values so we can check if init was called
  dtTemp   = -1.0;
  dx       = -1.0;
  dy       = -1.0;
  dzEQ     = -1.0;
  dzbEQ    = -1.0;
  bedR     = -1.0;
  nuEQ     = -1.0;
  iceRcold = -1.0;
  iceRtemp = -1.0;
  lambda   = -1.0;
  Enth_ks  = -1.0;
  Ghf      = -1.0;
  Fb       = -1.0;

  ice_rho = config.get("ice_density");
  ice_c   = config.get("ice_specific_heat_capacity");
  ice_k   = config.get("ice_thermal_conductivity");
  // diffusion in temperate ice; ice at base is cold in combinedSystem, but
  //    ice further up ice column could be temperate
  ice_nu  = config.get("enthalpy_temperate_diffusivity");
  bed_rho = config.get("bedrock_thermal_density");                // kg m-3
  bed_c   = config.get("bedrock_thermal_specific_heat_capacity"); // J kg-1 K-1
  bed_k   = config.get("bedrock_thermal_conductivity");           // W m-1 K-1

  iceK = ice_k / (ice_rho * ice_c);
  bedK = bed_k / (bed_rho * bed_c);

  u      = new PetscScalar[Mz];
  v      = new PetscScalar[Mz];
  w      = new PetscScalar[Mz];
  Sigma  = new PetscScalar[Mz];
  Enth_s = new PetscScalar[Mz];  // enthalpy of pressure-melting-point
  Enth   = new PetscScalar[Mz];
  Tb     = new PetscScalar[Mbz];  // bedrock temps at prev step

  Enth3 = &my_Enth3;  // points to IceModelVec3

  // invalidate column inputs
  u[0]      = GSL_NAN;
  v[0]      = GSL_NAN;
  w[0]      = GSL_NAN;
  Sigma[0]  = GSL_NAN;
  Enth_s[0] = GSL_NAN;
  Enth[0]   = GSL_NAN;
  Tb[0]     = GSL_NAN;
}


combinedSystemCtx::~combinedSystemCtx() {
  delete [] u;
  delete [] v;
  delete [] w;
  delete [] Sigma;
  delete [] Enth_s;
  delete [] Enth;
  delete [] Tb;
}


PetscErrorCode combinedSystemCtx::initAllColumns(
      const PetscScalar my_dx, const PetscScalar my_dy, 
      const PetscScalar my_dtTemp,
      const PetscScalar my_dzEQ, const PetscScalar my_dzbEQ) {
  dx       = my_dx;
  dy       = my_dy;
  dtTemp   = my_dtTemp;
  dzEQ     = my_dzEQ;
  dzbEQ    = my_dzbEQ;
  nuEQ     = dtTemp / dzEQ;
  iceRcold = iceK * dtTemp / PetscSqr(dzEQ);
  iceRtemp = ice_nu * dtTemp / PetscSqr(dzEQ);
  bedR     = bedK * dtTemp / PetscSqr(dzbEQ);
  return 0;
}


PetscErrorCode combinedSystemCtx::setSchemeParamsThisColumn(
                     bool my_ismarginal, const PetscScalar my_lambda) {
  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (bedR < 0.0)) {
    SETERRQ(2, "setSchemeParamsThisColumn() should only be called after\n"
               "  initAllColumns() in combinedSystemCtx"); }
  if (lambda >= 0.0) {
    SETERRQ(3, "setSchemeParamsThisColumn() called twice (?) in combinedSystemCtx"); }
  ismarginal = my_ismarginal;
  lambda = my_lambda;
  return 0;
}


PetscErrorCode combinedSystemCtx::setBoundaryValuesThisColumn(
           const PetscScalar my_Enth_surface, const PetscScalar my_Ghf,
          const PetscScalar my_Fb) {
  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (bedR < 0.0)) {
    SETERRQ(2,"setBoundaryValuesThisColumn() should only be called after\n"
              "  initAllColumns() in combinedSystemCtx"); }
  if ((Enth_ks >= 0.0) || (Ghf >= 0.0) || (Fb >= 0.0)) {
    SETERRQ(3,"setBoundaryValuesThisColumn() called twice (?) in combinedSystemCtx"); }
  Enth_ks = my_Enth_surface;
  Ghf     = my_Ghf;
  Fb      = my_Fb;
  return 0;
}


PetscErrorCode combinedSystemCtx::viewConstants(
                     PetscViewer viewer, bool show_col_dependent) {
  PetscErrorCode ierr;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }

  PetscTruth iascii;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for combinedSystemCtx::viewConstants()\n"); }
  
  ierr = PetscViewerASCIIPrintf(viewer,
                   "\n<<VIEWING combinedSystemCtx:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "for ALL columns:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  dx,dy,dtTemp,dzEQ,dzbEQ = %8.2f,%8.2f,%10.3e,%8.2f,%8.2f\n",
                     dx,dy,dtTemp,dzEQ,dzbEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  ice_rho,ice_c,ice_k,ice_nu = %10.3e,%10.3e,%10.3e,%10.3e\n",
                     ice_rho,ice_c,ice_k,ice_nu); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  bed_rho,bed_c,bed_k = %10.3e,%10.3e,%10.3e\n",
                     bed_rho,bed_c,bed_k); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  nuEQ = %10.3e\n",
                     nuEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  iceK,iceRcold,iceRtemp,bedK,bedR = %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n",
		     iceK,iceRcold,iceRtemp,bedK,bedR); CHKERRQ(ierr);
  if (show_col_dependent) {
    ierr = PetscViewerASCIIPrintf(viewer,
                     "for THIS column:\n"
                     "  i,j,ks = %d,%d,%d\n",
                     i,j,ks); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  ismarginal,lambda = %d,%10.3f\n",
                     (int)ismarginal,lambda); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  Enth_ks,Ghf,Fb = %10.3e,%10.3e,%10.3e\n",
                     Enth_ks,Ghf,Fb); CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,
                     ">>\n\n"); CHKERRQ(ierr);
  return 0;
}


/*! \brief Solve the tridiagonal system, in a single column, which determines
the new values of the ice enthalpy and the bedrock temperature. */
PetscErrorCode combinedSystemCtx::solveThisColumn(PetscScalar **x) {

  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (bedR < 0.0)) {
    SETERRQ(2, "solveThisColumn() should only be called after\n"
               "  initAllColumns() in combinedSystemCtx"); }
  if (lambda < 0.0) {
    SETERRQ(3, "solveThisColumn() should only be called after\n"
               "  setSchemeParamsThisColumn() in combinedSystemCtx"); }
  if ((Enth_ks < 0.0) || (Ghf < 0.0) || (Fb < 0.0)) {
    SETERRQ(4, "solveThisColumn() should only be called after\n"
               "  setBoundaryValuesThisColumn() in combinedSystemCtx"); }
  if (Mbz <= 1) {
    SETERRQ(5, "solveThisColumn() should only be called\n"
               "  if Mbz > 1 (in combinedSystemCtx)"); }

  if (Mbz == 2) {
    SETERRQ(6, "solveThisColumn() NOT IMPLEMENTED FOR\n"
               "  SPECIAL CASE Mbz == 2,  in combinedSystemCtx"); }

  if (gsl_isnan(u[0])) {
    SETERRQ(60, "solveThisColumn() called with invalid u[] in\n"
               "  combinedSystemCtx"); }
  if (gsl_isnan(v[0])) {
    SETERRQ(61, "solveThisColumn() called with invalid v[] in\n"
               "  combinedSystemCtx"); }
  if (gsl_isnan(w[0])) {
    SETERRQ(62, "solveThisColumn() called with invalid w[] in\n"
               "  combinedSystemCtx"); }
  if (gsl_isnan(Sigma[0])) {
    SETERRQ(63, "solveThisColumn() called with invalid Sigma[] in\n"
               "  combinedSystemCtx"); }
  if (gsl_isnan(Enth_s[0])) {
    SETERRQ(64, "solveThisColumn() called with invalid Enth_s[] in\n"
               "  combinedSystemCtx"); }
  if (gsl_isnan(Enth[0])) {
    SETERRQ(65, "solveThisColumn() called with invalid Enth[] in\n"
               "  combinedSystemCtx"); }
  if (gsl_isnan(Tb[0])) {
    SETERRQ(66, "solveThisColumn() called with invalid Tb[] in\n"
               "  combinedSystemCtx"); }

  // eqn:  - k_b (d Tb / d zb) = G + (heat equation); uses "add a point
  //       past the end" trick (Morton & Mayers)
  // L[0] is not allocated
  D[0] = (1.0 + 2.0 * bedR);
  U[0] = - 2.0 * bedR;  
  rhs[0] = Tb[0] + 2.0 * dtTemp * Ghf / (bed_rho * bed_c * dzbEQ);

  // k=1:Mbz-3  bedrock only; pure vertical conduction problem
  for (PetscInt k = 1; k <= Mbz-3; k++) {
    L[k] = -bedR;
    D[k] = 1.0 + 2.0 * bedR;
    U[k] = -bedR;
    rhs[k] = Tb[k];
  }

  // because k=Mbz-1 level is actually E[0], the k=Mbz-2 equation has a special U
  L[Mbz-2] = -bedR;
  D[Mbz-2] = 1.0 + 2.0 * bedR;
  U[Mbz-2] = -bedR / ice_c;
  rhs[Mbz-2] = Tb[Mbz-2];

  // k=Mbz-1 equation says flux is continuous at bedrock-ice transition, up to 
  //         heat from sliding; eqn is  L Tb[Mbz-2] + D E[0] + U E[1] = rhs
  //         so different scaling for L versus D and U
  L[Mbz-1] = - bed_k / dzbEQ; 
  D[Mbz-1] = (1.0 / ice_c) * ( (bed_k / dzbEQ) + (ice_k / dzEQ) );
  U[Mbz-1] = - (1.0 / ice_c) * (ice_k / dzEQ);
  rhs[Mbz-1] = Fb;

  // generic ice segment in k location (if any; only runs if ks >= 2)
  for (PetscInt k = 1; k < ks; k++) {
    const PetscScalar AA = nuEQ * w[k],
                      R = (Enth[k] > Enth_s[k]) ? iceRtemp : iceRcold;
    if (w[k] >= 0.0) {  // velocity upward
      L[Mbz-1+k] = - R - AA * (1.0 - lambda/2.0);
      D[Mbz-1+k] = 1.0 + 2.0 * R + AA * (1.0 - lambda);
      U[Mbz-1+k] = - R + AA * (lambda/2.0);
    } else {            // velocity downward
      L[Mbz-1+k] = - R - AA * (lambda/2.0);
      D[Mbz-1+k] = 1.0 + 2.0 * R - AA * (1.0 - lambda);
      U[Mbz-1+k] = - R + AA * (1.0 - lambda/2.0);
    }
    rhs[Mbz-1+k] = Enth[k];
    if (!ismarginal) {
      planeStar ss;
      Enth3->getPlaneStar_fine(i,j,k,&ss);
      const PetscScalar UpEnthu = (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx :
                                               u[k] * (ss.ij  - ss.im1) / dx;
      const PetscScalar UpEnthv = (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy :
                                               v[k] * (ss.ij  - ss.jm1) / dy;
      rhs[Mbz-1+k] += dtTemp * ((Sigma[k] / ice_rho) - UpEnthu - UpEnthv);
    }
  }

  // air above
  for (PetscInt k = ks; k < Mz; k++) {
    L[Mbz-1+k] = 0.0;
    D[Mbz-1+k] = 1.0;
    if (k < Mz-1) U[Mbz-1+k] = 0.0;
    rhs[Mbz-1+k] = Enth_ks;
  }

  // solve it; note drainage is not addressed yet and post-processing may occur
  PetscErrorCode retval = solveTridiagonalSystem(Mbz-1+Mz, x);

  if (!retval) { // mark column as done by making scheme params and b.c.s invalid
    lambda  = -1.0;
    Enth_ks = -1.0;
    Ghf     = -1.0;
    Fb      = -1.0;
    // invalidate column inputs
    u[0]      = GSL_NAN;
    v[0]      = GSL_NAN;
    w[0]      = GSL_NAN;
    Sigma[0]  = GSL_NAN;
    Enth_s[0] = GSL_NAN;
    Enth[0]   = GSL_NAN;
    Tb[0]     = GSL_NAN;
  }
  return retval;
}



// Copyright (C) 2009-2011 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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
#include "enthSystem.hh"
#include <gsl/gsl_math.h>


enthSystemCtx::enthSystemCtx(
    const NCConfigVariable &config, IceModelVec3 &my_Enth3, int my_Mz)
      : columnSystemCtx(my_Mz) {  // <- critical: sets size of sys
  Mz = my_Mz;

  // set some values so we can check if init was called
  nuEQ     = -1.0;
  iceRcold = -1.0;
  iceRtemp = -1.0;
  lambda   = -1.0;
  a0 = GSL_NAN;
  a1 = GSL_NAN;
  b  = GSL_NAN;

  ice_rho = config.get("ice_density");
  ice_c   = config.get("ice_specific_heat_capacity");
  ice_k   = config.get("ice_thermal_conductivity");

  ice_K   = ice_k / ice_c;
  ice_K0  = ice_K * config.get("enthalpy_temperate_conductivity_ratio");

  u      = new PetscScalar[Mz];
  v      = new PetscScalar[Mz];
  w      = new PetscScalar[Mz];
  Sigma  = new PetscScalar[Mz];
  Enth   = new PetscScalar[Mz];
  Enth_s = new PetscScalar[Mz];  // enthalpy of pressure-melting-point

  Enth3 = &my_Enth3;  // points to IceModelVec3

  // invalidate column inputs
  u[0]      = GSL_NAN;
  v[0]      = GSL_NAN;
  w[0]      = GSL_NAN;
  Sigma[0]  = GSL_NAN;
  Enth_s[0] = GSL_NAN;
  Enth[0]   = GSL_NAN;
}


enthSystemCtx::~enthSystemCtx() {
  delete [] u;
  delete [] v;
  delete [] w;
  delete [] Sigma;
  delete [] Enth_s;
  delete [] Enth;
}


PetscErrorCode enthSystemCtx::initAllColumns(
      const PetscScalar my_dx, const PetscScalar my_dy, 
      const PetscScalar my_dtTemp, const PetscScalar my_dzEQ) {
  dx     = my_dx;
  dy     = my_dy;
  dtTemp = my_dtTemp;
  dzEQ   = my_dzEQ;
  nuEQ     = dtTemp / dzEQ;
  iceRcold = (ice_K / ice_rho) * dtTemp / PetscSqr(dzEQ);
  iceRtemp = (ice_K0 / ice_rho) * dtTemp / PetscSqr(dzEQ);
  return 0;
}


PetscErrorCode enthSystemCtx::setSchemeParamsThisColumn(
                     bool my_ismarginal, const PetscScalar my_lambda) {
#ifdef PISM_DEBUG
  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (iceRtemp < 0.0)) {  SETERRQ(2,
     "setSchemeParamsThisColumn() should only be called after\n"
     "  initAllColumns() in enthSystemCtx"); }
  if (lambda >= 0.0) {  SETERRQ(3,
     "setSchemeParamsThisColumn() called twice (?) in enthSystemCtx"); }
#endif
  ismarginal = my_ismarginal;
  lambda = my_lambda;
  return 0;
}


PetscErrorCode enthSystemCtx::setBoundaryValuesThisColumn(
           const PetscScalar my_Enth_surface) {
#ifdef PISM_DEBUG
  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (iceRtemp < 0.0)) {  SETERRQ(2,
     "setBoundaryValuesThisColumn() should only be called after\n"
     "  initAllColumns() in enthSystemCtx"); }
#endif
  Enth_ks = my_Enth_surface;
  return 0;
}


//! Set coefficients in equation \f$a_0 E_0^{n+1} + a_1 E_1^{n+1} = b\f$ at base of ice.
PetscErrorCode enthSystemCtx::setLevel0EqnThisColumn(
      const PetscScalar my_a0, const PetscScalar my_a1, const PetscScalar my_b) {
#ifdef PISM_DEBUG
  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (iceRtemp < 0.0)) {
    SETERRQ(1, "setLevel0EqnThisColumn() should only be called after\n"
               "  initAllColumns() in enthSystemCtx");
  }
  if ((!gsl_isnan(a0)) || (!gsl_isnan(a1)) || (!gsl_isnan(b))) {
    SETERRQ(2, "setLevel0EqnThisColumn() called twice ? in enthSystemCtx");
  }
#endif
  a0 = my_a0;
  a1 = my_a1;
  b  = my_b;
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
  
  ierr = PetscViewerASCIIPrintf(viewer,
                   "\n<<VIEWING enthSystemCtx:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "for ALL columns:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  dx,dy,dtTemp,dzEQ = %8.2f,%8.2f,%10.3e,%8.2f\n",
                     dx,dy,dtTemp,dzEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  ice_rho,ice_c,ice_k,ice_K,ice_K0 = %10.3e,%10.3e,%10.3e,%10.3e,%10.3e\n",
                     ice_rho,ice_c,ice_k,ice_K,ice_K0); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  nuEQ = %10.3e\n",
                     nuEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  iceRcold,iceRtemp = %10.3e,%10.3e,\n",
		     iceRcold,iceRtemp); CHKERRQ(ierr);
  if (show_col_dependent) {
    ierr = PetscViewerASCIIPrintf(viewer,
                     "for THIS column:\n"
                     "  i,j,ks = %d,%d,%d\n",
                     i,j,ks); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  ismarginal,lambda = %d,%10.3f\n",
                     (int)ismarginal,lambda); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  Enth_ks = %10.3e\n",
                     Enth_ks); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  a0,a1,b = %10.3e,%10.3e\n",
                     a0,a1,b); CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,
                     ">>\n\n"); CHKERRQ(ierr);
  return 0;
}


/*! \brief Solve the tridiagonal system, in a single column, which determines
the new values of the ice enthalpy. */
PetscErrorCode enthSystemCtx::solveThisColumn(PetscScalar **x, PetscErrorCode &pivoterrorindex) {
#ifdef PISM_DEBUG
  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (iceRtemp < 0.0)) {
    SETERRQ(2, "solveThisColumn() should only be called after\n"
               "  initAllColumns() in enthSystemCtx"); }
  if (lambda < 0.0) {
    SETERRQ(3, "solveThisColumn() should only be called after\n"
               "  setSchemeParamsThisColumn() in enthSystemCtx"); }
  if ((gsl_isnan(a0)) || (gsl_isnan(a1)) || (gsl_isnan(b))) {
    SETERRQ(5, "solveThisColumn() should only be called after\n"
               "  setLevel0EqnThisColumn() in enthSystemCtx"); }
  if (gsl_isnan(u[0])) {
    SETERRQ(60, "solveThisColumn() called with invalid u[] in\n"
               "  enthSystemCtx"); }
  if (gsl_isnan(v[0])) {
    SETERRQ(61, "solveThisColumn() called with invalid v[] in\n"
               "  enthSystemCtx"); }
  if (gsl_isnan(w[0])) {
    SETERRQ(62, "solveThisColumn() called with invalid w[] in\n"
               "  enthSystemCtx"); }
  if (gsl_isnan(Sigma[0])) {
    SETERRQ(63, "solveThisColumn() called with invalid Sigma[] in\n"
               "  enthSystemCtx"); }
  if (gsl_isnan(Enth_s[0])) {
    SETERRQ(64, "solveThisColumn() called with invalid Enth_s[] in\n"
               "  enthSystemCtx"); }
  if (gsl_isnan(Enth[0])) {
    SETERRQ(65, "solveThisColumn() called with invalid Enth[] in\n"
               "  enthSystemCtx"); }
#endif
  // k=0 equation must be supplied to avoid making decisions here
  // L[0] = 0.0;  // not allocated
  D[0] = a0;
  U[0] = a1;
  rhs[0] = b;

  // generic ice segment in k location (if any; only runs if ks >= 2)
  for (PetscInt k = 1; k < ks; k++) {
    const PetscScalar
        Rl = (Enth[k-1] < Enth_s[k-1]) ? iceRcold : iceRtemp,
        Rc = (Enth[k] < Enth_s[k]) ? iceRcold : iceRtemp,
        Rr = (Enth[k+1] < Enth_s[k+1]) ? iceRcold : iceRtemp,
        Rminus = 0.5 * (Rl + Rc),
        Rplus  = 0.5 * (Rc + Rr);
    L[k] = - Rminus;
    D[k] = 1.0 + Rminus + Rplus;
    U[k] = - Rplus;
    const PetscScalar AA = nuEQ * w[k];
    if (w[k] >= 0.0) {  // velocity upward
      L[k] -= AA * (1.0 - lambda/2.0);
      D[k] += AA * (1.0 - lambda);
      U[k] += AA * (lambda/2.0);
    } else {            // velocity downward
      L[k] -= AA * (lambda/2.0);
      D[k] -= AA * (1.0 - lambda);
      U[k] += AA * (1.0 - lambda/2.0);
    }
    rhs[k] = Enth[k];
    if (!ismarginal) {
      planeStar ss;
      Enth3->getPlaneStar_fine(i,j,k,&ss);
      const PetscScalar UpEnthu = (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx :
                                               u[k] * (ss.ij  - ss.im1) / dx;
      const PetscScalar UpEnthv = (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy :
                                               v[k] * (ss.ij  - ss.jm1) / dy;
      rhs[k] += dtTemp * ((Sigma[k] / ice_rho) - UpEnthu - UpEnthv);
    }
  }

  // air above
  for (PetscInt k = ks; k < Mz; k++) {
    if (k > 0) L[k] = 0.0;
    D[k] = 1.0;
    if (k < Mz-1) U[k] = 0.0;
    rhs[k] = Enth_ks;
  }

  // solve it; note drainage is not addressed yet and post-processing may occur
  pivoterrorindex = solveTridiagonalSystem(Mz, x);

#ifdef PISM_DEBUG
  if (pivoterrorindex == 0) {
    // if success, mark column as done by making scheme params and b.c.s invalid
    lambda  = -1.0;
    a0 = GSL_NAN;
    a1 = GSL_NAN;
    b  = GSL_NAN;
    // invalidate column inputs
    u[0]      = GSL_NAN;
    v[0]      = GSL_NAN;
    w[0]      = GSL_NAN;
    Sigma[0]  = GSL_NAN;
    Enth_s[0] = GSL_NAN;
    Enth[0]   = GSL_NAN;
  }
#endif
  return 0;
}


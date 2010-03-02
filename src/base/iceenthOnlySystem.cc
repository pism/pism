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

#include "enthalpyConverter.hh"
#include "iceenthOnlySystem.hh"
#include <gsl/gsl_math.h>


iceenthOnlySystemCtx::iceenthOnlySystemCtx(
    const NCConfigVariable &config, IceModelVec3 &my_Enth3, int my_Mz)
      : columnSystemCtx(my_Mz) {  // <- critical: sets size of sys
  Mz = my_Mz;

  // set values so we can check if init was called
  nuEQ     = -1.0;
  iceRcold = -1.0;
  iceRtemp = -1.0;
  lambda   = -1.0;
  Enth_ks  = -1.0;
  a0 = GSL_NAN;
  a1 = GSL_NAN;
  b  = GSL_NAN;

  ice_rho = config.get("ice_density");
  ice_c   = config.get("ice_specific_heat_capacity");
  ice_k   = config.get("ice_thermal_conductivity");
  ice_nu  = config.get("enthalpy_temperate_diffusivity"); // diffusion in temperate ice

  iceK = ice_k / (ice_rho * ice_c);

  u      = new PetscScalar[Mz];
  v      = new PetscScalar[Mz];
  w      = new PetscScalar[Mz];
  Sigma  = new PetscScalar[Mz];
  Enth_s = new PetscScalar[Mz];  // enthalpy of pressure-melting-point
  Enth   = new PetscScalar[Mz];

  Enth3 = &my_Enth3;  // points to IceModelVec3
}


iceenthOnlySystemCtx::~iceenthOnlySystemCtx() {
  delete [] u;
  delete [] v;
  delete [] w;
  delete [] Sigma;
  delete [] Enth_s;
  delete [] Enth;
}


PetscErrorCode iceenthOnlySystemCtx::initAllColumns(
      const PetscScalar my_dx, const PetscScalar my_dy, 
      const PetscScalar my_dtTemp, const PetscScalar my_dzEQ) {
  dx     = my_dx;
  dy     = my_dy;
  dtTemp = my_dtTemp;
  dzEQ   = my_dzEQ;
  nuEQ     = dtTemp / dzEQ;
  iceRcold = iceK * dtTemp / PetscSqr(dzEQ);
  iceRtemp = ice_nu * dtTemp / PetscSqr(dzEQ);
  return 0;
}


PetscErrorCode iceenthOnlySystemCtx::setSchemeParamsThisColumn(
                     bool my_ismarginal, const PetscScalar my_lambda) {
  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (iceRtemp < 0.0)) {  SETERRQ(2,
     "setSchemeParamsThisColumn() should only be called after\n"
     "  initAllColumns() in iceenthOnlySystemCtx"); }
  if (lambda >= 0.0) {  SETERRQ(3,
     "setSchemeParamsThisColumn() called twice (?) in iceenthOnlySystemCtx"); }
  ismarginal = my_ismarginal;
  lambda = my_lambda;
  return 0;
}


PetscErrorCode iceenthOnlySystemCtx::setBoundaryValuesThisColumn(
           const PetscScalar my_Enth_surface) {
  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (iceRtemp < 0.0)) {  SETERRQ(2,
     "setBoundaryValuesThisColumn() should only be called after\n"
     "  initAllColumns() in iceenthOnlySystemCtx"); }
  if (Enth_ks >= 0.0) {  SETERRQ(3,
     "setBoundaryValuesThisColumn() called twice (?) in iceenthOnlySystemCtx"); }
  Enth_ks = my_Enth_surface;
  return 0;
}


//! Set coefficients in equation \f$a_0 E_0^{n+1} + a_1 E_1^{n+1} = b\f$ at base of ice.
PetscErrorCode iceenthOnlySystemCtx::setLevel0EqnThisColumn(
      const PetscScalar my_a0, const PetscScalar my_a1, const PetscScalar my_b) {
  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (iceRtemp < 0.0)) {
    SETERRQ(1, "setLevel0EqnThisColumn() should only be called after\n"
               "  initAllColumns() in iceenthOnlySystemCtx");
  }
  if ((!gsl_isnan(a0)) || (!gsl_isnan(a1)) || (!gsl_isnan(b))) {
    SETERRQ(2, "setLevel0EqnThisColumn() called twice ? in iceenthOnlySystemCtx");
  }
  a0 = my_a0;
  a1 = my_a1;
  b  = my_b;
  return 0;
}


PetscErrorCode iceenthOnlySystemCtx::viewConstants(
                     PetscViewer viewer, bool show_col_dependent) {
  PetscErrorCode ierr;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }

  PetscTruth iascii;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for iceenthOnlySystemCtx::viewConstants()\n"); }
  
  ierr = PetscViewerASCIIPrintf(viewer,
                   "\n<<VIEWING iceenthOnlySystemCtx:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "for ALL columns:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  dx,dy,dtTemp,dzEQ = %8.2f,%8.2f,%10.3e,%8.2f\n",
                     dx,dy,dtTemp,dzEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  ice_rho,ice_c,ice_k,ice_nu = %10.3e,%10.3e,%10.3e,%10.3e\n",
                     ice_rho,ice_c,ice_k,ice_nu); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  nuEQ = %10.3e\n",
                     nuEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  iceK,iceRcold,iceRtemp = %10.3e,%10.3e,%10.3e,\n",
		     iceK,iceRcold,iceRtemp); CHKERRQ(ierr);
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
PetscErrorCode iceenthOnlySystemCtx::solveThisColumn(PetscScalar **x) {

  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (iceRtemp < 0.0)) {
    SETERRQ(2, "solveThisColumn() should only be called after\n"
               "  initAllColumns() in iceenthOnlySystemCtx"); }
  if (lambda < 0.0) {
    SETERRQ(3, "solveThisColumn() should only be called after\n"
               "  setSchemeParamsThisColumn() in iceenthOnlySystemCtx"); }
  if (Enth_ks < 0.0) {
    SETERRQ(4, "solveThisColumn() should only be called after\n"
               "  setBoundaryValuesThisColumn() in iceenthOnlySystemCtx"); }
  if ((gsl_isnan(a0)) || (gsl_isnan(a1)) || (gsl_isnan(b))) {
    SETERRQ(5, "solveThisColumn() should only be called after\n"
               "  setLevel0EqnThisColumn() in iceenthOnlySystemCtx"); }

  // k=0 equation must be supplied to avoid making decisions here
  // L[0] = 0.0;  // not allocated
  D[0] = a0;
  U[0] = a1;
  rhs[0] = b;

  // generic ice segment in k location (if any; only runs if ks >= 2)
  for (PetscInt k = 1; k < ks; k++) {
    const PetscScalar AA = nuEQ * w[k],
                      R = (Enth[k] > Enth_s[k]) ? iceRtemp : iceRcold;
    if (w[k] >= 0.0) {  // velocity upward
      L[k] = - R - AA * (1.0 - lambda/2.0);
      D[k] = 1.0 + 2.0 * R + AA * (1.0 - lambda);
      U[k] = - R + AA * (lambda/2.0);
    } else {            // velocity downward
      L[k] = - R - AA * (lambda/2.0);
      D[k] = 1.0 + 2.0 * R - AA * (1.0 - lambda);
      U[k] = - R + AA * (1.0 - lambda/2.0);
    }
    rhs[k] = Enth[k];
    if (!ismarginal) {
      planeStar ss;
      Enth3->getPlaneStarZ(i,j,k * dzEQ,&ss);
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
  PetscErrorCode retval = solveTridiagonalSystem(Mz, x);

  if (!retval) {
    // if success, mark column as done by making scheme params and b.c.s invalid
    lambda  = -1.0;
    Enth_ks = -1.0;
    a0 = GSL_NAN;
    a1 = GSL_NAN;
    b  = GSL_NAN;
  }
  return retval;
}


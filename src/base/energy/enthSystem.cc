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
    const NCConfigVariable &config, IceModelVec3 &my_Enth3, int my_Mz, string my_prefix)
      : columnSystemCtx(my_Mz, my_prefix) {  // <- critical: sets size of sys
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


PetscErrorCode enthSystemCtx::checkReadyToSolve() {
  if ((nuEQ < 0.0) || (iceRcold < 0.0) || (iceRtemp < 0.0)) {
    SETERRQ(2,  "not ready to solve: need initAllColumns() in enthSystemCtx"); }
  if (lambda < 0.0) {
    SETERRQ(3,  "not ready to solve: need setSchemeParamsThisColumn() in enthSystemCtx"); }
  if (gsl_isnan(u[0])) {
    SETERRQ(60, "not ready to solve: invalid u[] in enthSystemCtx"); }
  if (gsl_isnan(v[0])) {
    SETERRQ(61, "not ready to solve: invalid v[] in enthSystemCtx"); }
  if (gsl_isnan(w[0])) {
    SETERRQ(62, "not ready to solve: invalid w[] in enthSystemCtx"); }
  if (gsl_isnan(Sigma[0])) {
    SETERRQ(63, "not ready to solve: invalid Sigma[] in enthSystemCtx"); }
  if (gsl_isnan(Enth_s[0])) {
    SETERRQ(64, "not ready to solve: invalid Enth_s[] in enthSystemCtx"); }
  if (gsl_isnan(Enth[0])) {
    SETERRQ(65, "not ready to solve: invalid Enth[] in enthSystemCtx"); }
  return 0;
}


//! Set coefficients in discrete equation for \f$E = Y\f$ at base of ice.
/*!
This method should only be called if everything but the basal boundary condition
is already set.
 */
PetscErrorCode enthSystemCtx::setDirichletBasal(PetscScalar Y) {
#ifdef PISM_DEBUG
  PetscErrorCode ierr;
  ierr = checkReadyToSolve(); CHKERRQ(ierr);
  if ((!gsl_isnan(a0)) || (!gsl_isnan(a1)) || (!gsl_isnan(b))) {
    SETERRQ(1, "setting basal boundary conditions twice in enthSystemCtx");
  }
#endif
  a0 = 1.0;
  a1 = 0.0;
  b  = Y;
  return 0;
}


//! Set coefficients in discrete equation for \f$\partial E / \partial z = Y\f$ at base of ice.
/*!
This method combines the Neumann boundary condition with the differential equation.
The vertical advection part is zeroed-out.  The error in the pure conductive and
smooth conductivity case is \f$O(\Delta z^2)\f$.

This code is near-duplication of code in solveThisColumn() below.

This method should only be called if everything but the basal boundary condition
is already set.
 */
PetscErrorCode enthSystemCtx::setNeumannBasal(PetscScalar Y) {
 PetscErrorCode ierr;
#ifdef PISM_DEBUG

  ierr = checkReadyToSolve(); CHKERRQ(ierr);
  if ((!gsl_isnan(a0)) || (!gsl_isnan(a1)) || (!gsl_isnan(b))) {
    SETERRQ(1, "setting basal boundary conditions twice in enthSystemCtx");
  }
#endif
  const PetscScalar
      Rc = (Enth[0] < Enth_s[0]) ? iceRcold : iceRtemp,
      Rr = (Enth[1] < Enth_s[1]) ? iceRcold : iceRtemp,
      Rminus = Rc,
      Rplus  = 0.5 * (Rc + Rr);
  a0 = 1.0 + Rminus + Rplus;  // = D[0]
  a1 = - Rminus - Rplus;      // = U[0]
  const PetscScalar X = - 2.0 * dzEQ * Y;  // E(-dz) = E(+dz) + X
  // zero vertical velocity contribution
  b = Enth[0] + Rminus * X;   // = rhs[0]
  if (!ismarginal) {
    planeStar<PetscScalar> ss;
    ierr = Enth3->getPlaneStar_fine(i,j,0,&ss); CHKERRQ(ierr);
    const PetscScalar UpEnthu = (u[0] < 0) ? u[0] * (ss.e -  ss.ij) / dx :
                                             u[0] * (ss.ij  - ss.w) / dx;
    const PetscScalar UpEnthv = (v[0] < 0) ? v[0] * (ss.n -  ss.ij) / dy :
                                             v[0] * (ss.ij  - ss.s) / dy;
    b += dtTemp * ((Sigma[0] / ice_rho) - UpEnthu - UpEnthv);  // = rhs[0]
  }
  return 0;
}


/*! \brief Solve the tridiagonal system, in a single column, which determines
the new values of the ice enthalpy. */
PetscErrorCode enthSystemCtx::solveThisColumn(PetscScalar **x, PetscErrorCode &pivoterrorindex) {
  PetscErrorCode ierr;
#ifdef PISM_DEBUG
  ierr = checkReadyToSolve(); CHKERRQ(ierr);
  if ((gsl_isnan(a0)) || (gsl_isnan(a1)) || (gsl_isnan(b))) {
    SETERRQ(1, "solveThisColumn() should only be called after\n"
               "  setting basal boundary condition in enthSystemCtx"); }
#endif
  // k=0 equation is already established
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
      planeStar<PetscScalar> ss;
      ierr = Enth3->getPlaneStar_fine(i,j,k,&ss); CHKERRQ(ierr);
      const PetscScalar UpEnthu = (u[k] < 0) ? u[k] * (ss.e -  ss.ij) / dx :
                                               u[k] * (ss.ij  - ss.w) / dx;
      const PetscScalar UpEnthv = (v[k] < 0) ? v[k] * (ss.n -  ss.ij) / dy :
                                               v[k] * (ss.ij  - ss.s) / dy;
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


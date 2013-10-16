// Copyright (C) 2009-2013 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "enthSystem.hh"
#include <gsl/gsl_math.h>
#include "NCVariable.hh"
#include "iceModelVec.hh"
#include "enthalpyConverter.hh"

enthSystemCtx::enthSystemCtx(const NCConfigVariable &config,
                             IceModelVec3 &my_Enth3,
                             PetscScalar my_dx,  PetscScalar my_dy,
                             PetscScalar my_dt,  PetscScalar my_dz,
                             int my_Mz, std::string my_prefix,
                             EnthalpyConverter *my_EC)
  : columnSystemCtx(my_Mz, my_prefix), EC(my_EC) {  // <- critical: sets size of sys
  Mz = my_Mz;

  // set some values so we can check if init was called
  nu       = -1.0;
  R_cold   = -1.0;
  R_temp   = -1.0;
  m_lambda = -1.0;
  a0 = GSL_NAN;
  a1 = GSL_NAN;
  b  = GSL_NAN;

  ice_rho = config.get("ice_density");
  ice_c   = config.get("ice_specific_heat_capacity");
  ice_k   = config.get("ice_thermal_conductivity");
  p_air   = config.get("surface_pressure");

  ice_K   = ice_k / ice_c;
  ice_K0  = ice_K * config.get("enthalpy_temperate_conductivity_ratio");

  u = new PetscScalar[Mz];
  v = new PetscScalar[Mz];
  w = new PetscScalar[Mz];
  Enth   = new PetscScalar[Mz];
  Enth_s = new PetscScalar[Mz]; // enthalpy of pressure-melting-point
  strain_heating = new PetscScalar[Mz];

  R.resize(Mz);

  Enth3 = &my_Enth3;  // points to IceModelVec3

  dx = my_dx;
  dy = my_dy;
  dt = my_dt;
  dz = my_dz;
  nu = dt / dz;

  R_factor = dt / (PetscSqr(dz) * ice_rho);
  R_cold = ice_K * R_factor;
  R_temp = ice_K0 * R_factor;

  if (config.get_flag("use_temperature_dependent_thermal_conductivity"))
    k_depends_on_T = true;
  else
    k_depends_on_T = false;

  // check if c(T) is a constant function:
  if (EC->c_from_T(260) != EC->c_from_T(270))
    c_depends_on_T = true;
  else
    c_depends_on_T = false;
}


enthSystemCtx::~enthSystemCtx() {
  delete [] u;
  delete [] v;
  delete [] w;
  delete [] strain_heating;
  delete [] Enth_s;
  delete [] Enth;
}

/*!
  In this implementation \f$k\f$ does not depend on temperature.
 */
PetscScalar enthSystemCtx::k_from_T(PetscScalar T) {

  if (k_depends_on_T)
    return 9.828 * exp(-0.0057 * T);

  return ice_k;
}

PetscErrorCode enthSystemCtx::initThisColumn(int my_i, int my_j, bool my_ismarginal,
                                             PetscReal my_ice_thickness,
                                             PetscReal till_water_thickness,
                                             IceModelVec3 *u3,
                                             IceModelVec3 *v3,
                                             IceModelVec3 *w3,
                                             IceModelVec3 *strain_heating3) {
  PetscErrorCode ierr;

  ice_thickness = my_ice_thickness;
  ismarginal    = my_ismarginal;

  m_ks = static_cast<PetscInt>(floor(ice_thickness/dz));
  ierr = setIndicesAndClearThisColumn(my_i, my_j, m_ks); CHKERRQ(ierr);

#if (PISM_DEBUG==1)
  // check if ks is valid
  if ((m_ks < 0) || (m_ks >= Mz)) {
    PetscPrintf(PETSC_COMM_SELF,
                "ERROR: ks = %d computed at i = %d, j = %d is invalid,"
                " possibly because of invalid ice thickness (%f meters) or dz (%f meters).\n",
                m_ks, i, j, ice_thickness, dz);
    SETERRQ(PETSC_COMM_SELF, 1, "invalid ks");
  }
#endif

  if (m_ks == 0)
    return 0;

  ierr = u3->getValColumn(i, j, m_ks, u); CHKERRQ(ierr);
  ierr = v3->getValColumn(i, j, m_ks, v); CHKERRQ(ierr);
  ierr = w3->getValColumn(i, j, m_ks, w); CHKERRQ(ierr);
  ierr = strain_heating3->getValColumn(i, j, m_ks, strain_heating); CHKERRQ(ierr);

  ierr = Enth3->getValColumn(i, j, m_ks, Enth); CHKERRQ(ierr);
  ierr = compute_enthalpy_CTS(); CHKERRQ(ierr);
  // if there is subglacial water, don't allow ice base enthalpy to be below
  // pressure-melting; that is, assume subglacial water is at the pressure-
  // melting temperature and enforce continuity of temperature
  if (till_water_thickness > 0.0 && Enth[0] < Enth_s[0]) {
    Enth[0] = Enth_s[0];
  }

  m_lambda = compute_lambda();

  ierr = assemble_R(); CHKERRQ(ierr);
  return 0;
}

//! Compute the CTS value of enthalpy in an ice column.
/*!
Return argument Enth_s[Mz] has the enthalpy value for the pressure-melting
temperature at the corresponding z level.
 */
PetscErrorCode enthSystemCtx::compute_enthalpy_CTS() {

  for (PetscInt k = 0; k <= m_ks; k++) {
    const PetscScalar
      depth = ice_thickness - k * dz,
      p = EC->getPressureFromDepth(depth); // FIXME issue #15
    Enth_s[k] = EC->getEnthalpyCTS(p);
  }

  const PetscScalar Es_air = EC->getEnthalpyCTS(p_air);
  for (PetscInt k = m_ks+1; k < Mz; k++) {
    Enth_s[k] = Es_air;
  }
  return 0;
}

//! Compute the lambda for BOMBPROOF.
/*!
See page \ref bombproofenth.
 */
PetscReal enthSystemCtx::compute_lambda() {
  PetscReal result = 1.0; // start with centered implicit for more accuracy
  const double epsilon = 1e-6 / 3.15569259747e7;

  for (PetscInt k = 0; k <= m_ks; k++) {
    if (Enth[k] > Enth_s[k]) { // lambda = 0 if temperate ice present in column
      result = 0.0;
    } else {
      const PetscScalar denom = (PetscAbs(w[k]) + epsilon) * ice_rho * ice_c * dz;
      result = PetscMin(result, 2.0 * ice_k / denom);
    }
  }
  return result;
}


PetscErrorCode enthSystemCtx::setDirichletSurface(PetscScalar my_Enth_surface) {
#if (PISM_DEBUG==1)
  if ((nu < 0.0) || (R_cold < 0.0) || (R_temp < 0.0)) {
    SETERRQ(PETSC_COMM_SELF, 2, "setDirichletSurface() should only be called after\n"
            "  initAllColumns() in enthSystemCtx");
  }
#endif
  Enth_ks = my_Enth_surface;
  return 0;
}


PetscErrorCode enthSystemCtx::viewConstants(PetscViewer viewer, bool show_col_dependent) {
  PetscErrorCode ierr;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }

  PetscBool iascii;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(PETSC_COMM_SELF, 1,"Only ASCII viewer for enthSystemCtx::viewConstants()\n"); }
  
  ierr = PetscViewerASCIIPrintf(viewer,
                   "\n<<VIEWING enthSystemCtx with prefix '%s':\n",prefix.c_str()); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "for ALL columns:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  dx,dy,dt,dz = %8.2f,%8.2f,%10.3e,%8.2f\n",
                     dx,dy,dt,dz); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  ice_rho,ice_c,ice_k,ice_K,ice_K0 = %10.3e,%10.3e,%10.3e,%10.3e,%10.3e\n",
                     ice_rho,ice_c,ice_k,ice_K,ice_K0); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  nu = %10.3e\n",
                     nu); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  R_cold,R_temp = %10.3e,%10.3e,\n",
		     R_cold,R_temp); CHKERRQ(ierr);
  if (show_col_dependent) {
    ierr = PetscViewerASCIIPrintf(viewer,
                     "for THIS column:\n"
                     "  i,j,ks = %d,%d,%d\n",
                     i,j,m_ks); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  ismarginal,lambda = %d,%10.3f\n",
                     (int)ismarginal,m_lambda); CHKERRQ(ierr);
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
  if (nu < 0.0 || R_cold < 0.0 || R_temp < 0.0) {
    SETERRQ(PETSC_COMM_SELF, 2,
            "not ready to solve: need initAllColumns() in enthSystemCtx");
  }
  if (m_lambda < 0.0) {
    SETERRQ(PETSC_COMM_SELF, 3,
            "not ready to solve: need setSchemeParamsThisColumn() in enthSystemCtx");
  }
  return 0;
}


//! Set coefficients in discrete equation for \f$E = Y\f$ at base of ice.
/*!
This method should only be called if everything but the basal boundary condition
is already set.
 */
PetscErrorCode enthSystemCtx::setDirichletBasal(PetscScalar Y) {
#if (PISM_DEBUG==1)
  PetscErrorCode ierr;
  ierr = checkReadyToSolve(); CHKERRQ(ierr);
  if (gsl_isnan(a0) == 0 || gsl_isnan(a1) == 0 || gsl_isnan(b) == 0) {
    SETERRQ(PETSC_COMM_SELF, 1, "setting basal boundary conditions twice in enthSystemCtx");
  }
#endif
  a0 = 1.0;
  a1 = 0.0;
  b  = Y;
  return 0;
}


//! Set coefficients in discrete equation for Neumann condition at base of ice.
/*!
This method generates the Neumann boundary condition for the linear system.

The Neumann boundary condition is
   \f[ \frac{\partial E}{\partial z} = - \frac{\phi}{K} \f]
where \f$\phi\f$ is the heat flux.  Here \f$K\f$ is allowed to vary, and takes
its value from the value computed in assemble_R().

The boundary condition is combined with the partial differential equation by the
technique of introducing an imaginary point at \f$z=-\Delta z\f$ and then
eliminating it.

The error in the pure conductive and smooth conductivity case is \f$O(\Delta z^2)\f$.

This method should only be called if everything but the basal boundary condition
is already set.
 */
PetscErrorCode enthSystemCtx::setBasalHeatFlux(PetscScalar hf) {
 PetscErrorCode ierr;
#if (PISM_DEBUG==1)
  ierr = checkReadyToSolve(); CHKERRQ(ierr);
  if (gsl_isnan(a0) == 0 || gsl_isnan(a1) == 0 || gsl_isnan(b) == 0) {
    SETERRQ(PETSC_COMM_SELF, 1, "setting basal boundary conditions twice in enthSystemCtx");
  }
#endif
  // extract K from R[0], so this code works even if K=K(T)
  // recall:   R = (ice_K / ice_rho) * dt / PetscSqr(dz)
  const PetscScalar
    K = (ice_rho * PetscSqr(dz) * R[0]) / dt,
    Y = - hf / K;
  const PetscScalar
    Rc = R[0],
    Rr = R[1],
    Rminus = Rc,
    Rplus  = 0.5 * (Rc + Rr);
  a0 = 1.0 + Rminus + Rplus;  // = D[0]
  a1 = - Rminus - Rplus;      // = U[0]
  // next line says 
  //   (E(+dz) - E(-dz)) / (2 dz) = Y
  // or equivalently
  //   E(-dz) = E(+dz) + X
  const PetscScalar X = - 2.0 * dz * Y;
  // zero vertical velocity contribution
  b = Enth[0] + Rminus * X;   // = rhs[0]
  if (!ismarginal) {
    planeStar<PetscScalar> ss;
    ierr = Enth3->getPlaneStar_fine(i,j,0,&ss); CHKERRQ(ierr);
    const PetscScalar UpEnthu = (u[0] < 0) ? u[0] * (ss.e -  ss.ij) / dx :
                                             u[0] * (ss.ij  - ss.w) / dx;
    const PetscScalar UpEnthv = (v[0] < 0) ? v[0] * (ss.n -  ss.ij) / dy :
                                             v[0] * (ss.ij  - ss.s) / dy;
    b += dt * ((strain_heating[0] / ice_rho) - UpEnthu - UpEnthv);  // = rhs[0]
  }
  return 0;
}


//! \brief Assemble the R array.  The R value switches at the CTS.
/*!  In a simple abstract diffusion
  \f[ \frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial z^2}, \f]
with time steps \f$\Delta t\f$ and spatial steps \f$\Delta z\f$ we define
  \f[ R = \frac{D \Delta t}{\Delta z^2}. \f]
This is used in an implicit method to write each line in the linear system, for
example [\ref MortonMayers]:
  \f[ -R U_{j-1}^{n+1} + (1+2R) U_j^{n+1} - R U_{j+1}^{n+1} = U_j^n. \f]
  
In the case of conservation of energy [\ref AschwandenBuelerKhroulevBlatter],
  \f[ u=E \qquad \text{ and } \qquad D = \frac{K}{\rho} \qquad \text{ and } \qquad K = \frac{k}{c}. \f]
Thus
  \f[ R = \frac{k \Delta t}{\rho c \Delta z^2}. \f]
 */
PetscErrorCode enthSystemCtx::assemble_R() {
  PetscErrorCode ierr;

  if (k_depends_on_T == false && c_depends_on_T == false) {

    for (PetscInt k = 0; k <= m_ks; k++)
      R[k] = (Enth[k] < Enth_s[k]) ? R_cold : R_temp;

  } else {

    for (PetscInt k = 0; k <= m_ks; k++) {
      if (Enth[k] < Enth_s[k]) {
        // cold case
        const PetscScalar depth = ice_thickness - k * dz;
        PetscScalar T;
        ierr = EC->getAbsTemp(Enth[k], EC->getPressureFromDepth(depth), // FIXME: issue #15
                              T); CHKERRQ(ierr);

        R[k] = ((k_depends_on_T ? k_from_T(T) : ice_k) / EC->c_from_T(T)) * R_factor;
      } else {
        // temperate case
        R[k] = R_temp;
      }
    }

  }

  // R[k] for k > ks are never used
#if (PISM_DEBUG==1)
  for (int k = m_ks + 1; k < Mz; ++k)
    R[k] = GSL_NAN;
#endif

  return 0;
}


/*! \brief Solve the tridiagonal system, in a single column, which
 *  determines the new values of the ice enthalpy.
 */
PetscErrorCode enthSystemCtx::solveThisColumn(PetscScalar **x) {
  PetscErrorCode ierr;
#if (PISM_DEBUG==1)
  ierr = checkReadyToSolve(); CHKERRQ(ierr);
  if (gsl_isnan(a0) || gsl_isnan(a1) || gsl_isnan(b)) {
    SETERRQ(PETSC_COMM_SELF, 1,
            "solveThisColumn() should only be called after\n"
            "  setting basal boundary condition in enthSystemCtx"); }
#endif

  // k=0 equation is already established
  // L[0] = 0.0;  // not allocated
  D[0]   = a0;
  U[0]   = a1;
  rhs[0] = b;

  // generic ice segment in k location (if any; only runs if ks >= 2)
  for (PetscInt k = 1; k < m_ks; k++) {
    const PetscScalar
        Rminus = 0.5 * (R[k-1] + R[k]  ),
        Rplus  = 0.5 * (R[k]   + R[k+1]);
    L[k] = - Rminus;
    D[k] = 1.0 + Rminus + Rplus;
    U[k] = - Rplus;
    const PetscScalar AA = nu * w[k];
    if (w[k] >= 0.0) {  // velocity upward
      L[k] -= AA * (1.0 - m_lambda/2.0);
      D[k] += AA * (1.0 - m_lambda);
      U[k] += AA * (m_lambda/2.0);
    } else {            // velocity downward
      L[k] -= AA * (m_lambda/2.0);
      D[k] -= AA * (1.0 - m_lambda);
      U[k] += AA * (1.0 - m_lambda/2.0);
    }
    rhs[k] = Enth[k];
    if (!ismarginal) {
      planeStar<PetscScalar> ss;
      ierr = Enth3->getPlaneStar_fine(i,j,k,&ss); CHKERRQ(ierr);
      const PetscScalar UpEnthu = (u[k] < 0) ? u[k] * (ss.e -  ss.ij) / dx :
                                               u[k] * (ss.ij  - ss.w) / dx;
      const PetscScalar UpEnthv = (v[k] < 0) ? v[k] * (ss.n -  ss.ij) / dy :
                                               v[k] * (ss.ij  - ss.s) / dy;
      rhs[k] += dt * ((strain_heating[k] / ice_rho) - UpEnthu - UpEnthv);
    }
  }

  // set Dirichlet boundary condition at top
  if (m_ks > 0) L[m_ks] = 0.0;
  D[m_ks] = 1.0;
  if (m_ks < Mz-1) U[m_ks] = 0.0;
  rhs[m_ks] = Enth_ks;

  // solve it; note drainage is not addressed yet and post-processing may occur
  int pivoterr = solveTridiagonalSystem(m_ks+1, x);

  if (pivoterr != 0) {
    ierr = PetscPrintf(PETSC_COMM_SELF,
                       "\n\ntridiagonal solve of enthSystemCtx in enthalpyAndDrainageStep() FAILED at (%d,%d)\n"
                       " with zero pivot position %d; viewing system to m-file ... \n",
                       i, j, pivoterr); CHKERRQ(ierr);
    ierr = reportColumnZeroPivotErrorMFile(pivoterr); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_SELF, 1,"PISM ERROR in enthalpyDrainageStep()\n");
  }

  // air above
  for (PetscInt k = m_ks+1; k < Mz; k++) {
    (*x)[k] = Enth_ks;
  }

#if (PISM_DEBUG==1)
  if (pivoterr == 0) {
    // if success, mark column as done by making scheme params and b.c. coeffs invalid
    m_lambda  = -1.0;
    a0 = GSL_NAN;
    a1 = GSL_NAN;
    b  = GSL_NAN;
  }
#endif
  return 0;
}

//! View the tridiagonal system A x = b to a PETSc viewer, both A as a full matrix and b as a vector.
PetscErrorCode enthSystemCtx::viewSystem(PetscViewer viewer) const {
  PetscErrorCode ierr;
  std::string info;
  info = prefix + "_A";
  ierr = viewMatrix(viewer,info.c_str()); CHKERRQ(ierr);
  info = prefix + "_rhs";
  ierr = viewVectorValues(viewer,rhs,nmax,info.c_str()); CHKERRQ(ierr);
  info = prefix + "_R";
  ierr = viewVectorValues(viewer,&R[0],Mz,info.c_str()); CHKERRQ(ierr);
  return 0;
}

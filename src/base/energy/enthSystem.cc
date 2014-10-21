// Copyright (C) 2009-2014 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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
#include "PISMConfig.hh"
#include "iceModelVec.hh"
#include "enthalpyConverter.hh"

#include "error_handling.hh"

namespace pism {

enthSystemCtx::enthSystemCtx(const Config &config,
                             IceModelVec3 &my_Enth3,
                             double my_dx,  double my_dy,
                             double my_dt,  double my_dz,
                             int my_Mz, const std::string &my_prefix,
                             const EnthalpyConverter &my_EC)
  : columnSystemCtx(my_Mz, my_prefix), EC(my_EC) {  // <- critical: sets size of sys
  Mz = my_Mz;

  // set some values so we can check if init was called
  R_cold   = -1.0;
  R_temp   = -1.0;
  m_lambda = -1.0;
  D0 = GSL_NAN;
  U0 = GSL_NAN;
  B0 = GSL_NAN;

  ice_rho = config.get("ice_density");
  ice_c   = config.get("ice_specific_heat_capacity");
  ice_k   = config.get("ice_thermal_conductivity");
  p_air   = config.get("surface_pressure");

  ice_K  = ice_k / ice_c;
  ice_K0 = ice_K * config.get("enthalpy_temperate_conductivity_ratio");

  u.resize(Mz);
  v.resize(Mz);
  w.resize(Mz);
  Enth.resize(Mz);
  Enth_s.resize(Mz);
  strain_heating.resize(Mz);
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

  if (config.get_flag("use_temperature_dependent_thermal_conductivity")) {
    k_depends_on_T = true;
  } else {
    k_depends_on_T = false;
  }

  // check if c(T) is a constant function:
  if (EC.c_from_T(260) != EC.c_from_T(270)) {
    c_depends_on_T = true;
  } else {
    c_depends_on_T = false;
  }
}


enthSystemCtx::~enthSystemCtx() {
}

/*!
  In this implementation \f$k\f$ does not depend on temperature.
 */
double enthSystemCtx::k_from_T(double T) {

  if (k_depends_on_T) {
    return 9.828 * exp(-0.0057 * T);
  }

  return ice_k;
}

PetscErrorCode enthSystemCtx::initThisColumn(int my_i, int my_j, bool my_ismarginal,
                                             double my_ice_thickness,
                                             IceModelVec3 *u3,
                                             IceModelVec3 *v3,
                                             IceModelVec3 *w3,
                                             IceModelVec3 *strain_heating3) {
  PetscErrorCode ierr;

  ice_thickness = my_ice_thickness;
  ismarginal    = my_ismarginal;

  setIndicesAndClearThisColumn(my_i, my_j, ice_thickness, dz, Mz);

  if (m_ks == 0) {
    return 0;
  }

  ierr = u3->getValColumn(m_i, m_j, m_ks, &u[0]); CHKERRQ(ierr);
  ierr = v3->getValColumn(m_i, m_j, m_ks, &v[0]); CHKERRQ(ierr);
  ierr = w3->getValColumn(m_i, m_j, m_ks, &w[0]); CHKERRQ(ierr);
  ierr = strain_heating3->getValColumn(m_i, m_j, m_ks, &strain_heating[0]); CHKERRQ(ierr);

  ierr = Enth3->getValColumn(m_i, m_j, m_ks, &Enth[0]); CHKERRQ(ierr);
  ierr = compute_enthalpy_CTS(); CHKERRQ(ierr);

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

  for (unsigned int k = 0; k <= m_ks; k++) {
    const double
      depth = ice_thickness - k * dz,
      p = EC.getPressureFromDepth(depth); // FIXME issue #15
    Enth_s[k] = EC.getEnthalpyCTS(p);
  }

  const double Es_air = EC.getEnthalpyCTS(p_air);
  for (unsigned int k = m_ks+1; k < Enth_s.size(); k++) {
    Enth_s[k] = Es_air;
  }
  return 0;
}

//! Compute the lambda for BOMBPROOF.
/*!
See page \ref bombproofenth.
 */
double enthSystemCtx::compute_lambda() {
  double result = 1.0; // start with centered implicit for more accuracy
  const double epsilon = 1e-6 / 3.15569259747e7;

  for (unsigned int k = 0; k <= m_ks; k++) {
    if (Enth[k] > Enth_s[k]) { // lambda = 0 if temperate ice present in column
      result = 0.0;
    } else {
      const double denom = (PetscAbs(w[k]) + epsilon) * ice_rho * ice_c * dz;
      result = PetscMin(result, 2.0 * ice_k / denom);
    }
  }
  return result;
}


PetscErrorCode enthSystemCtx::setDirichletSurface(double my_Enth_surface) {
#if (PISM_DEBUG==1)
  if ((nu < 0.0) || (R_cold < 0.0) || (R_temp < 0.0)) {
    throw RuntimeError("setDirichletSurface() should only be called after\n"
                       "initAllColumns() in enthSystemCtx");
  }
#endif
  Enth_ks = my_Enth_surface;
  return 0;
}


PetscErrorCode enthSystemCtx::viewConstants(PetscViewer viewer, bool show_col_dependent) {
  PetscErrorCode ierr;

  if (not viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }

  PetscBool iascii;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii); CHKERRQ(ierr);
  if (not iascii) {
    throw RuntimeError("Only ASCII viewer for enthSystemCtx::viewConstants()");
  }
  
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
                     m_i,m_j,m_ks); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  ismarginal,lambda = %d,%10.3f\n",
                     (int)ismarginal,m_lambda); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  Enth_ks = %10.3e\n",
                     Enth_ks); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
                     "  D0,U0,B0 = %10.3e,%10.3e\n",
                     D0,U0,B0); CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,
                     ">>\n\n"); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode enthSystemCtx::checkReadyToSolve() {
  if (nu < 0.0 || R_cold < 0.0 || R_temp < 0.0) {
    throw RuntimeError("not ready to solve: need initAllColumns() in enthSystemCtx");
  }
  if (m_lambda < 0.0) {
    throw RuntimeError("not ready to solve: need setSchemeParamsThisColumn() in enthSystemCtx");
  }
  return 0;
}


//! Set coefficients in discrete equation for \f$E = Y\f$ at base of ice.
/*!
This method should only be called if everything but the basal boundary condition
is already set.
 */
PetscErrorCode enthSystemCtx::setDirichletBasal(double Y) {
#if (PISM_DEBUG==1)
  PetscErrorCode ierr;
  ierr = checkReadyToSolve(); CHKERRQ(ierr);
  if (gsl_isnan(D0) == 0 || gsl_isnan(U0) == 0 || gsl_isnan(B0) == 0) {
    throw RuntimeError("setting basal boundary conditions twice in enthSystemCtx");
  }
#endif
  D0 = 1.0;
  U0 = 0.0;
  B0  = Y;
  return 0;
}


//! Set coefficients in discrete equation for Neumann condition at base of ice.
/*!
This method generates the Neumann boundary condition for the linear system.

The Neumann boundary condition is
   @f[ \frac{\partial E}{\partial z} = - \frac{\phi}{K} @f]
where \f$\phi\f$ is the heat flux.  Here \f$K\f$ is allowed to vary, and takes
its value from the value computed in assemble_R().

The boundary condition is combined with the partial differential equation by the
technique of introducing an imaginary point at \f$z=-\Delta z\f$ and then
eliminating it.

In other words, we combine the centered finite difference approximation
@f[ \frac{ E_{1} - E_{-1} }{2\dz}  = -\frac{\phi}{K} @f]
with

@f[ -R_{k-\frac12} E_{k-1} + \left( 1 + R_{k-\frac12} + R_{k+\frac12} \right) E_{k} - R_{k+\frac12} E_{k+1} + \text{advective terms} = \dots @f]

to get

@f{align*}{
   \frac{E_{1}-E_{-1}}{2\,\Delta z} & = -\frac{\phi}{K_{0}}, \\
   E_{1}-E_{-1} & = -\frac{2\,\Delta z\,\phi}{K_{0}}, \\
    E_{-1}\,R_{-\frac12}-R_{-\frac12}\,E_{1} & = \frac{2\,R_{-\frac12}\,\Delta z\,\phi}{K_{0}}, \\
    -R_{\frac12}\,E_{1}+E_{0}\,\left(R_{\frac12}+R_{-\frac12}+1\right)-E_{-1}\,R_{-\frac12} + \text{advective terms} & = \dots, \\
    \left(-R_{\frac12}-R_{-\frac12}\right)\,E_{1}+E_{0}\,\left(R_{\frac12}+R_{-\frac12}+1\right) + \text{advective terms} & = \frac{2\,R_{-\frac12}\,\Delta z\,\phi}{K_{0}}+\dots.
@f}

The error in the pure conductive and smooth conductivity case is @f$ O(\dz^2) @f$.

This method should only be called if everything but the basal boundary condition
is already set.

 */
PetscErrorCode enthSystemCtx::setBasalHeatFlux(double heat_flux) {
 PetscErrorCode ierr;
#if (PISM_DEBUG==1)
  ierr = checkReadyToSolve(); CHKERRQ(ierr);
  if (gsl_isnan(D0) == 0 || gsl_isnan(U0) == 0 || gsl_isnan(B0) == 0) {
    throw RuntimeError("setting basal boundary conditions twice in enthSystemCtx");
  }
#endif
  // extract K from R[0], so this code works even if K=K(T)
  // recall:   R = (ice_K / ice_rho) * dt / PetscSqr(dz)
  const double
    K      = (ice_rho * PetscSqr(dz) * R[0]) / dt,
    Rc     = R[0],
    Rr     = R[1],
    Rminus = Rc,
    Rplus  = 0.5 * (Rc + Rr);
  D0 = 1.0 + Rminus + Rplus;
  // modified upper-diagonal term:
  U0 = - Rminus - Rplus;
  // Enth[0] (below) is there due to the fully-implicit discretization
  // in time, the second term is the modification of the right-hand
  // side implementing the Neumann B.C. (see the doxygen comment)
  B0 = Enth[0] + 2.0 * Rminus * dz * heat_flux / K;
  // treat vertical velocity using first-order upwinding:
  if (not ismarginal) {
    planeStar<double> ss;
    ierr = Enth3->getPlaneStar_fine(m_i,m_j,0,&ss); CHKERRQ(ierr);
    const double UpEnthu = (u[0] < 0) ? u[0] * (ss.e -  ss.ij) / dx :
                                             u[0] * (ss.ij  - ss.w) / dx;
    const double UpEnthv = (v[0] < 0) ? v[0] * (ss.n -  ss.ij) / dy :
                                             v[0] * (ss.ij  - ss.s) / dy;
    B0 += dt * ((strain_heating[0] / ice_rho) - UpEnthu - UpEnthv);  // = rhs[0]
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

    for (unsigned int k = 0; k <= m_ks; k++)
      R[k] = (Enth[k] < Enth_s[k]) ? R_cold : R_temp;

  } else {

    for (unsigned int k = 0; k <= m_ks; k++) {
      if (Enth[k] < Enth_s[k]) {
        // cold case
        const double depth = ice_thickness - k * dz;
        double T;
        ierr = EC.getAbsTemp(Enth[k], EC.getPressureFromDepth(depth), // FIXME: issue #15
                              T); CHKERRQ(ierr);

        R[k] = ((k_depends_on_T ? k_from_T(T) : ice_k) / EC.c_from_T(T)) * R_factor;
      } else {
        // temperate case
        R[k] = R_temp;
      }
    }

  }

  // R[k] for k > m_ks are never used
#if (PISM_DEBUG==1)
  for (unsigned int k = m_ks + 1; k < R.size(); ++k)
    R[k] = GSL_NAN;
#endif

  return 0;
}


/*! \brief Solve the tridiagonal system, in a single column, which
 *  determines the new values of the ice enthalpy.
 *
 * To discretize
 * @f[ \diff{}{z} \left( K(E) \diff{E}{z}\right) = \diff{E}{t} @f]
 *
 * at a location @f$ k @f$ of the vertical grid, we use centered
 * finite differences and evaluate @f$ K(E) @f$ at
 * staggered-grid locations:
 *
 * @f[ \frac{K_{k+\frac12}\frac{E_{k+1} - E_{k}}{\dz} - K_{k-\frac12}\frac{E_{k} - E_{k-1}}{\dz}}{\dz}, @f]
 *
 * where @f$ K_{k\pm \frac12} = K(E_{k\pm \frac12}) @f$.
 *
 * We define
 *
 * @f[ R_i = \frac{\dt\, K_i}{\dz^2}, @f]
 *
 * and the discretization takes form
 *
 * @f[ -R_{k-\frac12} E_{k-1} + \left( 1 + R_{k-\frac12} + R_{k+\frac12} \right) E_{k} - R_{k+\frac12} E_{k+1} = @f].
 *
 * In the assembly of the tridiagonal system, this corresponds to
 *
 * @f{align*}{
 * L_i &= - \frac12 (R_{i} + R_{i-1}),\\
 * D_i &= 1 + \frac12 (R_{i} + R_{i-1}) + \frac12 (R_{i} + R_{i+1}),\\
 * U_i &= - \frac12 (R_{i} + R_{i+1}),
 * @f}
 *
 * where @f$ L_i, D_i, U_i @f$ are lower-diagonal, diagonal, and
 * upper-diagonal entries corresponding to an equation @f$ i @f$.
 * (Staggered-grid values are approximated by interpolating from the
 * regular grid).
 */
PetscErrorCode enthSystemCtx::solveThisColumn(std::vector<double> &x) {
  PetscErrorCode ierr;
#if (PISM_DEBUG==1)
  ierr = checkReadyToSolve(); CHKERRQ(ierr);
  if (gsl_isnan(D0) || gsl_isnan(U0) || gsl_isnan(B0)) {
    throw RuntimeError("solveThisColumn() should only be called after\n"
                       "  setting basal boundary condition in enthSystemCtx");
  }
#endif

  // k=0 equation is already established
  // L[0] = 0.0;  // not allocated
  D[0]   = D0;
  U[0]   = U0;
  rhs[0] = B0;

  // generic ice segment in k location (if any; only runs if m_ks >= 2)
  for (unsigned int k = 1; k < m_ks; k++) {
    const double
        Rminus = 0.5 * (R[k-1] + R[k]),
        Rplus  = 0.5 * (R[k]   + R[k+1]);
    L[k] = - Rminus;
    D[k] = 1.0 + Rminus + Rplus;
    U[k] = - Rplus;
    const double AA = nu * w[k];
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
    if (not ismarginal) {
      planeStar<double> ss;
      ierr = Enth3->getPlaneStar_fine(m_i,m_j,k,&ss); CHKERRQ(ierr);
      const double UpEnthu = (u[k] < 0) ? u[k] * (ss.e -  ss.ij) / dx :
                                               u[k] * (ss.ij  - ss.w) / dx;
      const double UpEnthv = (v[k] < 0) ? v[k] * (ss.n -  ss.ij) / dy :
                                               v[k] * (ss.ij  - ss.s) / dy;
      rhs[k] += dt * ((strain_heating[k] / ice_rho) - UpEnthu - UpEnthv);
    }
  }

  // set Dirichlet boundary condition at top
  if (m_ks > 0) {
    L[m_ks] = 0.0;
  }
  D[m_ks] = 1.0;
  if (m_ks < Mz-1) {
    U[m_ks] = 0.0;
  }
  rhs[m_ks] = Enth_ks;

  // Solve it; note drainage is not addressed yet and post-processing may occur
  try {
    solveTridiagonalSystem(m_ks+1, x);
  }
  catch (RuntimeError &e) {
    e.add_context("solving the tri-diagonal system (enthSystemCtx) at (%d,%d)\n"
                  "viewing system to m-file... ", m_i, m_j);
    reportColumnZeroPivotErrorMFile(m_ks + 1);
    throw;
  }

  // air above
  for (unsigned int k = m_ks+1; k < x.size(); k++) {
    x[k] = Enth_ks;
  }

#if (PISM_DEBUG==1)
  // if success, mark column as done by making scheme params and b.c. coeffs invalid
  m_lambda = -1.0;
  D0       = GSL_NAN;
  U0       = GSL_NAN;
  B0       = GSL_NAN;
#endif
  return 0;
}

//! View the tridiagonal system A x = b to a PETSc viewer, both A as a full matrix and b as a vector.
PetscErrorCode enthSystemCtx::viewSystem(PetscViewer viewer,
                                         unsigned int M) const {
  PetscErrorCode ierr;
  std::string info;
  info = prefix + "_A";
  ierr = viewMatrix(viewer, M, info.c_str()); CHKERRQ(ierr);
  info = prefix + "_rhs";
  viewVectorValues(viewer, rhs, M, info.c_str());
  info = prefix + "_R";
  viewVectorValues(viewer, R, M, info.c_str());
  return 0;
}

} // end of namespace pism

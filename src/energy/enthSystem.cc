// Copyright (C) 2009-2018 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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
#include <gsl/gsl_math.h>       // GSL_NAN, gsl_isnan()
#include "pism/util/ConfigInterface.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/EnthalpyConverter.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/ColumnInterpolation.hh"

namespace pism {
namespace energy {

enthSystemCtx::enthSystemCtx(const std::vector<double>& storage_grid,
                             const std::string &prefix,
                             double dx,  double dy, double dt,
                             const Config &config,
                             const IceModelVec3 &Enth3,
                             const IceModelVec3 &u3,
                             const IceModelVec3 &v3,
                             const IceModelVec3 &w3,
                             const IceModelVec3 &strain_heating3,
                             EnthalpyConverter::Ptr EC)
: columnSystemCtx(storage_grid, prefix, dx, dy, dt, u3, v3, w3),
  m_Enth3(Enth3),
  m_strain_heating3(strain_heating3),
  m_EC(EC) {

  // set some values so we can check if init was called
  m_R_cold   = -1.0;
  m_R_temp   = -1.0;
  m_lambda = -1.0;
  m_D0 = GSL_NAN;
  m_U0 = GSL_NAN;
  m_B0 = GSL_NAN;
  m_L_ks = GSL_NAN;
  m_D_ks = GSL_NAN;
  m_U_ks = GSL_NAN;
  m_B_ks = GSL_NAN;

  m_ice_density = config.get_double("constants.ice.density");
  m_ice_c   = config.get_double("constants.ice.specific_heat_capacity");
  m_ice_k   = config.get_double("constants.ice.thermal_conductivity");
  m_p_air   = config.get_double("surface.pressure");

  m_exclude_horizontal_advection = config.get_boolean("energy.margin_exclude_horizontal_advection");
  m_exclude_vertical_advection   = config.get_boolean("energy.margin_exclude_vertical_advection");
  m_exclude_strain_heat          = config.get_boolean("energy.margin_exclude_strain_heating");

  size_t Mz = m_z.size();
  m_Enth.resize(Mz);
  m_Enth_s.resize(Mz);
  m_strain_heating.resize(Mz);
  m_R.resize(Mz);

  m_E_n.resize(Mz);
  m_E_e.resize(Mz);
  m_E_s.resize(Mz);
  m_E_w.resize(Mz);

  m_nu = m_dt / m_dz;

  double
    ratio = config.get_double(prefix + ".temperate_ice_thermal_conductivity_ratio"),
    K     = m_ice_k / m_ice_c,
    K0    = (ratio * m_ice_k) / m_ice_c;

  m_R_factor = m_dt / (PetscSqr(m_dz) * m_ice_density);
  m_R_cold = K * m_R_factor;
  m_R_temp = K0 * m_R_factor;

  if (config.get_boolean("energy.temperature_dependent_thermal_conductivity")) {
    m_k_depends_on_T = true;
  } else {
    m_k_depends_on_T = false;
  }
}


enthSystemCtx::~enthSystemCtx() {
}

/*!
  In this implementation \f$k\f$ does not depend on temperature.
 */
double enthSystemCtx::k_from_T(double T) const {

  if (m_k_depends_on_T) {
    return 9.828 * exp(-0.0057 * T);
  }

  return m_ice_k;
}

void enthSystemCtx::init(int i, int j, bool marginal, double ice_thickness) {
  m_ice_thickness = ice_thickness;

  m_marginal = marginal;

  init_column(i, j, m_ice_thickness);

  if (m_ks == 0) {
    return;
  }

  coarse_to_fine(m_u3, m_i, m_j, &m_u[0]);
  coarse_to_fine(m_v3, m_i, m_j, &m_v[0]);

  if (m_marginal and m_exclude_vertical_advection) {
    for (unsigned int k = 0; k < m_w.size(); ++k) {
      m_w[k] = 0.0;
    }
  } else {
    coarse_to_fine(m_w3, m_i, m_j, &m_w[0]);
  }

  coarse_to_fine(m_strain_heating3, m_i, m_j, &m_strain_heating[0]);
  coarse_to_fine(m_Enth3, m_i, m_j, &m_Enth[0]);

  coarse_to_fine(m_Enth3, m_i, m_j+1, &m_E_n[0]);
  coarse_to_fine(m_Enth3, m_i+1, m_j, &m_E_e[0]);
  coarse_to_fine(m_Enth3, m_i, m_j-1, &m_E_s[0]);
  coarse_to_fine(m_Enth3, m_i-1, m_j, &m_E_w[0]);

  compute_enthalpy_CTS();

  m_lambda = compute_lambda();

  assemble_R();
}

//! Compute the CTS value of enthalpy in an ice column.
/*! m_Enth_s is set to the enthalpy value for the pressure-melting
  temperature with zero water fraction at the corresponding z level.
 */
void enthSystemCtx::compute_enthalpy_CTS() {

  for (unsigned int k = 0; k <= m_ks; k++) {
    const double
      depth = m_ice_thickness - k * m_dz,
      p = m_EC->pressure(depth); // FIXME issue #15
    m_Enth_s[k] = m_EC->enthalpy_cts(p);
  }

  const double Es_air = m_EC->enthalpy_cts(m_p_air);
  for (unsigned int k = m_ks+1; k < m_Enth_s.size(); k++) {
    m_Enth_s[k] = Es_air;
  }
}

//! Compute the lambda for BOMBPROOF.
/*!
See page \ref bombproofenth.
 */
double enthSystemCtx::compute_lambda() {
  double result = 1.0; // start with centered implicit for more accuracy
  const double epsilon = 1e-6 / 3.15569259747e7;

  for (unsigned int k = 0; k <= m_ks; k++) {
    if (m_Enth[k] > m_Enth_s[k]) { // lambda = 0 if temperate ice present in column
      result = 0.0;
    } else {
      const double denom = (fabs(m_w[k]) + epsilon) * m_ice_density * m_ice_c * m_dz;
      result = std::min(result, 2.0 * m_ice_k / denom);
    }
  }
  return result;
}


void enthSystemCtx::set_surface_dirichlet_bc(double E_surface) {
#if (PISM_DEBUG==1)
  if ((m_nu < 0.0) || (m_R_cold < 0.0) || (m_R_temp < 0.0)) {
    throw RuntimeError(PISM_ERROR_LOCATION, "setDirichletSurface() should only be called after\n"
                       "initAllColumns() in enthSystemCtx");
  }
#endif
  m_L_ks = 0.0;
  m_D_ks = 1.0;
  m_U_ks = 0.0;
  m_B_ks = E_surface;
}

static inline double upwind(double u, double E_m, double E, double E_p, double delta_inverse) {
  return u * delta_inverse * (u < 0 ? (E_p -  E) : (E  - E_m));
}


//! Set the top surface heat flux *into* the ice.
/** @param[in] heat_flux prescribed heat flux (positive means flux into the ice)
 */
void enthSystemCtx::set_surface_heat_flux(double heat_flux) {
  // extract K from R[ks], so this code works even if K=K(T)
  // recall:   R = (ice_K / ice_density) * dt / PetscSqr(dz)
  const double
    K = (m_ice_density * PetscSqr(m_dz) * m_R[m_ks]) / m_dt,
    G = heat_flux / K;

  this->set_surface_neumann_bc(G);
}

//! Set enthalpy flux at the surface.
/*! This method should probably be used for debugging only. Its purpose is to allow setting the
    enthalpy flux even if K == 0, i.e. in a "pure advection" setup.
 */
void enthSystemCtx::set_surface_neumann_bc(double G) {
  const bool include_horizontal_advection = not (m_marginal and m_exclude_horizontal_advection);
  const bool include_strain_heating       = not (m_marginal and m_exclude_strain_heat);

  const double
    Rminus = 0.5 * (m_R[m_ks - 1] + m_R[m_ks]), // R_{ks-1/2}
    Rplus  = m_R[m_ks],                         // R_{ks+1/2}
    mu_w   = 0.5 * m_nu * m_w[m_ks];

  const double A_l = m_w[m_ks] < 0.0 ? 1.0 - m_lambda : m_lambda - 1.0;
  const double A_d = m_w[m_ks] < 0.0 ? m_lambda - 1.0 : 1.0 - m_lambda;
  const double A_b = m_w[m_ks] < 0.0 ? m_lambda - 2.0 : -m_lambda;

  // modified lower-diagonal entry:
  m_L_ks = - Rminus - Rplus + 2.0 * mu_w * A_l;
  // diagonal entry
  m_D_ks = 1.0 + Rminus + Rplus + 2.0 * mu_w * A_d;
  // upper-diagonal entry (not used)
  m_U_ks = 0.0;
  // m_Enth[m_ks] (below) is there due to the fully-implicit discretization in time, the second term is
  // the modification of the right-hand side implementing the Neumann B.C. (similar to
  // set_basal_heat_flux(); see that method for details)
  m_B_ks = m_Enth[m_ks] + 2.0 * G * m_dz * (Rplus + mu_w * A_b);

  // treat horizontal velocity using first-order upwinding:
  double upwind_u = 0.0;
  double upwind_v = 0.0;
  if (include_horizontal_advection) {
    upwind_u = upwind(m_u[m_ks], m_E_w[m_ks], m_Enth[m_ks], m_E_e[m_ks], 1.0 / m_dx);
    upwind_v = upwind(m_v[m_ks], m_E_s[m_ks], m_Enth[m_ks], m_E_n[m_ks], 1.0 / m_dy);
  }
  double Sigma    = 0.0;
  if (include_strain_heating) {
    Sigma = m_strain_heating[m_ks];
  }

  m_B_ks += m_dt * ((Sigma / m_ice_density) - upwind_u - upwind_v);  // = rhs[m_ks]
}

void enthSystemCtx::checkReadyToSolve() {
  if (m_nu < 0.0 || m_R_cold < 0.0 || m_R_temp < 0.0) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "not ready to solve: need initAllColumns() in enthSystemCtx");
  }
  if (m_lambda < 0.0) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "not ready to solve: need setSchemeParamsThisColumn() in enthSystemCtx");
  }
}


//! Set coefficients in discrete equation for \f$E = Y\f$ at base of ice.
/*!
This method should only be called if everything but the basal boundary condition
is already set.
 */
void enthSystemCtx::set_basal_dirichlet_bc(double Y) {
#if (PISM_DEBUG==1)
  checkReadyToSolve();
  if (gsl_isnan(m_D0) == 0 || gsl_isnan(m_U0) == 0 || gsl_isnan(m_B0) == 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "setting basal boundary conditions twice in enthSystemCtx");
  }
#endif
  m_D0 = 1.0;
  m_U0 = 0.0;
  m_B0  = Y;
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

@param[in] heat_flux prescribed heat flux (positive means flux into the ice)

 */
void enthSystemCtx::set_basal_heat_flux(double heat_flux) {
  // extract K from R[0], so this code works even if K=K(T)
  // recall:   R = (ice_K / ice_density) * dt / PetscSqr(dz)
  const double
    K = (m_ice_density * PetscSqr(m_dz) * m_R[0]) / m_dt,
    G = - heat_flux / K;

  this->set_basal_neumann_bc(G);
}

//! Set enthalpy flux at the base.
/*! This method should probably be used for debugging only. Its purpose is to allow setting the
    enthalpy flux even if K == 0, i.e. in a "pure advection" setup.
 */
void enthSystemCtx::set_basal_neumann_bc(double G) {
  const bool include_horizontal_advection = not (m_marginal and m_exclude_horizontal_advection);
  const bool include_strain_heating       = not (m_marginal and m_exclude_strain_heat);

  const double
    Rminus = m_R[0],                  // R_{-1/2}
    Rplus  = 0.5 * (m_R[0] + m_R[1]), // R_{+1/2}
    mu_w   = 0.5 * m_nu * m_w[0];

  const double A_d = m_w[0] < 0.0 ? m_lambda - 1.0 : 1.0 - m_lambda;
  const double A_u = m_w[0] < 0.0 ? 1.0 - m_lambda : m_lambda - 1.0;
  const double A_b = m_w[0] < 0.0 ? -m_lambda : m_lambda - 2.0;

  // diagonal entry
  m_D0 = 1.0 + Rminus + Rplus + 2.0 * mu_w * A_d;
  // upper-diagonal entry
  m_U0 = - Rminus - Rplus + 2.0 * mu_w * A_u;
  // right-hand side, excluding the strain heating term and the horizontal advection
  m_B0 = m_Enth[0] + 2.0 * G * m_dz * (-Rminus + mu_w * A_b);

  // treat horizontal velocity using first-order upwinding:
  double upwind_u = 0.0;
  double upwind_v = 0.0;
  if (include_horizontal_advection) {
    upwind_u = upwind(m_u[0], m_E_w[0], m_Enth[0], m_E_e[0], 1.0 / m_dx);
    upwind_v = upwind(m_v[0], m_E_s[0], m_Enth[0], m_E_n[0], 1.0 / m_dy);
  }
  double Sigma    = 0.0;
  if (include_strain_heating) {
    Sigma = m_strain_heating[0];
  }

  m_B0 += m_dt * ((Sigma / m_ice_density) - upwind_u - upwind_v);  // = rhs[m_ks]
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
void enthSystemCtx::assemble_R() {
  if (not m_k_depends_on_T) {

    for (unsigned int k = 1; k <= m_ks; k++) {
      m_R[k] = (m_Enth[k] < m_Enth_s[k]) ? m_R_cold : m_R_temp;
    }
    //still the cold ice value, if no temperate layer above
    m_R[0] = (m_Enth[1] < m_Enth_s[1]) ? m_R_cold : m_R_temp;

  } else {

    for (unsigned int k = 1; k <= m_ks; k++) {
      if (m_Enth[k] < m_Enth_s[k]) {
        // cold case
        const double depth = m_ice_thickness - k * m_dz;
        double T = m_EC->temperature(m_Enth[k],
                                     m_EC->pressure(depth)); // FIXME: issue #15

        m_R[k] = ((m_k_depends_on_T ? k_from_T(T) : m_ice_k) / m_EC->c()) * m_R_factor;
      } else {
        // temperate case
        m_R[k] = m_R_temp;
      }
    }
    // still the cold ice value, if no temperate layer above
    if (m_Enth[1] < m_Enth_s[1]) {
      double T = m_EC->temperature(m_Enth[0],
                                   m_EC->pressure(m_ice_thickness)); // FIXME: issue #15
      m_R[0] = ((m_k_depends_on_T ? k_from_T(T) : m_ice_k) / m_EC->c()) * m_R_factor;
    } else {
      // temperate layer case
      m_R[0] = m_R_temp;
    }

  }

  // R[k] for k > m_ks are never used
#if (PISM_DEBUG==1)
  for (unsigned int k = m_ks + 1; k < m_R.size(); ++k) {
    m_R[k] = GSL_NAN;
  }
#endif
}


/*! \brief Solve the tridiagonal system, in a single column, which
 *  determines the new values of the ice enthalpy.
 *
 * We are solving a convection-diffusion equation, treating the @f$ z @f$ direction implicitly and
 * @f$ x, y @f$ directions explicitly. See @ref bombproofenth for the documentation of the treatment
 * of convection terms. The notes below document the way we treat diffusion using a simplified PDE.
 *
 * To discretize
 * @f[ \diff{}{z} \left( K(E) \diff{E}{z}\right) = \diff{E}{t} @f]
 *
 * at a location @f$ k @f$ of the vertical grid, we use centered finite differences in space,
 * backward differences in time, and evaluate @f$ K(E) @f$ at staggered-grid locations:
 *
 * @f[ \frac{K_{k+\frac12}\frac{E_{k+1} - E_{k}}{\dz} - K_{k-\frac12}\frac{E_{k} - E_{k-1}}{\dz}}{\dz} = \frac{E_{k} - E^{n-1}_{k}}{\dt}, @f]
 *
 * where @f$ E = E^{n} @f$ for clarity and @f$ K_{k\pm \frac12} = K(E^{n-1}_{k\pm \frac12}) @f$,
 * %i.e. we linearize this PDE by evaluating @f$ K(E) @f$ at the _previous_ time-step.
 *
 * We define
 *
 * @f[ R_i = \frac{\dt\, K_i}{\dz^2}, @f]
 *
 * and the discretization takes form
 *
 * @f[ -R_{k-\frac12} E_{k-1} + \left( 1 + R_{k-\frac12} + R_{k+\frac12} \right) E_{k} - R_{k+\frac12} E_{k+1} = E^{n-1}_{k}. @f]
 *
 * In the assembly of the tridiagonal system this corresponds to
 *
 * @f{align*}{
 * L_i &= - \frac12 (R_{i} + R_{i-1}),\\
 * D_i &= 1 + \frac12 (R_{i} + R_{i-1}) + \frac12 (R_{i} + R_{i+1}),\\
 * U_i &= - \frac12 (R_{i} + R_{i+1}),\\
 * b_i &= E^{n-1}_{i},
 * @f}
 *
 * where @f$ L_i, D_i, U_i @f$ are lower-diagonal, diagonal, and upper-diagonal entries
 * corresponding to an equation @f$ i @f$ and @f$ b_i @f$ is the corresponding right-hand side.
 * (Staggered-grid values are approximated by interpolating from the regular grid).
 *
 * This method is _unconditionally stable_ and has a maximum principle (see [@ref MortonMayers,
 * section 2.11]).
 */
void enthSystemCtx::solve(std::vector<double> &x) {

  TridiagonalSystem &S = *m_solver;

#if (PISM_DEBUG==1)
  checkReadyToSolve();
  if (gsl_isnan(m_D0) || gsl_isnan(m_U0) || gsl_isnan(m_B0)) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "solveThisColumn() should only be called after\n"
                       "  setting basal boundary condition in enthSystemCtx");
  }
#endif

  // k=0 equation is already established
  // L[0] = 0.0;  // not used
  S.D(0)   = m_D0;
  S.U(0)   = m_U0;
  S.RHS(0) = m_B0;

  const double
    one_over_rho = 1.0 / m_ice_density,
    Dx = 1.0 / m_dx,
    Dy = 1.0 / m_dy;

  const bool include_horizontal_advection = not (m_marginal and m_exclude_horizontal_advection);
  const bool include_strain_heating       = not (m_marginal and m_exclude_strain_heat);

  // generic ice segment in k location (if any; only runs if m_ks >= 2)
  for (unsigned int k = 1; k < m_ks; k++) {
    const double
      Rminus = 0.5 * (m_R[k-1] + m_R[k]),   // R_{k-1/2}
      Rplus  = 0.5 * (m_R[k]   + m_R[k+1]), // R_{k+1/2}
      nu_w   = m_nu * m_w[k];

    const double
      A_l = m_w[k] >= 0.0 ? 0.5 * m_lambda - 1.0 : -0.5 * m_lambda,
      A_d = m_w[k] >= 0.0 ? 1.0 - m_lambda : m_lambda - 1.0,
      A_u = m_w[k] >= 0.0 ? 0.5 * m_lambda : 1.0 - 0.5 * m_lambda;

    S.L(k) = - Rminus + nu_w * A_l;
    S.D(k) = 1.0 + Rminus + Rplus + nu_w * A_d;
    S.U(k) = - Rplus + nu_w * A_u;

    // horizontal velocity and strain heating
    double upwind_u = 0.0;
    double upwind_v = 0.0;
    if (include_horizontal_advection) {
      upwind_u = upwind(m_u[k], m_E_w[k], m_Enth[k], m_E_e[k], Dx);
      upwind_v = upwind(m_v[k], m_E_s[k], m_Enth[k], m_E_n[k], Dy);
    }
    double Sigma    = 0.0;
    if (include_strain_heating) {
      Sigma    = m_strain_heating[k];
    }

    S.RHS(k) = m_Enth[k] + m_dt * (one_over_rho * Sigma - upwind_u - upwind_v);
  }

  // Assemble the top surface equation. Values m_{L,D,U,B}_ks are set using set_surface_dirichlet()
  // or set_surface_heat_flux().
  if (m_ks > 0) {
    S.L(m_ks) = m_L_ks;
  }
  S.D(m_ks) = m_D_ks;
  if (m_ks < m_z.size() - 1) {
    S.U(m_ks) = m_U_ks;
  }
  S.RHS(m_ks) = m_B_ks;

  // Solve it; note drainage is not addressed yet and post-processing may occur
  try {
    S.solve(m_ks + 1, x);
  }
  catch (RuntimeError &e) {
    e.add_context("solving the tri-diagonal system (enthSystemCtx) at (%d,%d)\n"
                  "saving system to m-file... ", m_i, m_j);
    reportColumnZeroPivotErrorMFile(m_ks + 1);
    throw;
  }

  // air above
  for (unsigned int k = m_ks+1; k < x.size(); k++) {
    x[k] = m_B_ks;
  }

#if (PISM_DEBUG==1)
  // if success, mark column as done by making scheme params and b.c. coeffs invalid
  m_lambda = -1.0;
  m_D0     = GSL_NAN;
  m_U0     = GSL_NAN;
  m_B0     = GSL_NAN;
  m_L_ks = GSL_NAN;
  m_D_ks = GSL_NAN;
  m_U_ks = GSL_NAN;
  m_B_ks = GSL_NAN;
#endif
}

void enthSystemCtx::save_system(std::ostream &output, unsigned int system_size) const {
  m_solver->save_system(output, system_size);
  m_solver->save_vector(output, m_R, system_size, m_solver->prefix() + "_R");
}

} // end of namespace energy
} // end of namespace pism

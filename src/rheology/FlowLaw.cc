// Copyright (C) 2004-2018 Jed Brown, Ed Bueler, and Constantine Khroulev
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

#include "FlowLaw.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/iceModelVec.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"

#include "pism/util/error_handling.hh"

namespace pism {
namespace rheology {

FlowLaw::FlowLaw(const std::string &prefix, const Config &config,
                 EnthalpyConverter::Ptr ec)
  : m_EC(ec), m_e(1) {

  if (not m_EC) {
    throw RuntimeError(PISM_ERROR_LOCATION, "EC is NULL in FlowLaw::FlowLaw()");
  }

  m_standard_gravity   = config.get_double("constants.standard_gravity");
  m_ideal_gas_constant = config.get_double("constants.ideal_gas_constant");

  m_rho                = config.get_double("constants.ice.density");
  m_beta_CC_grad       = config.get_double("constants.ice.beta_Clausius_Clapeyron") * m_rho * m_standard_gravity;
  m_melting_point_temp = config.get_double("constants.fresh_water.melting_point_temperature");
  m_e                  = config.get_double(prefix + "enhancement_factor");
  m_e_interglacial     = config.get_double(prefix + "enhancement_factor_interglacial");
  m_n                  = config.get_double(prefix + "Glen_exponent");
  m_viscosity_power    = (1.0 - m_n) / (2.0 * m_n);
  m_hardness_power     = -1.0 / m_n;

  m_A_cold = config.get_double("flow_law.Paterson_Budd.A_cold");
  m_A_warm = config.get_double("flow_law.Paterson_Budd.A_warm");
  m_Q_cold = config.get_double("flow_law.Paterson_Budd.Q_cold");
  m_Q_warm = config.get_double("flow_law.Paterson_Budd.Q_warm");
  m_crit_temp = config.get_double("flow_law.Paterson_Budd.T_critical");
  m_schoofLen = config.get_double("flow_law.Schoof_regularizing_length", "m"); // convert to meters
  m_schoofVel = config.get_double("flow_law.Schoof_regularizing_velocity", "m second-1"); // convert to m second-1
  m_schoofReg = PetscSqr(m_schoofVel/m_schoofLen);
}

FlowLaw::~FlowLaw() {
  // empty
}

std::string FlowLaw::name() const {
  return m_name;
}

EnthalpyConverter::Ptr FlowLaw::EC() const {
  return m_EC;
}

double FlowLaw::exponent() const {
  return m_n;
}

double FlowLaw::enhancement_factor() const {
  return m_e;
}

double FlowLaw::enhancement_factor_interglacial() const {
  return m_e_interglacial;
}

//! Return the softness parameter A(T) for a given temperature T.
/*! This is not a natural part of all FlowLaw instances.   */
double FlowLaw::softness_paterson_budd(double T_pa) const {
  const double A = T_pa < m_crit_temp ? m_A_cold : m_A_warm;
  const double Q = T_pa < m_crit_temp ? m_Q_cold : m_Q_warm;

  return A * exp(-Q / (m_ideal_gas_constant * T_pa));
}

//! The flow law itself.
double FlowLaw::flow(double stress, double enthalpy,
                     double pressure, double gs) const {
  return this->flow_impl(stress, enthalpy, pressure, gs);
}

double FlowLaw::flow_impl(double stress, double enthalpy,
                          double pressure, double /* gs */) const {
  return softness(enthalpy, pressure) * pow(stress, m_n-1);
}

void FlowLaw::flow_n(const double *stress, const double *enthalpy,
                     const double *pressure, const double *grainsize,
                     unsigned int n, double *result) const {
  this->flow_n_impl(stress, enthalpy, pressure, grainsize, n, result);
}

void FlowLaw::flow_n_impl(const double *stress, const double *enthalpy,
                          const double *pressure, const double *grainsize,
                          unsigned int n, double *result) const {
  for (unsigned int k = 0; k < n; ++k) {
    result[k] = this->flow(stress[k], enthalpy[k], pressure[k], grainsize[k]);
  }
}


double FlowLaw::softness(double E, double p) const {
  return this->softness_impl(E, p);
}

double FlowLaw::hardness(double E, double p) const {
  return this->hardness_impl(E, p);
}

void FlowLaw::hardness_n(const double *enthalpy, const double *pressure,
                         unsigned int n, double *result) const {
  this->hardness_n_impl(enthalpy, pressure, n, result);
}

void FlowLaw::hardness_n_impl(const double *enthalpy, const double *pressure,
                              unsigned int n, double *result) const {
  for (unsigned int k = 0; k < n; ++k) {
    result[k] = this->hardness(enthalpy[k], pressure[k]);
  }
}

double FlowLaw::hardness_impl(double E, double p) const {
  return pow(softness(E, p), m_hardness_power);
}

//! \brief Computes the regularized effective viscosity and its derivative with respect to the
//! second invariant \f$ \gamma \f$.
/*!
 *
 * @f{align*}{
 * \nu &= \frac{1}{2} B \left( \epsilon + \gamma \right)^{(1-n)/(2n)},\\
 * \diff{\nu}{\gamma} &= \frac{1}{2} B \cdot \frac{1-n}{2n} \cdot \left(\epsilon + \gamma \right)^{(1-n)/(2n) - 1}, \\
 * &= \frac{1-n}{2n} \cdot \frac{1}{2} B \left( \epsilon + \gamma \right)^{(1-n)/(2n)} \cdot \frac{1}{\epsilon + \gamma}, \\
 * &= \frac{1-n}{2n} \cdot \frac{\nu}{\epsilon + \gamma}.
 * @f}
 * Here @f$ \gamma @f$ is the second invariant
 * @f{align*}{
 * \gamma &= \frac{1}{2} D_{ij} D_{ij}\\
 * &= \frac{1}{2}\, ((u_x)^2 + (v_y)^2 + (u_x + v_y)^2 + \frac{1}{2}\, (u_y + v_x)^2) \\
 * @f}
 * and
 * @f[ D_{ij}(\mathbf{u}) = \frac{1}{2}\left(\diff{u_{i}}{x_{j}} + \diff{u_{j}}{x_{i}}\right). @f]
 *
 * Either one of \c nu and \c dnu can be NULL if the corresponding output is not needed.
 *
 * \param[in] B ice hardness
 * \param[in] gamma the second invariant
 * \param[out] nu effective viscosity
 * \param[out] dnu derivative of \f$ \nu \f$ with respect to \f$ \gamma \f$
 */
void FlowLaw::effective_viscosity(double B, double gamma,
                                  double *nu, double *dnu) const {
  const double
    my_nu = 0.5 * B * pow(m_schoofReg + gamma, m_viscosity_power);

  if (PetscLikely(nu != NULL)) {
    *nu = my_nu;
  }

  if (PetscLikely(dnu != NULL)) {
    *dnu = m_viscosity_power * my_nu / (m_schoofReg + gamma);
  }
}

void averaged_hardness_vec(const FlowLaw &ice,
                           const IceModelVec2S &thickness,
                           const IceModelVec3  &enthalpy,
                           IceModelVec2S &result) {

  const IceGrid &grid = *thickness.grid();

  IceModelVec::AccessList list{&thickness, &result, &enthalpy};

  ParallelSection loop(grid.com);
  try {
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // Evaluate column integrals in flow law at every quadrature point's column
      double H = thickness(i,j);
      const double *enthColumn = enthalpy.get_column(i, j);
      result(i,j) = averaged_hardness(ice, H, grid.kBelowHeight(H),
                                      &(grid.z()[0]), enthColumn);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result.update_ghosts();
}

//! Computes vertical average of `B(E, p)` ice hardness, namely @f$\bar B(E, p)@f$.
/*!
 * See comment for hardness(). Note `E[0], ..., E[kbelowH]` must be valid.
 */
double averaged_hardness(const FlowLaw &ice,
                         double thickness, int kbelowH,
                         const double *zlevels,
                         const double *enthalpy) {
  double B = 0;

  EnthalpyConverter &EC = *ice.EC();

  // Use trapezoidal rule to integrate from 0 to zlevels[kbelowH]:
  if (kbelowH > 0) {
    double
      p0 = EC.pressure(thickness),
      E0 = enthalpy[0],
      h0 = ice.hardness(E0, p0); // ice hardness at the left endpoint

    for (int i = 1; i <= kbelowH; ++i) { // note the "1" and the "<="
      const double
        p1 = EC.pressure(thickness - zlevels[i]), // pressure at the right endpoint
        E1 = enthalpy[i], // enthalpy at the right endpoint
        h1 = ice.hardness(E1, p1); // ice hardness at the right endpoint

      // The trapezoid rule sans the "1/2":
      B += (zlevels[i] - zlevels[i-1]) * (h0 + h1);

      h0 = h1;
    }
  }

  // Add the "1/2":
  B *= 0.5;

  // use the "rectangle method" to integrate from
  // zlevels[kbelowH] to thickness:
  double
    depth = thickness - zlevels[kbelowH],
    p = EC.pressure(depth);

  B += depth * ice.hardness(enthalpy[kbelowH], p);

  // Now B is an integral of ice hardness; next, compute the average:
  if (thickness > 0) {
    B = B / thickness;
  } else {
    B = 0;
  }

  return B;
}

bool FlowLawUsesGrainSize(const FlowLaw &flow_law) {
  static const double gs[] = {1e-4, 1e-3, 1e-2, 1}, s=1e4, E=400000, p=1e6;
  double ref = flow_law.flow(s, E, p, gs[0]);
  for (int i=1; i<4; i++) {
    if (flow_law.flow(s, E, p, gs[i]) != ref) {
      return true;
    }
  }
  return false;
}

} // end of namespace rheology
} // end of namespace pism

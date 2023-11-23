// Copyright (C) 2004-2018, 2020, 2021, 2022 Jed Brown, Ed Bueler, and Constantine Khroulev
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

#ifndef __flowlaws_hh
#define __flowlaws_hh

#include <string>

#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/Vector2d.hh"

namespace pism {

namespace array {
class Scalar;
class Array3D;
}
class Config;

/*!
 * This uses the definition of squared second invariant from Hutter and several others, namely the
 * output is @f$ D^2 = \frac 1 2 D_{ij} D_{ij} @f$ where incompressibility is used to compute
 * @f$ D_{zz}. @f$
 *
 * This is the approximation of the full second invariant corresponding to the shallow shelf
 * approximation. In particular, we assume that @f$ u @f$ and @f$ v @f$ are depth-independent (@f$
 * u_z = v_z = 0 @f$) and neglect horizontal derivatives of the vertical velocity (@f$ w_x = w_y = 0
 * @f$).
 */
static inline double secondInvariant_2D(const Vector2d &U_x, const Vector2d &U_y) {
  const double
    u_x = U_x.u,
    u_y = U_y.u,
    v_x = U_x.v,
    v_y = U_y.v,
    w_z = -(u_x + v_y);         // w_z is computed assuming incompressibility of ice
  return 0.5 * (u_x * u_x + v_y * v_y + w_z * w_z + 0.5*(u_y + v_x)*(u_y + v_x));
}

//! Ice flow laws.
namespace rheology {

//! Abstract class containing the constitutive relation for the flow of ice (of
//! the Paterson-Budd type).
/*!
  This is the interface which most of PISM uses for rheology.

  The current implementation of stress-balance computations in PISM restrict
  possible choices of rheologies to ones that

  - are power laws

  - allow factoring out a temperature- (or enthalpy-) dependent ice hardness
    factor

  - can be represented in the viscosity form

  @note FlowLaw derived classes should implement hardness() in
  terms of softness(). That way in many cases we only need to
  re-implement softness... to turn one flow law into another.
*/
class FlowLaw {
public:
  FlowLaw(const std::string &prefix, const Config &config,
          EnthalpyConverter::Ptr EC);
  virtual ~FlowLaw() = default;

  void effective_viscosity(double hardness, double gamma,
                           double *nu, double *dnu) const;

  void effective_viscosity(double hardness, double gamma, double eps,
                           double *nu, double *dnu) const;

  std::string name() const;
  double exponent() const;

  EnthalpyConverter::Ptr EC() const;

  double hardness(double E, double p) const;
  void hardness_n(const double *enthalpy, const double *pressure,
                  unsigned int n, double *result) const;

  double softness(double E, double p) const;

  double flow(double stress, double enthalpy, double pressure, double grain_size) const;
  void flow_n(const double *stress, const double *E,
              const double *pressure, const double *grainsize,
              unsigned int n, double *result) const;

protected:
  virtual double flow_impl(double stress, double E,
                           double pressure, double grainsize) const;
  virtual void flow_n_impl(const double *stress, const double *E,
                           const double *pressure, const double *grainsize,
                           unsigned int n, double *result) const;
  virtual double hardness_impl(double E, double p) const;
  virtual void hardness_n_impl(const double *enthalpy, const double *pressure,
                               unsigned int n, double *result) const;
  virtual double softness_impl(double E, double p) const = 0;

protected:
  std::string m_name;

  //! ice density
  double m_rho;
  //! Clausius-Clapeyron gradient
  double m_beta_CC_grad;
  //! melting point temperature (for water, 273.15 K)
  double m_melting_point_temp;

  EnthalpyConverter::Ptr m_EC;

  double softness_paterson_budd(double T_pa) const;

  //! regularization parameter for @f$ \gamma @f$
  double m_schoofReg;

  //! @f$ (1 - n) / (2n) @f$; used to compute viscosity
  double m_viscosity_power;
  //! @f$ - 1 / n @f$; used to compute hardness
  double m_hardness_power;

  //! Paterson-Budd softness, cold case
  double m_A_cold;
  //! Paterson-Budd softness, warm case
  double m_A_warm;
  //! Activation energy, cold case
  double m_Q_cold;
  //! Activation energy, warm case
  double m_Q_warm;
  //! critical temperature (cold -- warm transition)
  double m_crit_temp;

  //! acceleration due to gravity
  double m_standard_gravity;
  //! ideal gas constant
  double m_ideal_gas_constant;
  //! power law exponent
  double m_n;
};

double averaged_hardness(const FlowLaw &ice,
                         double ice_thickness,
                         unsigned int kbelowH,
                         const double *zlevels,
                         const double *enthalpy);

void averaged_hardness_vec(const FlowLaw &ice,
                           const array::Scalar &ice_thickness,
                           const array::Array3D  &enthalpy,
                           array::Scalar &result);

bool FlowLawUsesGrainSize(const FlowLaw &flow_law);

} // end of namespace rheology
} // end of namespace pism

#endif // __flowlaws_hh

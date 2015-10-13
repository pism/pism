/* Copyright (C) 2015 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _GPBLD_OPTIMIZED_H_
#define _GPBLD_OPTIMIZED_H_

#include <string>

#include "base/enthalpyConverter.hh"

namespace pism {

class Config;
class IceModelVec2S;
class IceModelVec3;

namespace rheology {

class GPBLD_Optimized {
public:
  GPBLD_Optimized(const std::string &prefix,
                  const Config &config,
                  EnthalpyConverter::Ptr EC);

  inline void effective_viscosity(double hardness, double gamma,
                                  double *nu, double *dnu) const;

  double averaged_hardness(double thickness,
                           int kbelowH,
                           const double *zlevels,
                           const double *enthalpy) const;

  void averaged_hardness_vec(const IceModelVec2S &thickness,
                             const IceModelVec3& enthalpy,
                             IceModelVec2S &hardav) const;

  std::string name() const;
  inline double exponent() const;
  inline double enhancement_factor() const;

  double hardness_parameter(double E, double p) const;
  double softness_parameter(double E, double p) const;
  double flow(double stress, double E,
              double pressure, double grainsize) const;
private:
  std::string m_name;

  double m_rho,
    m_beta_CC_grad,
    m_melting_point_temp;
  EnthalpyConverter::Ptr m_EC;

  double softness_parameter_paterson_budd(double T_pa) const;

  double m_schoofLen, m_schoofVel, m_schoofReg, m_viscosity_power,
    m_hardness_power,
    m_A_cold, m_A_warm, m_Q_cold, m_Q_warm,
    m_crit_temp;

  double m_standard_gravity,
    m_ideal_gas_constant,
    m_e,                          // flow enhancement factor
    m_n;                          // power law exponent

  double m_T_0, m_water_frac_coeff, m_water_frac_observed_limit;
};

} // end of namespace rheology
} // end of namespace pism

#endif /* _GPBLD_OPTIMIZED_H_ */

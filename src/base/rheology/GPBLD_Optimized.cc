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

#include <petscsys.h>

#include "GPBLD_Optimized.hh"

#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"

namespace pism {
namespace rheology {

GPBLD_Optimized::GPBLD_Optimized(const std::string &prefix,
                                 const Config &config,
                                 EnthalpyConverter::Ptr EC)
  : m_EC(EC), m_e(1) {

  m_name = "Glen-Paterson-Budd-Lliboutry-Duval (optimized for n == 3)";

  if (not m_EC) {
    throw RuntimeError("EC is NULL in FlowLaw::FlowLaw()");
  }

  m_standard_gravity   = config.get_double("standard_gravity");
  m_ideal_gas_constant = config.get_double("ideal_gas_constant");

  m_rho                = config.get_double("ice_density");
  m_beta_CC_grad       = config.get_double("beta_CC") * m_rho * m_standard_gravity;
  m_melting_point_temp = config.get_double("water_melting_point_temperature");
  m_e                  = config.get_double(prefix + "enhancement_factor");
  m_n                  = config.get_double(prefix + "Glen_exponent");

  if (m_n != 3.0) {
    throw RuntimeError::formatted("GPBLD_Optimized does not support n=%3.3f", m_n);
  }

  m_viscosity_power    = -1.0 / 3.0;
  m_hardness_power     = -1.0 / 3.0;

  m_A_cold = config.get_double("Paterson_Budd_A_cold");
  m_A_warm = config.get_double("Paterson_Budd_A_warm");
  m_Q_cold = config.get_double("Paterson_Budd_Q_cold");
  m_Q_warm = config.get_double("Paterson_Budd_Q_warm");
  m_crit_temp = config.get_double("Paterson_Budd_critical_temperature");
  m_schoofLen = config.get_double("Schoof_regularizing_length", "m"); // convert to meters
  m_schoofVel = config.get_double("Schoof_regularizing_velocity", "m/s"); // convert to m/s
  m_schoofReg = PetscSqr(m_schoofVel/m_schoofLen);

  m_T_0              = config.get_double("water_melting_point_temperature");    // K
  m_water_frac_coeff = config.get_double("gpbld_water_frac_coeff");
  m_water_frac_observed_limit
    = config.get_double("gpbld_water_frac_observed_limit");
}

double GPBLD_Optimized::averaged_hardness(double thickness,
                                          int kbelowH,
                                          const double *zlevels,
                                          const double *enthalpy) const {

}

void GPBLD_Optimized::averaged_hardness_vec(const IceModelVec2S &thickness,
                                            const IceModelVec3& enthalpy,
                                            IceModelVec2S &hardav) const {

}

std::string name() const;
double exponent() const;
inline double enhancement_factor() const;


double GPBLD_Optimized::hardness_parameter(double E, double p) const {

}

double GPBLD_Optimized::softness_parameter(double E, double p) const {

}

double GPBLD_Optimized::flow(double stress, double E,
                             double pressure, double grainsize) const {

}

double GPBLD_Optimized::softness_parameter_paterson_budd(double T_pa) const {

}

} // end of namespace rheology
} // end of namespace pism

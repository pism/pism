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

#include <cmath>
#include <petscsys.h>

#include "GPBLD3.hh"

#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"

namespace pism {
namespace rheology {

GPBLD3::GPBLD3(const std::string &prefix,
               const Config &config,
               EnthalpyConverter::Ptr EC)
  : m_EC(EC), m_e(1) {

  if (not m_EC) {
    throw RuntimeError("EC is NULL in FlowLaw::FlowLaw()");
  }

  m_ideal_gas_constant = config.get_double("ideal_gas_constant");
  m_e                  = config.get_double(prefix + "enhancement_factor");
  m_n                  = config.get_double(prefix + "Glen_exponent");

  if (m_n != 3.0) {
    throw RuntimeError::formatted("GPBLD3 does not support n=%3.3f", m_n);
  }

  m_A_cold    = config.get_double("Paterson_Budd_A_cold");
  m_A_warm    = config.get_double("Paterson_Budd_A_warm");
  m_Q_cold    = config.get_double("Paterson_Budd_Q_cold");
  m_Q_warm    = config.get_double("Paterson_Budd_Q_warm");
  m_crit_temp = config.get_double("Paterson_Budd_critical_temperature");
  m_schoofLen = config.get_double("Schoof_regularizing_length", "m"); // convert to meters
  m_schoofVel = config.get_double("Schoof_regularizing_velocity", "m/s"); // convert to m/s
  m_schoofReg = pow(m_schoofVel/m_schoofLen, 2.0);

  m_T_0              = config.get_double("water_melting_point_temperature");    // K
  m_water_frac_coeff = config.get_double("gpbld_water_frac_coeff");
  m_water_frac_observed_limit
    = config.get_double("gpbld_water_frac_observed_limit");
}

std::string GPBLD3::name() const {
  return "Glen-Paterson-Budd-Lliboutry-Duval (optimized for n == 3)";
}

EnthalpyConverter::Ptr GPBLD3::EC() const {
  return m_EC;
}

double GPBLD3::exponent() const {
  return m_n;
}

double GPBLD3::enhancement_factor() const {
  return m_e;
}

void GPBLD3::effective_viscosity(double hardness, double gamma,
                                 double *nu, double *dnu) const {
  const double
    my_nu = 0.5 * hardness / cbrt(m_schoofReg + gamma);

  if (PetscLikely(nu != NULL)) {
    *nu = my_nu;
  }

  if (PetscLikely(dnu != NULL)) {
    *dnu = (-1.0/3.0) * my_nu / (m_schoofReg + gamma);
  }
}

double GPBLD3::hardness_parameter(double E, double p) const {
  return 1.0 / cbrt(softness_parameter(E, p));
}

double GPBLD3::softness_parameter(double enthalpy, double pressure) const {
  const double E_s = m_EC->enthalpy_cts(pressure);
  if (PetscLikely(enthalpy < E_s)) {       // cold ice
    double T_pa = m_EC->pressure_adjusted_temperature(enthalpy, pressure);
    return softness_parameter_paterson_budd(T_pa);
  } else { // temperate ice
    double omega = m_EC->water_fraction(enthalpy, pressure);
    // as stated in \ref AschwandenBuelerBlatter, cap omega at max of observations:
    omega = std::min(omega, m_water_frac_observed_limit);
    // next line implements eqn (23) in \ref AschwandenBlatter2009
    return softness_parameter_paterson_budd(m_T_0) * (1.0 + m_water_frac_coeff * omega);
  }
}

double GPBLD3::flow(double stress, double enthalpy,
                    double pressure, double grainsize) const {
  (void) grainsize;

  return softness_parameter(enthalpy, pressure) * stress * stress;
}

double GPBLD3::softness_parameter_paterson_budd(double T_pa) const {
  if (PetscLikely(T_pa < m_crit_temp)) {
    return m_A_cold * exp(-m_Q_cold/(m_ideal_gas_constant * T_pa));
  }
  return m_A_warm * exp(-m_Q_warm/(m_ideal_gas_constant * T_pa));
}

} // end of namespace rheology
} // end of namespace pism

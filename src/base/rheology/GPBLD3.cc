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

GPBLD3::GPBLD3(const std::string &prefix, const Config &config, EnthalpyConverter::Ptr EC)
  : GPBLD(prefix, config, EC) {
  m_name = "Glen-Paterson-Budd-Lliboutry-Duval (optimized for n == 3)";

  if (this->exponent() != 3.0) {
    throw RuntimeError::formatted("Can't use GPBLD3 with Glen exponent %f",
                                  this->exponent());
  }
  m_softness_T0 = softness_paterson_budd(m_T_0);
}

double GPBLD3::hardness_impl(double E, double p) const {
  return 1.0 / cbrt(softness_impl(E, p));
}

double GPBLD3::softness_paterson_budd(double T_pa) const {
  const double A = T_pa < m_crit_temp ? m_A_cold : m_A_warm;
  const double Q = T_pa < m_crit_temp ? m_Q_cold : m_Q_warm;

  return A * exp(-Q / (m_ideal_gas_constant * T_pa));
}

double GPBLD3::softness_impl(double enthalpy, double pressure) const {
  const double E_s = m_EC->enthalpy_cts(pressure);
  if (PetscLikely(enthalpy < E_s)) {       // cold ice
    double T_pa = m_EC->pressure_adjusted_temperature(enthalpy, pressure);
    return softness_paterson_budd(T_pa);
  } else { // temperate ice
    double omega = m_EC->water_fraction(enthalpy, pressure);
    // as stated in \ref AschwandenBuelerBlatter, cap omega at max of observations:
    omega = std::min(omega, m_water_frac_observed_limit);
    // next line implements eqn (23) in \ref AschwandenBlatter2009
    return m_softness_T0 * (1.0 + m_water_frac_coeff * omega);
  }
}

double GPBLD3::flow_impl(double stress, double enthalpy, double pressure, double grainsize) const {
  (void) grainsize;

  return softness_impl(enthalpy, pressure) * stress * stress;
}

} // end of namespace rheology
} // end of namespace pism

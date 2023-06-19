/* Copyright (C) 2015, 2023 PISM Authors
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

#include "pism/rheology/PatersonBudd.hh"
#include <cmath>   // for pow
#include <memory>  // for __shared_ptr_access

namespace pism {
namespace rheology {

// PatersonBudd

PatersonBudd::PatersonBudd(const std::string &prefix,
                           const Config &config,
                           EnthalpyConverter::Ptr ec)
  : FlowLaw(prefix, config, ec) {
  m_name = "Paterson-Budd";
}

/*! Converts enthalpy to temperature and uses the Paterson-Budd formula. */
double PatersonBudd::softness_impl(double E, double pressure) const {
  double T_pa = m_EC->pressure_adjusted_temperature(E, pressure);
  return softness_from_temp(T_pa);
}

/*! Converts enthalpy to temperature and calls flow_from_temp. */
double PatersonBudd::flow_impl(double stress, double E,
                               double pressure, double gs) const {
  double temp = m_EC->temperature(E, pressure);
  return flow_from_temp(stress, temp, pressure, gs);
}

//! The flow law (temperature-dependent version).
double PatersonBudd::flow_from_temp(double stress, double temp,
                                    double pressure, double /*gs*/) const {
  // pressure-adjusted temperature:
  const double T_pa = temp + (m_beta_CC_grad / (m_rho * m_standard_gravity)) * pressure;
  return softness_from_temp(T_pa) * pow(stress, m_n-1);
}

double PatersonBudd::softness_from_temp(double T_pa) const {
  return softness_paterson_budd(T_pa);
}

double PatersonBudd::hardness_from_temp(double T_pa) const {
  return pow(softness_from_temp(T_pa), m_hardness_power);
}

} // end of namespace rheology
} // end of namespace pism

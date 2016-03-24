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

#include "GPBLD3.hh"

#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"

#include "approximate/gpbld.hh"
#include "approximate/enthalpy_converter.h"

namespace pism {
namespace rheology {

/*! Compare provided enthalpy converter to the one hard-wired in the implementation.
 */
static void check_enthalpy_converter(EnthalpyConverter::Ptr EC,
                                     const Config &config) {
  struct enth_constants c = enth_get_constants();

  if (c.T_melting != config.get_double("water_melting_point_temperature")) {
    throw RuntimeError("water_melting_point_temperature mismatch");
  }

  if (c.c_i != config.get_double("ice_specific_heat_capacity")) {
    throw RuntimeError("ice_specific_heat_capacity mismatch");
  }

  if (c.c_w != config.get_double("water_specific_heat_capacity")) {
    throw RuntimeError("water_specific_heat_capacity mismatch");
  }

  if (c.T_0 != config.get_double("enthalpy_converter_reference_temperature")) {
    throw RuntimeError("enthalpy_converter_reference_temperature mismatch");
  }

  if (c.beta != config.get_double("beta_CC")) {
    throw RuntimeError("beta_CC mismatch");
  }

  if (c.L0 != config.get_double("water_latent_heat_fusion")) {
    throw RuntimeError("water_latent_heat_fusion mismatch");
  }

  // check that L depends on T_m
  const double depth = 5000.0;
  double T_m = EC->melting_temperature(EC->pressure(depth));
  if (EC->L(T_m) != enth_L(T_m)) {
    throw RuntimeError("selected enthalpy converter is not compatible with GPBLD3: "
                       "different parameterizations of the latent heat of fusion");
  }
}

GPBLD3::GPBLD3(const std::string &prefix, const Config &config, EnthalpyConverter::Ptr ec)
  : GPBLD(prefix, config, ec) {
  m_name = "Glen-Paterson-Budd-Lliboutry-Duval (using a polynomial approximation, optimized for n == 3)";

  if (this->exponent() != 3.0) {
    throw RuntimeError::formatted("Can't use GPBLD3 with Glen exponent %f",
                                  this->exponent());
  }

  struct gpbld_constants c = gpbld_get_constants();

  if (c.ideal_gas_constant != m_ideal_gas_constant) {
    throw RuntimeError::formatted("ideal_gas_constant mismatch: %f != %f",
                                  c.ideal_gas_constant, m_ideal_gas_constant);
  }

  if (c.A_cold != m_A_cold) {
    throw RuntimeError::formatted("A_cold mismatch: %f != %f",
                                  c.A_cold, m_A_cold);
  }

  if (c.A_warm != m_A_warm) {
    throw RuntimeError::formatted("A_warm mismatch: %f != %f",
                                  c.A_warm, m_A_warm);
  }

  if (c.Q_cold != m_Q_cold) {
    throw RuntimeError::formatted("Q_cold mismatch: %f != %f",
                                  c.Q_cold, m_Q_cold);
  }

  if (c.Q_warm != m_Q_warm) {
    throw RuntimeError::formatted("Q_warm mismatch: %f != %f",
                                  c.Q_warm, m_Q_warm);
  }

  if (c.T_critical != m_crit_temp) {
    throw RuntimeError::formatted("T_critical mismatch: %f != %f",
                                  c.T_critical, m_crit_temp);
  }

  if (c.T_melting != m_T_0) {
    throw RuntimeError::formatted("T_melting mismatch: %f != %f",
                                  c.T_melting, m_T_0);
  }

  if (c.water_fraction_coeff != m_water_frac_coeff) {
    throw RuntimeError::formatted("water_fraction_coeff mismatch: %f != %f",
                                  c.water_fraction_coeff, m_water_frac_coeff);
  }

  if (c.water_frac_observed_limit != m_water_frac_observed_limit) {
    throw RuntimeError::formatted("water_frac_observed_limit mismatch: %f != %f",
                                  c.water_frac_observed_limit, m_water_frac_observed_limit);
  }

  check_enthalpy_converter(ec, config);

}

double GPBLD3::hardness_impl(double enthalpy, double pressure) const {
  return gpbld_hardness(enthalpy, pressure);
}

void GPBLD3::hardness_n_impl(const double *enthalpy, const double *pressure,
                             unsigned int n, double *result) const {
  gpbld_hardness_n(enthalpy, pressure, n, result);
}

double GPBLD3::softness_impl(double enthalpy, double pressure) const {
  return gpbld_softness(enthalpy, pressure);
}

double GPBLD3::flow_impl(double stress, double enthalpy, double pressure, double grainsize) const {
  (void) grainsize;

  return gpbld_flow(stress, enthalpy, pressure);
}

void GPBLD3::flow_n_impl(const double *stress, const double *enthalpy,
                         const double *pressure, const double *grainsize,
                         unsigned int n, double *result) const {
  (void) grainsize;

  gpbld_flow_n(stress, enthalpy, pressure, n, result);
}

} // end of namespace rheology
} // end of namespace pism

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

#include "enthalpyConverter.hh"
#include "PISMConfig.hh"
#include "error_handling.hh"

namespace pism {

KirchhoffEnthalpyConverter::KirchhoffEnthalpyConverter(const Config &config)
  : EnthalpyConverter(config) {
  m_c_w = config.get("water_specific_heat_capacity"); // J kg-1 K-1
}

KirchhoffEnthalpyConverter::~KirchhoffEnthalpyConverter() {
  // empty
}

// FIXME: We need a reference for this!
double KirchhoffEnthalpyConverter::L(double T_pm) const {
  return m_L + (m_c_w - m_c_i) * (T_pm - 273.15);
}

double KirchhoffEnthalpyConverter::water_fraction_impl(double E, double pressure) const {

#if (PISM_DEBUG==1)
  if (E >= enthalpy_liquid(pressure)) {
    throw RuntimeError::formatted("E=%f and pressure=%f correspond to liquid water",
                                  E, pressure);
  }
#endif

  double E_s = enthalpy_cts(pressure);
  if (E <= E_s) {
    return 0.0;
  } else {
    return (E - E_s) / L(melting_temperature(pressure));
  }
}

double KirchhoffEnthalpyConverter::enthalpy_impl(double T, double omega, double pressure) const {
  const double T_m = melting_temperature(pressure);

#if (PISM_DEBUG==1)
  if (T <= 0.0) {
    throw RuntimeError::formatted("T = %f <= 0 is not a valid absolute temperature",T);
  }
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    throw RuntimeError::formatted("water fraction omega=%f not in range [0,1]",omega);
  }
  if (T > T_m + 1.0e-6) {
    throw RuntimeError::formatted("T=%f exceeds T_melting=%f; not allowed", T, T_m);
  }
  if ((T < T_m - 1.0e-6) && (omega > 0.0 + 1.0e-6)) {
    throw RuntimeError::formatted("T < T_m AND omega > 0 is contradictory;"
                                  " got T=%f, T_melting=%f, omega=%f",
                                  T, T_m, omega);
  }
#endif

  if (T < T_m) {
    return m_c_i * (T - m_T_0);
  } else {
    return enthalpy_cts(pressure) + omega * L(T_m);
  }
}

double KirchhoffEnthalpyConverter::enthalpy_liquid_impl(double pressure) const {
  return enthalpy_cts(pressure) + L(melting_temperature(pressure));
}

} // end of namespace pism

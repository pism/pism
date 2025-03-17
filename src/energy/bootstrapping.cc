/* Copyright (C) 2016, 2017, 2023, 2025 PISM Authors
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

#include <algorithm>                       // for min
#include <cmath>                           // for sqrt, erf, M_PI
#include <memory>                          // for __shared_ptr_access
#include "pism/energy/bootstrapping.hh"    // for ice_temperature_guess, ice...
#include "pism/util/EnthalpyConverter.hh"  // for EnthalpyConverter, Enthalp...

namespace pism {
namespace energy {

double ice_temperature_guess(EnthalpyConverter &EC,
                             double H, double z, double T_surface,
                             double G, double ice_k) {

  const double
    depth = H - z,
    d2    = depth * depth,
    Tpmp  = EC.melting_temperature(EC.pressure(depth));

  const double
    beta = (4.0/21.0) * (G / (2.0 * ice_k * H * H * H)),
    alpha = (G / (2.0 * H * ice_k)) - 2.0 * H * H * beta;

  return std::min(Tpmp, T_surface + alpha * d2 + beta * d2 * d2);
}

double ice_temperature_guess_smb(EnthalpyConverter &EC, double H, double z, double T_surface,
                                 double G, double ice_k, double K, double SMB) {
  const double depth = H - z, Tpmp = EC.melting_temperature(EC.pressure(depth));

  if (SMB <= 0.0) {
    // negative or zero surface mass balance: linear temperature profile
    return std::min(Tpmp, G / ice_k * depth + T_surface);
  }

  // positive surface mass balance
  const double C0     = (G * sqrt(M_PI * H * K)) / (ice_k * sqrt(2.0 * SMB)),
               gamma0 = sqrt(SMB * H / (2.0 * K));

  return std::min(Tpmp, T_surface + C0 * (erf(gamma0) - erf(gamma0 * z / H)));
}

} // end of namespace energy
} // end of namespace pism

/* Copyright (C) 2015, 2016, 2017, 2018, 2023 PISM Authors
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

#include "pism/rheology/GPBLD.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace rheology {

/*!
  This constructor just sets flow law factor for nonzero water content, from
  \ref AschwandenBlatter and \ref LliboutryDuval1985.
*/
GPBLD::GPBLD(const std::string &prefix,
             const Config &config, EnthalpyConverter::Ptr ec)
  : FlowLaw(prefix, config, ec) {
  m_name = "Glen-Paterson-Budd-Lliboutry-Duval";

  m_T_0              = config.get_number("constants.fresh_water.melting_point_temperature"); // K
  m_water_frac_coeff = config.get_number("flow_law.gpbld.water_frac_coeff");

  m_water_frac_observed_limit = config.get_number("flow_law.gpbld.water_frac_observed_limit");
}

//! The softness factor in the Glen-Paterson-Budd-Lliboutry-Duval flow law.  For constitutive law form.
/*!
  This is a modification of Glen-Paterson-Budd ice, which is PatersonBudd.  In particular, if
  \f$A()\f$ is the softness factor for PatersonBudd, if \f$E\f$ is the enthalpy, and \f$p\f$ is
  the pressure then the softness we compute is
  \f[A = A(T_{pa}(E, p))(1+184\omega).\f]
  The pressure-melting temperature \f$T_{pa}(E, p)\f$ is computed by pressure_adjusted_temperature().
*/
double GPBLD::softness_impl(double enthalpy, double pressure) const {
  const double E_s = m_EC->enthalpy_cts(pressure);
  if (enthalpy < E_s) {       // cold ice
    double T_pa = m_EC->pressure_adjusted_temperature(enthalpy, pressure);
    return softness_paterson_budd(T_pa);
  } else { // temperate ice
    double omega = m_EC->water_fraction(enthalpy, pressure);
    // as stated in \ref AschwandenBuelerBlatter, cap omega at max of observations:
    omega = std::min(omega, m_water_frac_observed_limit);
    // next line implements eqn (23) in \ref AschwandenBlatter2009
    return softness_paterson_budd(m_T_0) * (1.0 + m_water_frac_coeff * omega);
  }
}

void GPBLD::flow_n_impl(const double *stress, const double *enthalpy,
                        const double *pressure, const double *grainsize,
                        unsigned int n, double *result) const {
  // optimize the common case of Glen n=3
  if (m_n == 3.0) {
    for (unsigned int k = 0; k < n; ++k) {
      result[k] = this->softness(enthalpy[k], pressure[k]) * (stress[k] * stress[k]);
    }

    return;
  }

  for (unsigned int k = 0; k < n; ++k) {
    result[k] = this->flow(stress[k], enthalpy[k], pressure[k], grainsize[k]);
  }
}

} // end of namespace rheology
} // end of namespace pism

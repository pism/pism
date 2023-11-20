/* Copyright (C) 2015, 2016, 2017, 2023 PISM Authors
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

#include "pism/rheology/IsothermalGlen.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace rheology {

IsothermalGlen::IsothermalGlen(const std::string &prefix,
                               const Config &config, EnthalpyConverter::Ptr ec)
  : PatersonBudd(prefix, config, ec) {
  m_name = "isothermal Glen";

  m_softness_A = config.get_number("flow_law.isothermal_Glen.ice_softness");
  m_hardness_B = pow(m_softness_A, m_hardness_power);
}

double IsothermalGlen::flow_impl(double stress, double, double, double) const {
  return m_softness_A * pow(stress, m_n-1);
}

double IsothermalGlen::softness_impl(double, double) const {
  return m_softness_A;
}

double IsothermalGlen::hardness_impl(double, double) const {
  return m_hardness_B;
}

double IsothermalGlen::flow_from_temp(double stress, double, double, double) const {
  return m_softness_A * pow(stress,m_n-1);
}

} // end of namespace rheology
} // end of namespace pism

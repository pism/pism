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

#include "PatersonBuddWarm.hh"

namespace pism {
namespace rheology {

PatersonBuddWarm::PatersonBuddWarm(const std::string &prefix,
                   const Config &config, EnthalpyConverter::Ptr ec)
  : PatersonBudd(prefix, config, ec) {
  m_name = "Paterson-Budd (warm case)";
}

double PatersonBuddWarm::tempFromSoftness(double A) const {
  return - m_Q_warm / (m_ideal_gas_constant * (log(A) - log(m_A_warm)));
}

// takes care of hardness...
double PatersonBuddWarm::softness_from_temp(double T_pa) const {
  return m_A_warm * exp(-m_Q_warm / (m_ideal_gas_constant * T_pa));
}

// ignores pressure and uses non-pressure-adjusted temperature
double PatersonBuddWarm::flow_from_temp(double stress, double temp,
                                        double , double) const {
  return softness_from_temp(temp) * pow(stress,m_n-1);
}


} // end of namespace rheology
} // end of namespace pism

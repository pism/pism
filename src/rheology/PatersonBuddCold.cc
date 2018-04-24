/* Copyright (C) 2015, 2016, 2018 PISM Authors
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

#include "PatersonBuddCold.hh"

namespace pism {
namespace rheology {

PatersonBuddCold::PatersonBuddCold(const std::string &prefix,
                                   const Config &config,
                                   EnthalpyConverter::Ptr ec)
  : PatersonBudd(prefix, config, ec) {
  m_name = "Paterson-Budd (cold case)";
}

double PatersonBuddCold::tempFromSoftness(double A) const {
  return - m_Q_cold / (m_ideal_gas_constant * (log(A) - log(m_A_cold)));
}

PatersonBuddCold::~PatersonBuddCold() {
  // empty
}

// takes care of hardness...
double PatersonBuddCold::softness_from_temp(double T_pa) const {
  return m_A_cold * exp(-m_Q_cold / (m_ideal_gas_constant * T_pa));
}

// ignores pressure and uses non-pressure-adjusted temperature
double PatersonBuddCold::flow_from_temp(double stress, double temp,
                                        double , double) const {
  return softness_from_temp(temp) * pow(stress,m_n-1);
}


// Rather than make this part of the base class, we just check at some reference values.
bool FlowLawIsPatersonBuddCold(const FlowLaw &flow_law, const Config &config,
                               EnthalpyConverter::Ptr EC) {
  static const struct {double s, E, p, gs;} v[] = {
    {1e3, 223, 1e6, 1e-3}, {450000, 475000, 500000, 525000}, {5e4, 268, 5e6, 3e-3}, {1e5, 273, 8e6, 5e-3}};
  PatersonBuddCold cpb("stress_balance.sia.", config, EC); // This is unmodified cold Paterson-Budd
  for (int i=0; i<4; i++) {
    const double left  = flow_law.flow(v[i].s, v[i].E, v[i].p, v[i].gs),
      right =  cpb.flow(v[i].s, v[i].E, v[i].p, v[i].gs);
    if (fabs((left - right)/left)>1.0e-15) {
      return false;
    }
  }
  return true;
}

} // end of namespace rheology
} // end of namespace pism

/* Copyright (C) 2018 PISM Authors
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
#include <cmath> // sqrt

#include "FrontalMeltPhysics.hh"

#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace frontalmelt {

FrontalMeltPhysics::FrontalMeltPhysics(const Config &config) {

  m_alpha         = config.get_double("frontal_melt.power_alpha");
  m_beta          = config.get_double("frontal_melt.power_beta");
  std::string a_units  = "m−" + std::to_string(m_alpha) + "s^(" + std::to_string(m_alpha) + "−1) Celsius−" + std::to_string(m_beta);
  std::string b_units  = "m s^(" + std::to_string(m_alpha) + "−1) Celsius−" + std::to_string(m_beta);
  m_A             = config.get_double("frontal_melt.parameter_a", a_units);
  m_B             = config.get_double("frontal_melt.parameter_b", b_units);
  m_water_density = config.get_double("constants.fresh_water.density");
  
}

double FrontalMeltPhysics::frontal_melt_from_undercutting(double h, double Qsg, double TF) const {

    return (m_A * h * pow(Qsg, m_alpha) + m_B) * pow(TF, m_beta);

}


} // end of namespace frontalmeltmodel
} // end of namespace pism

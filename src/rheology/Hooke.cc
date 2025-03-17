/* Copyright (C) 2015, 2016, 2017, 2023, 2025 PISM Authors
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

#include "pism/rheology/Hooke.hh"
#include "pism/util/Config.hh"

namespace pism {
namespace rheology {

// Hooke

Hooke::Hooke(const std::string &prefix,
             const Config &config, std::shared_ptr<EnthalpyConverter> ec)
  : PatersonBudd(prefix, config, ec) {
  m_name = "Hooke";

  m_Q_Hooke  = config.get_number("flow_law.Hooke.Q");
  m_A_Hooke  = config.get_number("flow_law.Hooke.A");
  m_C_Hooke  = config.get_number("flow_law.Hooke.C");
  m_K_Hooke  = config.get_number("flow_law.Hooke.k");
  m_Tr_Hooke = config.get_number("flow_law.Hooke.Tr");
}

double Hooke::softness_from_temp(double T_pa) const {
  return m_A_Hooke * exp(-m_Q_Hooke/(m_ideal_gas_constant * T_pa)
                         + 3.0 * m_C_Hooke * pow(m_Tr_Hooke - T_pa, -m_K_Hooke));
}

} // end of namespace rheology
} // end of namespace pism

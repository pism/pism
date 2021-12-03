/* Copyright (C) 2015, 2021 PISM Authors
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

#ifndef _HOOKE_H_
#define _HOOKE_H_

#include "PatersonBudd.hh"

namespace pism {
namespace rheology {

//! The Hooke flow law.
class Hooke : public PatersonBudd {
public:
  Hooke(const std::string &prefix, const Config &config, EnthalpyConverter::Ptr EC);
  virtual ~Hooke() = default;
protected:
  virtual double softness_from_temp(double T_pa) const;

  double m_A_Hooke, m_Q_Hooke, m_C_Hooke, m_K_Hooke, m_Tr_Hooke; // constants from Hooke (1981)
  // R_Hooke is the ideal_gas_constant.
};

} // end of namespace rheology
} // end of namespace pism


#endif /* _HOOKE_H_ */

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

#ifndef _PATERSONBUDDCOLD_H_
#define _PATERSONBUDDCOLD_H_

#include "PatersonBudd.hh"

namespace pism {
namespace rheology {

//! Cold case of Paterson-Budd
class PatersonBuddCold : public PatersonBudd {
public:
  PatersonBuddCold(const std::string &prefix, const Config &config, EnthalpyConverter::Ptr EC);
  virtual ~PatersonBuddCold();

  //! Return the temperature T corresponding to a given value A=A(T).
  double tempFromSoftness(double A) const;

protected:
  // takes care of hardness...
  double softness_from_temp(double T_pa) const;

  // ignores pressure and uses non-pressure-adjusted temperature
  double flow_from_temp(double stress, double temp,
                        double , double) const;
};

bool FlowLawIsPatersonBuddCold(const FlowLaw &flow_law,
                               const Config &, EnthalpyConverter::Ptr);

} // end of namespace rheology
} // end of namespace pism


#endif /* _PATERSONBUDDCOLD_H_ */

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

#ifndef _PATERSONBUDD_H_
#define _PATERSONBUDD_H_

#include "FlowLaw.hh"

namespace pism {
namespace rheology {

//! Derived class of FlowLaw for Paterson-Budd (1982)-Glen ice.
class PatersonBudd : public FlowLaw {
public:
  PatersonBudd(const std::string &prefix,
               const Config &config,
               EnthalpyConverter::Ptr EC);
  virtual ~PatersonBudd() = default;

protected:
  virtual double flow_impl(double stress, double E,
                           double pressure, double gs) const;
  // This also takes care of hardness
  virtual double softness_impl(double enthalpy, double pressure) const;

  virtual double softness_from_temp(double T_pa) const;
  virtual double hardness_from_temp(double T_pa) const;

  // special temperature-dependent method
  virtual double flow_from_temp(double stress, double temp,
                                double pressure, double gs) const;
};

} // end of namespace rheology
} // end of namespace pism


#endif /* _PATERSONBUDD_H_ */

/* Copyright (C) 2015, 2023, 2025 PISM Authors
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

#ifndef _ISOTHERMALGLEN_H_
#define _ISOTHERMALGLEN_H_

#include "pism/rheology/PatersonBudd.hh"

namespace pism {
namespace rheology {

//! Isothermal Glen ice allowing extra customization.
class IsothermalGlen : public PatersonBudd {
public:
  IsothermalGlen(const std::string &prefix, const Config &config, std::shared_ptr<EnthalpyConverter> EC);
protected:
  double flow_impl(double stress, double, double, double) const;
  double softness_impl(double, double) const;
  double hardness_impl(double, double) const;
  double flow_from_temp(double stress, double, double, double) const;
protected:
  double m_softness_A, m_hardness_B;
};

} // end of namespace rheology
} // end of namespace pism


#endif /* _ISOTHERMALGLEN_H_ */

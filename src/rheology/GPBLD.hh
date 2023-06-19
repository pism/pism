/* Copyright (C) 2015, 2023 PISM Authors
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

#ifndef _GPBLD_H_
#define _GPBLD_H_

#include "pism/rheology/FlowLaw.hh"

namespace pism {
namespace rheology {

//! Glen (1955) and Paterson-Budd (1982) flow law with additional water fraction factor from Lliboutry & Duval (1985).
/*!
  See [\ref AschwandenBlatter]. The basic references are [\ref Glen] and [\ref
  PatersonBudd] and [\ref LliboutryDuval1985].
*/
class GPBLD : public FlowLaw {
public:
  GPBLD(const std::string &prefix, const Config &config, EnthalpyConverter::Ptr EC);
protected:
  double softness_impl(double enthalpy, double pressure) const;
  void flow_n_impl(const double *stress, const double *enthalpy,
                   const double *pressure, const double *grainsize,
                   unsigned int n, double *result) const;
  double m_T_0, m_water_frac_coeff, m_water_frac_observed_limit;
};

} // end of namespace rheology
} // end of namespace pism


#endif /* _GPBLD_H_ */

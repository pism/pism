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

#ifndef _GPBLD_OPTIMIZED_H_
#define _GPBLD_OPTIMIZED_H_

#include "GPBLD.hh"

namespace pism {
namespace rheology {

class GPBLD3 : public GPBLD {
public:
  GPBLD3(const std::string &prefix, const Config &config, EnthalpyConverter::Ptr EC);
protected:
  double flow_impl(double stress, double E, double pressure, double grainsize) const;
  void flow_n_impl(const double *stress, const double *E,
                   const double *pressure, const double *grainsize,
                   unsigned int n, double *result) const;
  double hardness_impl(double E, double p) const;
  void hardness_n_impl(const double *enthalpy, const double *pressure,
                       unsigned int n, double *result) const;
  double softness_impl(double E, double p) const;
};

} // end of namespace rheology
} // end of namespace pism

#endif /* _GPBLD_OPTIMIZED_H_ */

/* Copyright (C) 2015, 2017, 2018 PISM Authors
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

#ifndef _REGIONALDEFAULTYIELDSTRESS_H_
#define _REGIONALDEFAULTYIELDSTRESS_H_

#include "pism/basalstrength/MohrCoulombYieldStress.hh"

namespace pism {

class RegionalDefaultYieldStress : public MohrCoulombYieldStress {
public:
  RegionalDefaultYieldStress(IceGrid::ConstPtr g);
  virtual ~RegionalDefaultYieldStress();
protected:
  virtual void init_impl(const Geometry &geometry,
                         const IceModelVec2S &till_water_thickness,
                         const IceModelVec2S &overburden_pressure);
  virtual void update_impl(const YieldStressInputs &inputs);
};

} // end of namespace pism

#endif /* _REGIONALDEFAULTYIELDSTRESS_H_ */

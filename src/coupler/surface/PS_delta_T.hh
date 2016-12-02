// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _PS_DELTA_T_H_
#define _PS_DELTA_T_H_

#include "coupler/util/PScalarForcing.hh"
#include "coupler/PISMSurface.hh"
#include "PSModifier.hh"

namespace pism {
namespace surface {

/** @brief Implements the scalar temperature offsets for the ice
 * surface temperature.
 *
 * Other fields are passed through without change.
 */
class Delta_T : public PScalarForcing<SurfaceModel,SurfaceModifier>
{
public:
  Delta_T(IceGrid::ConstPtr g, SurfaceModel* in);
  virtual ~Delta_T();
protected:
  virtual void init_impl();
  virtual void ice_surface_temperature_impl(IceModelVec2S &result) const;
  virtual MaxTimestep max_timestep_impl(double t) const;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PS_DELTA_T_H_ */

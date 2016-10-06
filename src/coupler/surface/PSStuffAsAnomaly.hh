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

#ifndef _PSSTUFFASANOMALY_H_
#define _PSSTUFFASANOMALY_H_

#include "coupler/PISMSurface.hh"
#include "PSModifier.hh"
#include "base/util/iceModelVec.hh"

namespace pism {
namespace surface {

//! \brief A surface modifier class applying its input as anomalies.
class StuffAsAnomaly : public SurfaceModifier
{
public:
  StuffAsAnomaly(IceGrid::ConstPtr g, SurfaceModel *input);
  virtual ~StuffAsAnomaly();
protected:
  virtual void init_impl();
  virtual void update_impl(double my_t, double my_dt);

  virtual void ice_surface_mass_flux_impl(IceModelVec2S &result) const;
  virtual void ice_surface_temperature_impl(IceModelVec2S &result) const;
  virtual MaxTimestep max_timestep_impl(double t) const;

protected:
  IceModelVec2S m_mass_flux;
  IceModelVec2S m_mass_flux_0;
  IceModelVec2S m_mass_flux_input;
  IceModelVec2S m_temp;
  IceModelVec2S m_temp_0;
  IceModelVec2S m_temp_input;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSSTUFFASANOMALY_H_ */

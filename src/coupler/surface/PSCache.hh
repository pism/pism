/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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

#ifndef _PSCACHE_H_
#define _PSCACHE_H_

#include "PSModifier.hh"
#include "base/util/iceModelVec.hh"

namespace pism {
namespace surface {

class Cache : public SurfaceModifier {
public:
  Cache(IceGrid::ConstPtr g, SurfaceModel* in);
  virtual ~Cache();
protected:
  virtual void init_impl();
  virtual void update_impl(double my_t, double my_dt);

  virtual void surface_layer_thickness_impl(IceModelVec2S &result) const;
  virtual void mass_held_in_surface_layer_impl(IceModelVec2S &result) const;
  virtual void ice_surface_mass_flux_impl(IceModelVec2S &result) const;
  virtual void ice_surface_temperature_impl(IceModelVec2S &result) const;
  virtual void ice_surface_liquid_water_fraction_impl(IceModelVec2S &result) const;

  virtual MaxTimestep max_timestep_impl(double t) const;
protected:
  IceModelVec2S m_mass_flux;
  IceModelVec2S m_temperature;
  IceModelVec2S m_liquid_water_fraction;
  IceModelVec2S m_mass_held_in_surface_layer;
  IceModelVec2S m_surface_layer_thickness;
  double m_next_update_time;
  unsigned int m_update_interval_years;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSCACHE_H_ */

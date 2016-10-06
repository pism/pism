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

#ifndef _POCACHE_H_
#define _POCACHE_H_

#include "POModifier.hh"
#include "base/util/iceModelVec.hh"

namespace pism {
namespace ocean {
class Cache : public OceanModifier {
public:
  Cache(IceGrid::ConstPtr g, OceanModel* in);
  virtual ~Cache();

protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(double my_t, double my_dt);

  virtual void init_impl();
  virtual void melange_back_pressure_fraction_impl(IceModelVec2S &result) const;
  virtual void sea_level_elevation_impl(double &result) const;
  virtual void shelf_base_temperature_impl(IceModelVec2S &result) const;
  virtual void shelf_base_mass_flux_impl(IceModelVec2S &result) const;
protected:
  IceModelVec2S m_shelf_base_temperature, m_shelf_base_mass_flux,
    m_melange_back_pressure_fraction;
  double m_sea_level;
  double m_next_update_time;
  unsigned int m_update_interval_years;
};

} // end of namespace ocean
} // end of namespace pism
#endif /* _POCACHE_H_ */

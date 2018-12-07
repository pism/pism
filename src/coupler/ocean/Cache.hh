/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "pism/coupler/OceanModel.hh"

namespace pism {
namespace ocean {

class Cache : public OceanModel {
public:
  Cache(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> in);
  virtual ~Cache();

protected:
  MaxTimestep max_timestep_impl(double t) const;

  void update_impl(const Geometry &geometry, double my_t, double my_dt);
  void init_impl(const Geometry &geometry);

  const IceModelVec2S& shelf_base_temperature_impl() const;
  const IceModelVec2S& shelf_base_mass_flux_impl() const;
  const IceModelVec2S& melange_back_pressure_fraction_impl() const;
private:
  double m_next_update_time;
  unsigned int m_update_interval_years;

  // storage for melange_back_pressure_fraction is inherited from OceanModel
  IceModelVec2S::Ptr m_shelf_base_temperature;
  IceModelVec2S::Ptr m_shelf_base_mass_flux;
};

} // end of namespace ocean
} // end of namespace pism
#endif /* _POCACHE_H_ */

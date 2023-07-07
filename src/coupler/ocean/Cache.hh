/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2021 PISM Authors
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
  Cache(std::shared_ptr<const Grid> g, std::shared_ptr<OceanModel> in);
  virtual ~Cache() = default;

protected:
  MaxTimestep max_timestep_impl(double t) const;

  void update_impl(const Geometry &geometry, double my_t, double my_dt);
  void init_impl(const Geometry &geometry);

  const array::Scalar& shelf_base_temperature_impl() const;
  const array::Scalar& shelf_base_mass_flux_impl() const;
  const array::Scalar& average_water_column_pressure_impl() const;
private:
  double m_next_update_time;
  double m_update_interval_years;

  // storage for melange_back_pressure_fraction is inherited from OceanModel
  std::shared_ptr<array::Scalar> m_shelf_base_temperature;
  std::shared_ptr<array::Scalar> m_shelf_base_mass_flux;
};

} // end of namespace ocean
} // end of namespace pism
#endif /* _POCACHE_H_ */

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

#ifndef _PSCACHE_H_
#define _PSCACHE_H_

#include "pism/coupler/SurfaceModel.hh"

namespace pism {
namespace surface {

class Cache : public SurfaceModel {
public:
  Cache(std::shared_ptr<const Grid> g, std::shared_ptr<SurfaceModel> in);
  virtual ~Cache() = default;
protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const array::Scalar &layer_mass_impl() const;
  const array::Scalar &liquid_water_fraction_impl() const;
  const array::Scalar &temperature_impl() const;
  const array::Scalar &mass_flux_impl() const;
  const array::Scalar &layer_thickness_impl() const;

  virtual const array::Scalar& accumulation_impl() const;
  virtual const array::Scalar& melt_impl() const;
  virtual const array::Scalar& runoff_impl() const;

  MaxTimestep max_timestep_impl(double t) const;
protected:
  // storage for the rest of the fields is inherited from SurfaceModel
  array::Scalar::Ptr m_mass_flux;
  array::Scalar::Ptr m_temperature;

  double m_next_update_time;
  double m_update_interval_years;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSCACHE_H_ */

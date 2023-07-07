// Copyright (C) 2011, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2022, 2023 Constantine Khroulev
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

#ifndef _PODIRECTFORCING_H_
#define _PODIRECTFORCING_H_

#include "pism/coupler/OceanModel.hh"

namespace pism {
namespace ocean {

class Given : public OceanModel
{
public:
  Given(std::shared_ptr<const Grid> g);
  virtual ~Given() = default;

private:
  void update_impl(const Geometry &geometry, double t, double dt);
  void init_impl(const Geometry &geometry);
  MaxTimestep max_timestep_impl(double t) const;

  const array::Scalar& shelf_base_temperature_impl() const;
  const array::Scalar& shelf_base_mass_flux_impl() const;

  std::shared_ptr<array::Forcing> m_shelfbtemp;
  std::shared_ptr<array::Forcing> m_shelfbmassflux;

  std::shared_ptr<array::Scalar> m_shelf_base_temperature;
  std::shared_ptr<array::Scalar> m_shelf_base_mass_flux;
};

} // end of namespace ocean
} // end of namespace pism
#endif /* _PODIRECTFORCING_H_ */

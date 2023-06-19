// Copyright (C) 2019, 2021, 2023 PISM Authors
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

#ifndef _PSISMIP6_H_
#define _PSISMIP6_H_

#include "pism/coupler/SurfaceModel.hh"

namespace pism {
namespace surface {

class ISMIP6 : public SurfaceModel {
public:
  ISMIP6(std::shared_ptr<const Grid> g, std::shared_ptr<atmosphere::AtmosphereModel> input);
  virtual ~ISMIP6() = default;
protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const array::Scalar &temperature_impl() const;
  const array::Scalar &mass_flux_impl() const;

  const array::Scalar& accumulation_impl() const;
  const array::Scalar& melt_impl() const;
  const array::Scalar& runoff_impl() const;
  MaxTimestep max_timestep_impl(double t) const;
protected:
  // time-dependent inputs
  std::shared_ptr<array::Forcing> m_mass_flux_anomaly;
  std::shared_ptr<array::Forcing> m_temperature_anomaly;
  std::shared_ptr<array::Forcing> m_mass_flux_gradient;
  std::shared_ptr<array::Forcing> m_temperature_gradient;

  // time-independent inputs
  array::Scalar m_mass_flux_reference;
  array::Scalar m_temperature_reference;
  array::Scalar m_surface_reference;

  // outputs; stored as shared_ptr to be able to use SurfaceModel::allocate_xxx()
  array::Scalar::Ptr m_mass_flux;
  array::Scalar::Ptr m_temperature;

};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSISMIP6_H_ */

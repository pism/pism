// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Andy Aschwanden and Constantine Khroulev
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

#ifndef _PSELEVATION_H_
#define _PSELEVATION_H_

#include "pism/coupler/SurfaceModel.hh"

namespace pism {
namespace surface {

//! \brief A class implementing a elevation-dependent temperature and mass balance model.
class Elevation : public SurfaceModel {
public:
  Elevation(std::shared_ptr<const IceGrid> grid, std::shared_ptr<atmosphere::AtmosphereModel> input);

protected:

  virtual const array::Scalar& accumulation_impl() const;
  virtual const array::Scalar& melt_impl() const;
  virtual const array::Scalar& runoff_impl() const;
  
private:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  MaxTimestep max_timestep_impl(double t) const;

  const array::Scalar& mass_flux_impl() const;
  const array::Scalar& temperature_impl() const;

  void compute_mass_flux(const array::Scalar &surface, array::Scalar &result) const;
  void compute_temperature(const array::Scalar &surface, array::Scalar &result) const;

  double m_T_min, m_T_max, m_z_T_min, m_z_T_max;
  double m_M_min, m_M_max, m_M_limit_min, m_M_limit_max, m_z_M_min, m_z_ELA, m_z_M_max;

  array::Scalar::Ptr m_mass_flux;
  array::Scalar::Ptr m_temperature;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSELEVATION_H_ */

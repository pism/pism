// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019 PISM Authors
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

#ifndef _PSCONSTANTPIK_H_
#define _PSCONSTANTPIK_H_

#include "pism/coupler/SurfaceModel.hh"

namespace pism {
namespace surface {

//! \brief A class implementing a constant-in-time surface model for the surface mass balance.
//!
//! Reads data from a PISM input file.
//!
//! Ice surface temperature is parameterized as in PISM-PIK, using a latitude
//! and surface elevation-dependent formula.

class PIK : public SurfaceModel {
public:
  PIK(std::shared_ptr<const Grid> g, std::shared_ptr<atmosphere::AtmosphereModel> atmosphere);
protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const array::Scalar& mass_flux_impl() const;
  const array::Scalar& temperature_impl() const;

  const array::Scalar& accumulation_impl() const;
  const array::Scalar& melt_impl() const;
  const array::Scalar& runoff_impl() const;

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  MaxTimestep max_timestep_impl(double t) const;

protected:
  std::shared_ptr<array::Scalar> m_mass_flux;
  std::shared_ptr<array::Scalar> m_temperature;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSCONSTANTPIK_H_ */

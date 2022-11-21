// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2022 PISM Authors
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

#ifndef PISM_OCEAN_ISMIP6_H
#define PISM_OCEAN_ISMIP6_H

#include "CompleteOceanModel.hh"

#include "pism/util/array/Forcing.hh"

namespace pism {
namespace ocean {
//! \brief Implements the ocean model used for ISMIP6.
//!
//! Uses a parameterization of sub-shelf melting with respect to
//! sub-shelf heat flux like in [@ref Jourdain et al., 2020].
//!
//! Models heat flux into the base of the shelf as
//!
//! @f[ Q_{\text{heat}} = \rho_{o} c_{p_{o}} \gamma_{0} (max((T_{o} - T_{f})+\delta T_{sector},0))**2, @f]
//!
//! where @f$\rho_{o}@f$ is the density of ocean water, @f$c_{p_{o}}@f$ and
//! @f$T_{o}@f$ are the heat capacity and temperature of the ocean mixed
//! layer, @f$T_{f}@f$ is the freezing temperature of ocean water at the
//! shelf bottom.
class ISMIP6 : public CompleteOceanModel {
public:
  ISMIP6(IceGrid::ConstPtr g);
  virtual ~ISMIP6() = default;

private:
  MaxTimestep max_timestep_impl(double t) const;
  void update_impl(const Geometry &geometry, double t, double dt);
  void init_impl(const Geometry &geometry);

  // outputs variables from ISMIP6 routine
  const array::Scalar& shelf_base_temperature_impl() const;
  const array::Scalar& shelf_base_mass_flux_impl() const;

  // Variables to be read from input file
  std::shared_ptr<array::Forcing> m_shelfbtemp;
  std::shared_ptr<array::Forcing> m_salinity_ocean;
  
  void mass_flux(const array::Scalar &ice_thickness,
                 const array::Scalar &m_shelfbtemp,
                 const array::Scalar &m_salinity_ocean,
                 array::Scalar &result) const;
};

} // end of namespace ocean
} // end of namespace pism
#endif /* PISM_OCEAN_ISMIP6_H */

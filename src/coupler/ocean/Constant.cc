// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2022 PISM Authors
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

#include "Constant.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Grid.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace ocean {

Constant::Constant(std::shared_ptr<const Grid> g)
  : CompleteOceanModel(g) {
  // empty
}

void Constant::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  const array::Scalar &ice_thickness = geometry.ice_thickness;

  const double
    melt_rate     = m_config->get_number("ocean.constant.melt_rate", "m second-1"),
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity"),
    mass_flux     = melt_rate * ice_density;

  melting_point_temperature(ice_thickness, *m_shelf_base_temperature);

  m_shelf_base_mass_flux->set(mass_flux);

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                           *m_water_column_pressure);
}

void Constant::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2, "* Initializing the constant ocean model...\n");
  m_log->message(2, "  Sub-shelf melt rate set to %f m/year.\n",
                 m_config->get_number("ocean.constant.melt_rate", "m year-1"));

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                           *m_water_column_pressure);
}

MaxTimestep Constant::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("ocean constant");
}

/*!
 * Compute melting point temperature of ice at the depth `depth`.
 */
void Constant::melting_point_temperature(const array::Scalar& depth,
                                         array::Scalar &result) const {
  const double
    T0          = m_config->get_number("constants.fresh_water.melting_point_temperature"),
    beta_CC     = m_config->get_number("constants.ice.beta_Clausius_Clapeyron"),
    g           = m_config->get_number("constants.standard_gravity"),
    ice_density = m_config->get_number("constants.ice.density");

  array::AccessScope list{&depth, &result};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double pressure = ice_density * g * depth(i, j); // FIXME issue #15

    result(i, j) = T0 - beta_CC * pressure;
  }
}

} // end of namespape ocean
} // end of namespace pism

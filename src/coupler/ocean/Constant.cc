// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "pism/util/Vars.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace ocean {

Constant::Constant(IceGrid::ConstPtr g)
  : CompleteOceanModel(g) {
  // empty
}

Constant::~Constant() {
  // empty
}

void Constant::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  const IceModelVec2S &ice_thickness = geometry.ice_thickness;

  const double
    melt_rate   = m_config->get_double("ocean.constant.melt_rate", "m second-1"),
    ice_density = m_config->get_double("constants.ice.density"),
    mass_flux   = melt_rate * ice_density;

  melting_point_temperature(ice_thickness, *m_shelf_base_temperature);

  m_shelf_base_mass_flux->set(mass_flux);

  m_melange_back_pressure_fraction->set(m_config->get_double("ocean.constant.melange_back_pressure_fraction"));
}

void Constant::init_impl(const Geometry &geometry) {
  (void) geometry;

  if (not m_config->get_boolean("ocean.always_grounded")) {
    m_log->message(2, "* Initializing the constant ocean model...\n");
    m_log->message(2, "  Sub-shelf melt rate set to %f m/year.\n",
                   m_config->get_double("ocean.constant.melt_rate", "m year-1"));

  }
}

MaxTimestep Constant::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("ocean constant");
}

/*!
 * Compute melting point temperature of ice at the depth `depth`.
 */
void Constant::melting_point_temperature(const IceModelVec2S& depth,
                                         IceModelVec2S &result) const {
  const double
    T0          = m_config->get_double("constants.fresh_water.melting_point_temperature"),
    beta_CC     = m_config->get_double("constants.ice.beta_Clausius_Clapeyron"),
    g           = m_config->get_double("constants.standard_gravity"),
    ice_density = m_config->get_double("constants.ice.density");

  IceModelVec::AccessList list{&depth, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double pressure = ice_density * g * depth(i, j); // FIXME issue #15

    result(i, j) = T0 - beta_CC * pressure;
  }
}

} // end of namespape ocean
} // end of namespace pism

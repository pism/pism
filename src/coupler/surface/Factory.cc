/* Copyright (C) 2015, 2016, 2017 PISM Authors
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

#include "Factory.hh"

// surface models:
#include "Anomaly.hh"
#include "Elevation.hh"
#include "GivenClimate.hh"
#include "LapseRates.hh"
#include "StuffAsAnomaly.hh"
#include "Delta_T.hh"
#include "TemperatureIndex.hh"
#include "Simple.hh"
#include "ConstantPIK.hh"
#include "ForceThickness.hh"
#include "Cache.hh"


namespace pism {
namespace surface {
// Surface
Factory::Factory(IceGrid::ConstPtr  g)
  : PCFactory<SurfaceModel,SurfaceModifier>(g) {
  m_option = "surface";

  add_model<Elevation>("elevation");
  add_model<Given>("given");
  add_model<TemperatureIndex>("pdd");
  add_model<PIK>("pik");
  add_model<Simple>("simple");
  set_default("given");

  add_modifier<Anomaly>("anomaly");
  add_modifier<Cache>("cache");
  add_modifier<Delta_T>("delta_T");
  add_modifier<ForceThickness>("forcing");
  add_modifier<LapseRates>("lapse_rate");
  add_modifier<StuffAsAnomaly>("turn_into_anomaly");
}

Factory::~Factory() {
  // empty
}

} // end of namespace surface
} // end of namespace pism

/* Copyright (C) 2015, 2017, 2018, 2019, 2020, 2023 PISM Authors
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

#include "pism/coupler/atmosphere/Factory.hh"

// atmosphere models:
#include "pism/coupler/atmosphere/GivenClimate.hh"
#include "pism/coupler/atmosphere/ElevationChange.hh"
#include "pism/coupler/atmosphere/SeariseGreenland.hh"
#include "pism/coupler/atmosphere/Delta_T.hh"
#include "pism/coupler/atmosphere/Delta_P.hh"
#include "pism/coupler/atmosphere/Frac_P.hh"
#include "pism/coupler/atmosphere/PrecipitationScaling.hh"
#include "pism/coupler/atmosphere/PIK.hh"
#include "pism/coupler/atmosphere/Anomaly.hh"
#include "pism/coupler/atmosphere/CosineYearlyCycle.hh"
#include "pism/coupler/atmosphere/WeatherStation.hh"
#include "pism/coupler/atmosphere/Uniform.hh"
#include "pism/coupler/atmosphere/OrographicPrecipitation.hh"

namespace pism {
namespace atmosphere {

// Atmosphere

Factory::Factory(std::shared_ptr<const Grid> g)
  : PCFactory<AtmosphereModel>(g, "atmosphere.models") {

  add_model<PIK>("pik");
  add_model<Given>("given");
  add_model<SeaRISEGreenland>("searise_greenland");
  add_model<CosineYearlyCycle>("yearly_cycle");
  add_model<WeatherStation>("one_station");
  add_model<Uniform>("uniform");

  add_modifier<Anomaly>("anomaly");
  add_modifier<PrecipitationScaling>("precip_scaling");
  add_modifier<Frac_P>("frac_P");
  add_modifier<Delta_P>("delta_P");
  add_modifier<Delta_T>("delta_T");
  add_modifier<ElevationChange>("elevation_change");
  add_modifier<OrographicPrecipitation>("orographic_precipitation");
}

} // end of namespace atmosphere
} // end of namespace pism

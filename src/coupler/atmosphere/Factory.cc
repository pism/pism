/* Copyright (C) 2015, 2017, 2018 PISM Authors
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

// atmosphere models:
#include "GivenClimate.hh"
#include "LapseRates.hh"
#include "SeariseGreenland.hh"
#include "Delta_T.hh"
#include "Delta_P.hh"
#include "Frac_P.hh"
#include "Paleo_precip.hh"
#include "PIK.hh"
#include "Anomaly.hh"
#include "CosineYearlyCycle.hh"
#include "WeatherStation.hh"
#include "Uniform.hh"

namespace pism {
namespace atmosphere {

// Atmosphere

Factory::Factory(IceGrid::ConstPtr g)
  : PCFactory<AtmosphereModel>(g) {
  m_option = "atmosphere";

  add_model<PIK>("pik");
  add_model<Given>("given");
  add_model<SeaRISEGreenland>("searise_greenland");
  add_model<CosineYearlyCycle>("yearly_cycle");
  add_model<WeatherStation>("one_station");
  add_model<Uniform>("uniform");
  set_default("given");

  add_modifier<Anomaly>("anomaly");
  add_modifier<PaleoPrecip>("paleo_precip");
  add_modifier<Frac_P>("frac_P");
  add_modifier<Delta_P>("delta_P");
  add_modifier<Delta_T>("delta_T");
  add_modifier<LapseRates>("lapse_rate");
}

Factory::~Factory() {
  // empty
}

} // end of namespace atmosphere
} // end of namespace pism

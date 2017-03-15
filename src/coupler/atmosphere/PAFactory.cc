/* Copyright (C) 2015 PISM Authors
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

#include "PAFactory.hh"

// atmosphere models:
#include "PAGivenClimate.hh"
#include "PALapseRates.hh"
#include "PASeariseGreenland.hh"
#include "PA_delta_T.hh"
#include "PA_delta_P.hh"
#include "PA_frac_P.hh"
#include "PA_paleo_precip.hh"
#include "PATemperaturePIK.hh"
#include "PAConstantPIK.hh"
#include "PAAnomaly.hh"
#include "PACosineYearlyCycle.hh"
#include "PAWeatherStation.hh"

namespace pism {
namespace atmosphere {

// Atmosphere

Factory::Factory(IceGrid::ConstPtr g)
  : PCFactory<AtmosphereModel,PAModifier>(g) {
  m_option = "atmosphere";

  add_model<PIK>("pik");
  add_model<TemperaturePIK>("pik_temp");
  add_model<Given>("given");
  add_model<SeaRISEGreenland>("searise_greenland");
  add_model<CosineYearlyCycle>("yearly_cycle");
  add_model<WeatherStation>("one_station");
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

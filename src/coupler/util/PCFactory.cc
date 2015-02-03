// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 Constantine Khroulev and Torsten Albrecht
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

// Implementation of "surface", "atmosphere" and "ocean" model factories: C++
// classes processing -surface, -atmosphere and -ocean command-line options,
// creating corresponding models and stringing them together to get requested
// data-flow.

#include "PCFactory.hh"
#include "PAFactory.hh"
#include "POFactory.hh"
#include "PSFactory.hh"

#include "pism_const.hh"

// atmosphere models:
#include "PAGivenClimate.hh"
#include "PALapseRates.hh"
#include "PASeariseGreenland.hh"
#include "PA_delta_T.hh"
#include "PA_delta_P.hh"
#include "PA_frac_P.hh"
#include "PA_paleo_precip.hh"
#include "PAConstantPIK.hh"
#include "PAAnomaly.hh"
#include "PACosineYearlyCycle.hh"
#include "PAWeatherStation.hh"

// ocean models:
#include "POConstant.hh"
#include "POConstantPIK.hh"
#include "POGivenClimate.hh"
#include "PO_delta_SL.hh"
#include "PO_delta_T.hh"
#include "PO_delta_SMB.hh"
#include "PO_delta_MBP.hh"
#include "POCache.hh"
#include "POGivenTH.hh" 

// surface models:
#include "PSAnomaly.hh"
#include "PSElevation.hh"
#include "PSGivenClimate.hh"
#include "PSLapseRates.hh"
#include "PSStuffAsAnomaly.hh"
#include "PS_delta_T.hh"
#include "PSTemperatureIndex.hh"
#include "PSTemperatureIndex_Old.hh"
#include "PSSimple.hh"
#include "PSConstantPIK.hh"
#include "PSForceThickness.hh"
#include "PSCache.hh"

namespace pism {

// Atmosphere

PAFactory::PAFactory(const IceGrid& g)
  : PCFactory<AtmosphereModel,PAModifier>(g) {
  m_option = "atmosphere";

  add_model<PAPIK>("pik");
  add_model<PAGivenClimate>("given");
  add_model<PA_SeaRISE_Greenland>("searise_greenland");
  add_model<PACosineYearlyCycle>("yearly_cycle");
  add_model<PAWeatherStation>("one_station");
  set_default("given");

  add_modifier<PAAnomaly>("anomaly");
  add_modifier<PA_paleo_precip>("paleo_precip");
  add_modifier<PA_frac_P>("frac_P");
  add_modifier<PA_delta_P>("delta_P");
  add_modifier<PA_delta_T>("delta_T");
  add_modifier<PALapseRates>("lapse_rate");
}

PAFactory::~PAFactory() {
  // empty
}

// Ocean
POFactory::POFactory(const IceGrid& g)
  : PCFactory<ocean::OceanModel,ocean::OceanModifier>(g) {
  m_option = "ocean";

  add_model<ocean::GivenTH>("th");
  add_model<ocean::PIK>("pik");
  add_model<ocean::Constant>("constant");
  add_model<ocean::Given>("given");
  set_default("constant");

  add_modifier<ocean::Cache>("cache");
  add_modifier<ocean::Delta_SMB>("delta_SMB");
  add_modifier<ocean::Delta_T>("delta_T");
  add_modifier<ocean::Delta_MBP>("delta_MBP");
  add_modifier<ocean::Delta_SL>("delta_SL");
}

POFactory::~POFactory() {
  // empty
}

// Surface
PSFactory::PSFactory(const IceGrid& g)
  : PCFactory<SurfaceModel,PSModifier>(g) {
  m_option = "surface";

  add_model<PSElevation>("elevation");
  add_model<PSGivenClimate>("given");
  add_model<PSTemperatureIndex>("pdd");
  add_model<PSTemperatureIndex_Old>("pdd_old");
  add_model<PSPIK>("pik");
  add_model<PSSimple>("simple");
  set_default("given");

  add_modifier<PSAnomaly>("anomaly");
  add_modifier<PSCache>("cache");
  add_modifier<PS_delta_T>("delta_T");
  add_modifier<PSForceThickness>("forcing");
  add_modifier<PSLapseRates>("lapse_rate");
  add_modifier<PSStuffAsAnomaly>("turn_into_anomaly");
}

PSFactory::~PSFactory() {
  // empty
}

} // end of namespace pism

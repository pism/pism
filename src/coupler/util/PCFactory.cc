// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev and Torsten Albrecht
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
#include "PSSimple.hh"
#include "PSConstantPIK.hh"
#include "PSForceThickness.hh"
#include "PSCache.hh"

// Atmosphere
static void create_pa_constant_pik(IceGrid& g, const PISMConfig& conf, PISMAtmosphereModel* &result) {
  result = new PAConstantPIK(g, conf);
}

static void create_pa_given(IceGrid& g, const PISMConfig& conf, PISMAtmosphereModel* &result) {
  result = new PAGivenClimate(g, conf);
}

static void create_pa_searise_greenland(IceGrid& g, const PISMConfig& conf, PISMAtmosphereModel* &result) {
  result = new PA_SeaRISE_Greenland(g, conf);
}

static void create_pa_yearly_cycle(IceGrid& g, const PISMConfig& conf,
                                   PISMAtmosphereModel* &result) {
  result = new PACosineYearlyCycle(g, conf);
}

static void create_pa_one_station(IceGrid& g, const PISMConfig& conf,
                                  PISMAtmosphereModel* &result) {
  result = new PAWeatherStation(g, conf);
}

static void create_pa_lapse_rates(IceGrid& g, const PISMConfig& conf,
                                  PISMAtmosphereModel *input, PAModifier* &result) {
  result = new PALapseRates(g, conf, input);
}

static void create_pa_delta_T(IceGrid& g, const PISMConfig& conf,
                                PISMAtmosphereModel *input, PAModifier* &result) {
  result = new PA_delta_T(g, conf, input);
}

static void create_pa_delta_P(IceGrid& g, const PISMConfig& conf,
                              PISMAtmosphereModel *input, PAModifier* &result) {
  result = new PA_delta_P(g, conf, input);
}

static void create_pa_paleo_precip(IceGrid& g, const PISMConfig& conf,
                                   PISMAtmosphereModel *input, PAModifier* &result) {
  result = new PA_paleo_precip(g, conf, input);
}

static void create_pa_anomaly(IceGrid& g, const PISMConfig& conf,
                              PISMAtmosphereModel *input, PAModifier* &result) {
  result = new PAAnomaly(g, conf, input);
}

void PAFactory::add_standard_types() {
  add_model("given",             &create_pa_given);
  add_model("searise_greenland", &create_pa_searise_greenland);
  add_model("pik",               &create_pa_constant_pik);
  add_model("yearly_cycle",      &create_pa_yearly_cycle);
  add_model("one_station",       &create_pa_one_station);
  set_default("given");

  add_modifier("anomaly",      &create_pa_anomaly);
  add_modifier("delta_T",      &create_pa_delta_T);
  add_modifier("delta_P",      &create_pa_delta_P);
  add_modifier("paleo_precip", &create_pa_paleo_precip);
  add_modifier("lapse_rate",   &create_pa_lapse_rates);
}


// Ocean
static void create_po_given(IceGrid& g, const PISMConfig& conf, PISMOceanModel* &result) {
  result = new POGiven(g, conf);
}

static void create_po_constant(IceGrid& g, const PISMConfig& conf, PISMOceanModel* &result) {
  result = new POConstant(g, conf);
}

static void create_po_pik(IceGrid& g, const PISMConfig& conf, PISMOceanModel* &result) {
  result = new POConstantPIK(g, conf);
}

static void create_po_th(IceGrid& g, const PISMConfig& conf, PISMOceanModel* &result) {
  result = new POGivenTH(g, conf);
}

static void create_po_delta_SL(IceGrid& g, const PISMConfig& conf, PISMOceanModel *input, POModifier* &result) {
  result = new PO_delta_SL(g, conf, input);
}

static void create_po_delta_MBP(IceGrid& g, const PISMConfig& conf, PISMOceanModel *input, POModifier* &result) {
  result = new PO_delta_MBP(g, conf, input);
}

static void create_po_delta_T(IceGrid& g, const PISMConfig& conf, PISMOceanModel *input, POModifier* &result) {
  result = new PO_delta_T(g, conf, input);
}

static void create_po_delta_SMB(IceGrid& g, const PISMConfig& conf, PISMOceanModel *input, POModifier* &result) {
  result = new PO_delta_SMB(g, conf, input);
}

static void create_po_cache(IceGrid& g, const PISMConfig& conf, PISMOceanModel *input, POModifier* &result) {
  result = new POCache(g, conf, input);
}

void POFactory::add_standard_types() {
  add_model("constant", &create_po_constant);
  add_model("given",    &create_po_given);
  add_model("pik",      &create_po_pik);
  add_model("th",       &create_po_th);
  set_default("constant");

  add_modifier("delta_SL",  &create_po_delta_SL);
  add_modifier("delta_T",   &create_po_delta_T);
  add_modifier("delta_SMB", &create_po_delta_SMB);
  add_modifier("delta_MBP", &create_po_delta_MBP);
  add_modifier("cache",     &create_po_cache);
}

// Surface
static void create_ps_temperatureindex(IceGrid& g, const PISMConfig& conf, PISMSurfaceModel* &result) {
  result = new PSTemperatureIndex(g, conf);
}

static void create_ps_simple(IceGrid& g, const PISMConfig& conf, PISMSurfaceModel* &result) {
  result = new PSSimple(g, conf);
}

static void create_ps_constant_pik(IceGrid& g, const PISMConfig& conf, PISMSurfaceModel* &result) {
  result = new PSConstantPIK(g, conf);
}

static void create_ps_elevation(IceGrid& g, const PISMConfig& conf, PISMSurfaceModel* &result) {
  result = new PSElevation(g, conf);
}

static void create_ps_forcing(IceGrid& g, const PISMConfig& conf,
                              PISMSurfaceModel *input, PSModifier* &result) {
  result = new PSForceThickness(g, conf, input);
}

static void create_ps_lapse_rates(IceGrid& g, const PISMConfig& conf,
                                  PISMSurfaceModel *input, PSModifier* &result) {
  result = new PSLapseRates(g, conf, input);
}

static void create_ps_given(IceGrid& g, const PISMConfig& conf, PISMSurfaceModel* &result) {
  result = new PSGivenClimate(g, conf);
}

static void create_ps_delta_T(IceGrid& g, const PISMConfig& conf,
                              PISMSurfaceModel *input, PSModifier* &result) {
  result = new PS_delta_T(g, conf, input);
}

static void create_ps_cache(IceGrid& g, const PISMConfig& conf,
                            PISMSurfaceModel *input, PSModifier* &result) {
  result = new PSCache(g, conf, input);
}

static void create_ps_stuff_as_anomaly(IceGrid& g, const PISMConfig& conf,
                                       PISMSurfaceModel *input, PSModifier* &result) {
  result = new PSStuffAsAnomaly(g, conf, input);
}

static void create_ps_anomaly(IceGrid& g, const PISMConfig& conf,
                              PISMSurfaceModel *input, PSModifier* &result) {
  result = new PSAnomaly(g, conf, input);
}

void PSFactory::add_standard_types() {
  add_model("simple",    &create_ps_simple);           
  add_model("pdd",       &create_ps_temperatureindex); 
  add_model("given",     &create_ps_given);            
  add_model("pik",       &create_ps_constant_pik);     
  add_model("elevation", &create_ps_elevation);        
  set_default("given");

  add_modifier("anomaly",           &create_ps_anomaly);          
  add_modifier("cache",             &create_ps_cache);            
  add_modifier("forcing",           &create_ps_forcing);          
  add_modifier("delta_T",           &create_ps_delta_T);          
  add_modifier("lapse_rate",        &create_ps_lapse_rates);      
  add_modifier("turn_into_anomaly", &create_ps_stuff_as_anomaly); 
}

// Copyright (C) 2010, 2011 Constantine Khroulev and Torsten Albrecht
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
#include "PSExternal.hh"
#include "PASDirectForcing.hh"
#include "PSElevation.hh"
#include "pism_const.hh"
#include "PASLapseRates.hh"
#include "PASDirectForcing.hh"
#include "PScalarForcing.hh"

// Atmosphere
static void create_pa_constant(IceGrid& g, const NCConfigVariable& conf, PISMAtmosphereModel* &result) {
  result = new PAConstant(g, conf);
}

static void create_pa_given(IceGrid& g, const NCConfigVariable& conf, PISMAtmosphereModel* &result) {
  result = new PADirectForcing(g, conf);
}

static void create_pa_searise_greenland(IceGrid& g, const NCConfigVariable& conf, PISMAtmosphereModel* &result) {
  result = new PA_SeaRISE_Greenland(g, conf);
}

static void create_pa_lapse_rates(IceGrid& g, const NCConfigVariable& conf,
                                  PISMAtmosphereModel *input, PAModifier* &result) {
  result = new PALapseRates(g, conf, input);
}

static void create_pa_anomalies(IceGrid& g, const NCConfigVariable& conf,
                                PISMAtmosphereModel *input, PAModifier* &result) {
  result = new PAForcing(g, conf, input);
}

static void create_pa_dTforcing(IceGrid& g, const NCConfigVariable& conf,
                                PISMAtmosphereModel *input, PAModifier* &result) {
  result = new PAdTforcing(g, conf, input);
}

void PAFactory::add_standard_types() {
  add_model("constant",          &create_pa_constant);
  add_model("given",             &create_pa_given);
  add_model("searise_greenland", &create_pa_searise_greenland);
  set_default("constant");

  add_modifier("anomalies",    &create_pa_anomalies);
  add_modifier("dTforcing",    &create_pa_dTforcing);
  add_modifier("lapse_rate", &create_pa_lapse_rates);
}


// Ocean
static void create_po_constant(IceGrid& g, const NCConfigVariable& conf, PISMOceanModel* &result) {
  result = new POConstant(g, conf);
}

static void create_po_pik(IceGrid& g, const NCConfigVariable& conf, PISMOceanModel* &result) {
  result = new POConstantPIK(g, conf);
}

static void create_po_forcing(IceGrid& g, const NCConfigVariable& conf, PISMOceanModel *input, POModifier* &result) {
  result = new POdSLforcing(g, conf, input);
}

void POFactory::add_standard_types() {
  add_model("constant", &create_po_constant);
  add_model("pik",      &create_po_pik);
  set_default("constant");

  add_modifier("dSLforcing", &create_po_forcing);
}

// Surface
static void create_ps_temperatureindex(IceGrid& g, const NCConfigVariable& conf, PISMSurfaceModel* &result) {
  result = new PSTemperatureIndex(g, conf);
}

static void create_ps_simple(IceGrid& g, const NCConfigVariable& conf, PISMSurfaceModel* &result) {
  result = new PSSimple(g, conf);
}

static void create_ps_constant(IceGrid& g, const NCConfigVariable& conf, PISMSurfaceModel* &result) {
  result = new PSConstant(g, conf);
}

static void create_ps_constant_pik(IceGrid& g, const NCConfigVariable& conf, PISMSurfaceModel* &result) {
  result = new PSConstantPIK(g, conf);
}

static void create_ps_elevation(IceGrid& g, const NCConfigVariable& conf, PISMSurfaceModel* &result) {
  result = new PSElevation(g, conf);
}

static void create_ps_forcing(IceGrid& g, const NCConfigVariable& conf,
                              PISMSurfaceModel *input, PSModifier* &result) {
  result = new PSForceThickness(g, conf, input);
}

static void create_ps_lapse_rates(IceGrid& g, const NCConfigVariable& conf,
                                  PISMSurfaceModel *input, PSModifier* &result) {
  result = new PSLapseRates(g, conf, input);
}

static void create_ps_given(IceGrid& g, const NCConfigVariable& conf, PISMSurfaceModel* &result) {
  result = new PSDirectForcing(g, conf);
}

static void create_ps_dTforcing(IceGrid& g, const NCConfigVariable& conf,
                                PISMSurfaceModel *input, PSModifier* &result) {
  result = new PSdTforcing(g, conf, input);
}

void PSFactory::add_standard_types() {
  add_model("constant",      &create_ps_constant);
  add_model("simple",        &create_ps_simple);
  add_model("pdd",           &create_ps_temperatureindex); 
  add_model("given",         &create_ps_given); 
  add_model("pik",           &create_ps_constant_pik);
  add_model("elevation",     &create_ps_elevation);
  set_default("simple");

  add_modifier("forcing",    &create_ps_forcing);
  add_modifier("dTforcing",  &create_ps_dTforcing);
  add_modifier("lapse_rate", &create_ps_lapse_rates);
}

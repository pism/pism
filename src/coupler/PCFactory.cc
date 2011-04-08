#include "PCFactory.hh"
#include "PSExternal.hh"
#include "PSDirectForcing.hh"
#include "PADirectForcing.hh"
#include "PSElevation.hh"
#include "pism_const.hh"

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

static void create_pa_lapse_rates(IceGrid& g, const NCConfigVariable& conf, PISMAtmosphereModel* &result) {
  result = new PALapseRates(g, conf);
}

static void create_pa_forcing(IceGrid& g, const NCConfigVariable& conf, PAModifier* &result) {
  result = new PAForcing(g, conf);
}

void PAFactory::add_standard_types() {
  add_model("constant",          &create_pa_constant);
  add_model("given",          &create_pa_given);
  add_model("searise_greenland", &create_pa_searise_greenland);
  add_model("temp_lapse_rate",   &create_pa_lapse_rates);
  set_default("constant");

  add_modifier("forcing",        &create_pa_forcing);
}


// Ocean
static void create_po_constant(IceGrid& g, const NCConfigVariable& conf, PISMOceanModel* &result) {
  result = new POConstant(g, conf);
}

static void create_po_pik(IceGrid& g, const NCConfigVariable& conf, PISMOceanModel* &result) {
  result = new POConstantPIK(g, conf);
}

static void create_po_forcing(IceGrid& g, const NCConfigVariable& conf, POModifier* &result) {
  result = new POForcing(g, conf);
}

void POFactory::add_standard_types() {
  add_model("constant",     &create_po_constant);
  add_model("pik",     &create_po_pik);
  set_default("constant");

  add_modifier("forcing",   &create_po_forcing);
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

static void create_ps_forcing(IceGrid& g, const NCConfigVariable& conf, PSModifier* &result) {
  result = new PSForceThickness(g, conf);
}

static void create_ps_given(IceGrid& g, const NCConfigVariable& conf, PISMSurfaceModel* &result) {
  result = new PSDirectForcing(g, conf);
}

void PSFactory::add_standard_types() {
  add_model("constant",     &create_ps_constant);
  add_model("simple",       &create_ps_simple);
  add_model("pdd",          &create_ps_temperatureindex); 
  add_model("given",        &create_ps_given); 
  add_model("pik", &create_ps_constant_pik);
  add_model("elevation", &create_ps_elevation);
  set_default("simple");

  add_modifier("forcing",   &create_ps_forcing);
}

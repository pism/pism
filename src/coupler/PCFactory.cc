#include "PCFactory.hh"
#include "../base/pism_const.hh"

// Atmosphere
static void create_pa_greenland(IceGrid& g, const NCConfigVariable& conf, PISMAtmosphereModel* &result) {
  result = new PAFausto(g, conf);
}

static void create_pa_constant(IceGrid& g, const NCConfigVariable& conf, PISMAtmosphereModel* &result) {
  result = new PAConstant(g, conf);
}

static void create_pa_lapse_rates(IceGrid& g, const NCConfigVariable& conf, PISMAtmosphereModel* &result) {
  result = new PALapseRates(g, conf);
}

static void create_pa_forcing(IceGrid& g, const NCConfigVariable& conf, PAModifier* &result) {
  result = new PAForcing(g, conf);
}

void PAFactory::add_standard_types() {
  add_model("greenland",       &create_pa_greenland);
  add_model("constant",        &create_pa_constant);
  add_model("temp_lapse_rate", &create_pa_lapse_rates);
  add_modifier("forcing",      &create_pa_forcing);
  set_default("constant");
}


// Ocean
static void create_po_constant(IceGrid& g, const NCConfigVariable& conf, PISMOceanModel* &result) {
  result = new POConstant(g, conf);
}

static void create_po_forcing(IceGrid& g, const NCConfigVariable& conf, POModifier* &result) {
  result = new POForcing(g, conf);
}

void POFactory::add_standard_types() {
  add_model("constant",     &create_po_constant);
  add_modifier("forcing",   &create_po_forcing);
  set_default("constant");
}

// Surface
static void create_ps_localmassbalance(IceGrid& g, const NCConfigVariable& conf, PISMSurfaceModel* &result) {
  result = new PSLocalMassBalance(g, conf);
}

static void create_ps_simple(IceGrid& g, const NCConfigVariable& conf, PISMSurfaceModel* &result) {
  result = new PSSimple(g, conf);
}

static void create_ps_constant(IceGrid& g, const NCConfigVariable& conf, PISMSurfaceModel* &result) {
  result = new PSConstant(g, conf);
}

static void create_ps_forcing(IceGrid& g, const NCConfigVariable& conf, PSModifier* &result) {
  result = new PSForceThickness(g, conf);
}

void PSFactory::add_standard_types() {
  add_model("constant",     &create_ps_constant);
  add_model("simple",       &create_ps_simple);
  add_model("pdd",          &create_ps_localmassbalance); 
  add_modifier("forcing",   &create_ps_forcing);
  set_default("simple");
}

// Copyright (C) 2004-2015 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <petscsys.h>

#include "iceModel.hh"

#include "base/basalstrength/PISMYieldStress.hh"
#include "base/basalstrength/basal_resistance.hh"
#include "base/calving/PISMCalvingAtThickness.hh"
#include "base/calving/PISMEigenCalving.hh"
#include "base/calving/PISMFloatKill.hh"
#include "base/calving/PISMIcebergRemover.hh"
#include "base/calving/PISMOceanKill.hh"
#include "base/energy/bedrockThermalUnit.hh"
#include "base/hydrology/PISMHydrology.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "earth/PISMBedDef.hh"
#include "enthalpyConverter.hh"
#include "pism_signal.h"
#include "base/util/PISMVars.hh"
#include "base/util/Profiling.hh"

namespace pism {

IceModel::IceModel(IceGrid::Ptr g, Context::Ptr context)
  : m_grid(g),
    m_config(context->config()),
    m_ctx(context),
    m_sys(context->unit_system()),
    m_log(context->log()),
    m_time(context->time()),
    global_attributes("PISM_GLOBAL", m_sys),
    mapping("mapping", m_sys),
    run_stats("run_stats", m_sys),
    extra_bounds("time_bounds", m_config->get_string("time_dimension_name"), m_sys),
    timestamp("timestamp", m_config->get_string("time_dimension_name"), m_sys) {

  extra_bounds.set_string("units", m_time->units_string());

  timestamp.set_string("units", "hours");
  timestamp.set_string("long_name", "wall-clock time since the beginning of the run");

  pism_signal = 0;
  signal(SIGTERM, pism_signal_handler);
  signal(SIGUSR1, pism_signal_handler);
  signal(SIGUSR2, pism_signal_handler);

  subglacial_hydrology = NULL;
  basal_yield_stress_model = NULL;

  stress_balance = NULL;

  external_surface_model = false;
  external_ocean_model   = false;

  surface = NULL;
  ocean   = NULL;
  beddef  = NULL;

  btu = NULL;

  iceberg_remover             = NULL;
  ocean_kill_calving          = NULL;
  float_kill_calving          = NULL;
  thickness_threshold_calving = NULL;
  eigen_calving               = NULL;

  // initializr maximum |u|,|v|,|w| in ice
  gmaxu = 0;
  gmaxv = 0;
  gmaxw = 0;

  // set default locations of the column used by -view_system
  id = (m_grid->Mx() - 1)/2;
  jd = (m_grid->My() - 1)/2;

  global_attributes.set_string("Conventions", "CF-1.5");
  global_attributes.set_string("source", std::string("PISM ") + PISM_Revision);

  // Do not save snapshots by default:
  save_snapshots = false;
  // Do not save time-series by default:
  save_ts        = false;
  save_extra     = false;

  reset_counters();
}

void IceModel::reset_counters() {
  CFLmaxdt     = 0.0;
  CFLmaxdt2D   = 0.0;
  CFLviolcount = 0;
  dt_TempAge   = 0.0;
  dt_from_cfl  = 0.0;

  gmaxu = 0.0;
  gmaxv = 0.0;
  gmaxw = 0.0;

  maxdt_temporary = 0.0;
  dt              = 0.0;
  dt_force        = 0.0;
  skipCountDown   = 0;

  timestep_hit_multiples_last_time = m_time->current();

  grounded_basal_ice_flux_cumulative = 0;
  nonneg_rule_flux_cumulative        = 0;
  sub_shelf_ice_flux_cumulative      = 0;
  surface_ice_flux_cumulative        = 0;
  sum_divQ_SIA_cumulative            = 0;
  sum_divQ_SSA_cumulative            = 0;
  Href_to_H_flux_cumulative          = 0;
  H_to_Href_flux_cumulative          = 0;
  discharge_flux_cumulative          = 0;
}


IceModel::~IceModel() {

  // de-allocate time-series diagnostics
  std::map<std::string,TSDiagnostic*>::iterator i = ts_diagnostics.begin();
  while (i != ts_diagnostics.end()) {
    delete (i++)->second;
  }

  // de-allocate diagnostics
  std::map<std::string,Diagnostic*>::iterator j = diagnostics.begin();
  while (j != diagnostics.end()) {
    delete (j++)->second;
  }

  delete stress_balance;

  if (external_ocean_model == false) {
    delete ocean;
  }

  if (external_surface_model == false) {
    delete surface;
  }

  delete beddef;

  delete subglacial_hydrology;
  delete basal_yield_stress_model;
  delete btu;

  delete iceberg_remover;
  delete ocean_kill_calving;
  delete float_kill_calving;
  delete thickness_threshold_calving;
  delete eigen_calving;
}


//! Allocate all IceModelVecs defined in IceModel.
/*!
  This procedure allocates the memory used to store model state, diagnostic and
  work vectors and sets metadata.

  Default values should not be set here; please use set_vars_from_options().

  All the memory allocated here is freed by IceModelVecs' destructors.
*/
void IceModel::createVecs() {

  const unsigned int WIDE_STENCIL = m_config->get_double("grid_max_stencil_width");

  m_log->message(3,
                 "Allocating memory...\n");

  // get the list of selected calving methods:
  std::istringstream calving_methods_list(m_config->get_string("calving_methods"));
  std::string calving_method_name;
  std::set<std::string> calving_methods;

  while (getline(calving_methods_list, calving_method_name, ',')) {
    calving_methods.insert(calving_method_name);
  }

  // The following code creates (and documents -- to some extent) the
  // variables. The main (and only) principle here is using standard names from
  // the CF conventions; see
  // http://cf-pcmdi.llnl.gov/documents/cf-standard-names

  Enth3.create(m_grid, "enthalpy", WITH_GHOSTS, WIDE_STENCIL);
  // POSSIBLE standard name = land_ice_enthalpy
  Enth3.set_attrs("model_state",
                  "ice enthalpy (includes sensible heat, latent heat, pressure)",
                  "J kg-1", "");
  m_grid->variables().add(Enth3);

  if (m_config->get_boolean("do_cold_ice_methods")) {
    // ice temperature
    T3.create(m_grid, "temp", WITH_GHOSTS);
    T3.set_attrs("model_state",
                 "ice temperature", "K", "land_ice_temperature");
    T3.metadata().set_double("valid_min", 0.0);
    m_grid->variables().add(T3);

    if (m_config->get_boolean("do_energy") == true) {
      Enth3.metadata().set_string("pism_intent", "diagnostic");
    } else {
      T3.metadata().set_string("pism_intent", "diagnostic");
    }
  }

  // age of ice but only if age will be computed
  if (m_config->get_boolean("do_age")) {
    age3.create(m_grid, "age", WITH_GHOSTS, WIDE_STENCIL);
    // PROPOSED standard_name = land_ice_age
    age3.set_attrs("model_state", "age of ice",
                   "s", "");
    age3.metadata().set_string("glaciological_units", "years");
    age3.write_in_glaciological_units = true;
    age3.metadata().set_double("valid_min", 0.0);
    m_grid->variables().add(age3);
  }
  if (m_config->get_boolean("do_tracer")) {
    tracer_x.create(m_grid, "tracer_x", WITH_GHOSTS, WIDE_STENCIL);
    tracer_x.set_attrs("model_state", "x-origin of ice",
                   "m", "");
    m_grid->variables().add(tracer_x);

    tracer_y.create(m_grid, "tracer_y", WITH_GHOSTS, WIDE_STENCIL);
    tracer_y.set_attrs("model_state", "y-origin of ice",
                   "m", "");
    m_grid->variables().add(tracer_y);

    tracer_z.create(m_grid, "tracer_z", WITH_GHOSTS, WIDE_STENCIL);
    tracer_z.set_attrs("model_state", "z-origin of ice",
                   "m", "");
    m_grid->variables().add(tracer_z);

    tracer_t.create(m_grid, "tracer_t", WITH_GHOSTS, WIDE_STENCIL);
    tracer_t.set_attrs("model_state", "time of creation of ice",
                   "s", "");
    tracer_t.metadata().set_string("glaciological_units", "years");
    tracer_t.write_in_glaciological_units = true;
    m_grid->variables().add(tracer_t);
  }

  // ice upper surface elevation
  ice_surface_elevation.create(m_grid, "usurf", WITH_GHOSTS, WIDE_STENCIL);
  ice_surface_elevation.set_attrs("diagnostic", "ice upper surface elevation",
                                  "m", "surface_altitude");
  m_grid->variables().add(ice_surface_elevation);

  // land ice thickness
  ice_thickness.create(m_grid, "thk", WITH_GHOSTS, WIDE_STENCIL);
  ice_thickness.set_attrs("model_state", "land ice thickness",
                          "m", "land_ice_thickness");
  ice_thickness.metadata().set_double("valid_min", 0.0);
  m_grid->variables().add(ice_thickness);

  if (m_config->get_boolean("sub_groundingline")) {
    gl_mask.create(m_grid, "gl_mask", WITHOUT_GHOSTS);
    gl_mask.set_attrs("internal",
                      "fractional grounded/floating mask (floating=0, grounded=1)",
                      "", "");
    m_grid->variables().add(gl_mask);

    gl_mask_x.create(m_grid, "gl_mask_x", WITHOUT_GHOSTS);
    gl_mask_x.set_attrs("internal",
                        "fractional grounded/floating mask in x-direction (floating=0, grounded=1)",
                        "", "");
    m_grid->variables().add(gl_mask_x);

    gl_mask_y.create(m_grid, "gl_mask_y", WITHOUT_GHOSTS);
    gl_mask_y.set_attrs("internal",
                        "fractional grounded/floating mask in y-direction (floating=0, grounded=1)",
                        "", "");
    m_grid->variables().add(gl_mask_y);
  }

  // grounded_dragging_floating integer mask
  vMask.create(m_grid, "mask", WITH_GHOSTS, WIDE_STENCIL);
  vMask.set_attrs("diagnostic", "ice-type (ice-free/grounded/floating/ocean) integer mask",
                  "", "");
  std::vector<double> mask_values(4);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_GROUNDED;
  mask_values[2] = MASK_FLOATING;
  mask_values[3] = MASK_ICE_FREE_OCEAN;
  vMask.metadata().set_doubles("flag_values", mask_values);
  vMask.metadata().set_string("flag_meanings",
                              "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
  m_grid->variables().add(vMask);

  // upward geothermal flux at bedrock surface
  geothermal_flux.create(m_grid, "bheatflx", WITHOUT_GHOSTS);
  // PROPOSED standard_name = lithosphere_upward_heat_flux
  geothermal_flux.set_attrs("climate_steady", "upward geothermal flux at bedrock surface",
                            "W m-2", "");
  geothermal_flux.metadata().set_string("glaciological_units", "mW m-2");
  geothermal_flux.write_in_glaciological_units = true;
  geothermal_flux.set_time_independent(true);
  m_grid->variables().add(geothermal_flux);

  // temperature seen by top of bedrock thermal layer
  bedtoptemp.create(m_grid, "bedtoptemp", WITHOUT_GHOSTS);
  bedtoptemp.set_attrs("internal",
                       "temperature of top of bedrock thermal layer",
                       "K", "");
  bedtoptemp.metadata().set_string("glaciological_units", "K");
  m_grid->variables().add(bedtoptemp);

  // yield stress for basal till (plastic or pseudo-plastic model)
  {
    basal_yield_stress.create(m_grid, "tauc", WITH_GHOSTS, WIDE_STENCIL);
    // PROPOSED standard_name = land_ice_basal_material_yield_stress
    basal_yield_stress.set_attrs("diagnostic",
                                 "yield stress for basal till (plastic or pseudo-plastic model)",
                                 "Pa", "");
    m_grid->variables().add(basal_yield_stress);
  }

  // basal melt rate
  basal_melt_rate.create(m_grid, "bmelt", WITHOUT_GHOSTS);
  // ghosted to allow the "redundant" computation of tauc
  basal_melt_rate.set_attrs("model_state",
                            "ice basal melt rate from energy conservation and subshelf melt, in ice thickness per time",
                            "m s-1", "land_ice_basal_melt_rate");
  basal_melt_rate.metadata().set_string("glaciological_units", "m year-1");
  basal_melt_rate.write_in_glaciological_units = true;
  basal_melt_rate.metadata().set_string("comment", "positive basal melt rate corresponds to ice loss");
  m_grid->variables().add(basal_melt_rate);

  // longitude
  vLongitude.create(m_grid, "lon", WITH_GHOSTS);
  vLongitude.set_attrs("mapping", "longitude", "degree_east", "longitude");
  vLongitude.set_time_independent(true);
  vLongitude.metadata().set_string("coordinates", "");
  vLongitude.metadata().set_string("grid_mapping", "");
  vLongitude.metadata().set_double("valid_min", -180.0);
  vLongitude.metadata().set_double("valid_max",  180.0);
  m_grid->variables().add(vLongitude);

  // latitude
  vLatitude.create(m_grid, "lat", WITH_GHOSTS); // has ghosts so that we can compute cell areas
  vLatitude.set_attrs("mapping", "latitude", "degree_north", "latitude");
  vLatitude.set_time_independent(true);
  vLatitude.metadata().set_string("coordinates", "");
  vLatitude.metadata().set_string("grid_mapping", "");
  vLatitude.metadata().set_double("valid_min", -90.0);
  vLatitude.metadata().set_double("valid_max",  90.0);
  m_grid->variables().add(vLatitude);

  if (m_config->get_boolean("part_grid") == true) {
    // Href
    vHref.create(m_grid, "Href", WITH_GHOSTS);
    vHref.set_attrs("model_state", "temporary ice thickness at calving front boundary",
                    "m", "");
    m_grid->variables().add(vHref);
  }

  if (m_config->get_string("calving_methods").find("eigen_calving") != std::string::npos ||
      m_config->get_boolean("do_fracture_density") == true) {

    strain_rates.create(m_grid, "edot", WITH_GHOSTS,
                        2, // stencil width, has to match or exceed the "offset" in eigenCalving
                        2);

    strain_rates.metadata(0).set_name("edot_1");
    strain_rates.set_attrs("internal",
                           "major principal component of horizontal strain-rate",
                           "1/s", "", 0);

    strain_rates.metadata(1).set_name("edot_2");
    strain_rates.set_attrs("internal",
                           "minor principal component of horizontal strain-rate",
                           "1/s", "", 1);
  }

  if (m_config->get_boolean("do_fracture_density") == true) {

    deviatoric_stresses.create(m_grid, "sigma", WITH_GHOSTS,
                               2, // stencil width
                               3);

    deviatoric_stresses.metadata(0).set_name("sigma_xx");
    deviatoric_stresses.set_attrs("internal",
                                  "deviatoric stress in x direction",
                                  "Pa", "", 0);

    deviatoric_stresses.metadata(1).set_name("sigma_yy");
    deviatoric_stresses.set_attrs("internal",
                                  "deviatoric stress in y direction",
                                  "Pa", "", 1);

    deviatoric_stresses.metadata(2).set_name("sigma_xy");
    deviatoric_stresses.set_attrs("internal",
                                  "deviatoric shear stress",
                                  "Pa", "", 2);
  }

  if (m_config->get_boolean("ssa_dirichlet_bc") == true) {
    // bc_locations
    vBCMask.create(m_grid, "bc_mask", WITH_GHOSTS, WIDE_STENCIL);
    vBCMask.set_attrs("model_state", "Dirichlet boundary mask",
                      "", "");
    std::vector<double> bc_mask_values(2);
    bc_mask_values[0] = 0;
    bc_mask_values[1] = 1;
    vBCMask.metadata().set_doubles("flag_values", bc_mask_values);
    vBCMask.metadata().set_string("flag_meanings", "no_data bc_condition");
    m_grid->variables().add(vBCMask);


    Config::ConstPtr config = m_grid->ctx()->config();
    units::System::Ptr sys = m_grid->ctx()->unit_system();
    double fill_value = units::convert(sys, config->get_double("fill_value"), "m/year", "m/s");
    // vel_bc
    vBCvel.create(m_grid, "_ssa_bc", WITH_GHOSTS, WIDE_STENCIL); // u_ssa_bc and v_ssa_bc
    vBCvel.set_attrs("model_state",
                     "X-component of the SSA velocity boundary conditions",
                     "m s-1", "", 0);
    vBCvel.set_attrs("model_state",
                     "Y-component of the SSA velocity boundary conditions",
                     "m s-1", "", 1);
    for (int j = 0; j < 2; ++j) {
      vBCvel.metadata(j).set_string("glaciological_units", "m year-1");
      vBCvel.metadata(j).set_double("valid_min",  units::convert(m_sys, -1e6, "m/year", "m/second"));
      vBCvel.metadata(j).set_double("valid_max",  units::convert(m_sys,  1e6, "m/year", "m/second"));
      vBCvel.metadata(j).set_double("_FillValue", fill_value);
    }
    //just for diagnostics...
    m_grid->variables().add(vBCvel);
  }

  // fracture density field
  if (m_config->get_boolean("do_fracture_density")) {
    vFD.create(m_grid, "fracture_density", WITH_GHOSTS, WIDE_STENCIL);
    vFD.set_attrs("model_state", "fracture density in ice shelf", "", "");
    vFD.metadata().set_double("valid_max", 1.0);
    vFD.metadata().set_double("valid_min", 0.0);
    m_grid->variables().add(vFD);

    if (m_config->get_boolean("write_fd_fields")) {
      vFG.create(m_grid, "fracture_growth_rate", WITH_GHOSTS, WIDE_STENCIL);
      vFG.set_attrs("model_state", "fracture growth rate",       "1/s", "");
      vFG.metadata().set_double("valid_min", 0.0);
      m_grid->variables().add(vFG);

      vFH.create(m_grid, "fracture_healing_rate", WITH_GHOSTS, WIDE_STENCIL);
      vFH.set_attrs("model_state", "fracture healing rate",      "1/s", "");
      m_grid->variables().add(vFH);

      vFE.create(m_grid, "fracture_flow_enhancement", WITH_GHOSTS, WIDE_STENCIL);
      vFE.set_attrs("model_state", "fracture-induced flow enhancement", "", "");
      m_grid->variables().add(vFE);

      vFA.create(m_grid, "fracture_age", WITH_GHOSTS, WIDE_STENCIL);
      vFA.set_attrs("model_state", "age since fracturing",       "years", "");
      m_grid->variables().add(vFA);

      vFT.create(m_grid, "fracture_toughness", WITH_GHOSTS, WIDE_STENCIL);
      vFT.set_attrs("model_state", "fracture toughness", "Pa", "");
      m_grid->variables().add(vFT);
    }
  }

  // cell areas
  cell_area.create(m_grid, "cell_area", WITHOUT_GHOSTS);
  cell_area.set_attrs("diagnostic", "cell areas", "m2", "");
  cell_area.metadata().set_string("comment",
                                  "values are equal to dx*dy "
                                  "if projection parameters are not available; "
                                  "otherwise WGS84 ellipsoid is used");
  cell_area.set_time_independent(true);
  cell_area.metadata().set_string("glaciological_units", "km2");
  cell_area.write_in_glaciological_units = true;
  m_grid->variables().add(cell_area);

  // fields owned by IceModel but filled by SurfaceModel *surface:
  // mean annual net ice equivalent surface mass balance rate
  climatic_mass_balance.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  climatic_mass_balance.set_attrs("climate_from_SurfaceModel",  // FIXME: can we do better?
                                  "ice-equivalent surface mass balance (accumulation/ablation) rate",
                                  "kg m-2 s-1",
                                  "land_ice_surface_specific_mass_balance_flux");
  climatic_mass_balance.metadata().set_string("glaciological_units", "kg m-2 year-1");
  climatic_mass_balance.write_in_glaciological_units = true;
  climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // annual mean air temperature at "ice surface", at level below all firn
  //   processes (e.g. "10 m" or ice temperatures)
  ice_surface_temp.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
  ice_surface_temp.set_attrs("climate_from_SurfaceModel",  // FIXME: can we do better?
                             "annual average ice surface temperature, below firn processes",
                             "K",
                             "");

  liqfrac_surface.create(m_grid, "liqfrac_surface", WITHOUT_GHOSTS);
  liqfrac_surface.set_attrs("climate_from_SurfaceModel",
                            "liquid water fraction at the top surface of the ice",
                            "1", "");
  // grid.variables().add(liqfrac_surface);

  // ice mass balance rate at the base of the ice shelf; sign convention for
  //   vshelfbasemass matches standard sign convention for basal melt rate of
  //   grounded ice
  shelfbmassflux.create(m_grid, "shelfbmassflux", WITHOUT_GHOSTS); // no ghosts; NO HOR. DIFF.!
  shelfbmassflux.set_attrs("climate_state", "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                           "m s-1", "");
  // PROPOSED standard name = ice_shelf_basal_specific_mass_balance
  // rescales from m/s to m/year when writing to NetCDF and std out:
  shelfbmassflux.write_in_glaciological_units = true;
  shelfbmassflux.metadata().set_string("glaciological_units", "m year-1");
  // do not add; boundary models are in charge here
  // grid.variables().add(shelfbmassflux);

  // ice boundary tempature at the base of the ice shelf
  shelfbtemp.create(m_grid, "shelfbtemp", WITHOUT_GHOSTS); // no ghosts; NO HOR. DIFF.!
  shelfbtemp.set_attrs("climate_state", "absolute temperature at ice shelf base",
                       "K", "");
  // PROPOSED standard name = ice_shelf_basal_temperature
  // do not add; boundary models are in charge here
  // grid.variables().add(shelfbtemp);

  // take care of 2D cumulative fluxes: we need to allocate special storage if
  // the user asked for climatic_mass_balance_cumulative or some others (below).

  options::StringSet extras("-extra_vars",
                            "list of spatially-variable diagnostics to save",
                            ""); // don't save anything by default

  if (set_contains(extras, "climatic_mass_balance_cumulative")) {
    climatic_mass_balance_cumulative.create(m_grid,
                                            "climatic_mass_balance_cumulative",
                                            WITHOUT_GHOSTS);
    climatic_mass_balance_cumulative.set_attrs("diagnostic",
                                               "cumulative surface mass balance",
                                               "kg m-2", "");
  }

  std::string o_size = get_output_size("-o_size");

  if (set_contains(extras, "flux_divergence") || o_size == "big") {
    flux_divergence.create(m_grid, "flux_divergence", WITHOUT_GHOSTS);
    flux_divergence.set_attrs("diagnostic",
                              "flux divergence",
                              "m s-1", "");
    flux_divergence.metadata().set_string("glaciological_units", "m year-1");
    flux_divergence.write_in_glaciological_units = true;
  }

  if (set_contains(extras, "grounded_basal_flux_cumulative")) {
    grounded_basal_flux_2D_cumulative.create(m_grid, "grounded_basal_flux_cumulative", WITHOUT_GHOSTS);
    grounded_basal_flux_2D_cumulative.set_attrs("diagnostic",
                                                "cumulative basal flux into the ice "
                                                "in grounded areas (positive means ice gain)",
                                                "kg m-2", "");
    grounded_basal_flux_2D_cumulative.set_time_independent(false);
    grounded_basal_flux_2D_cumulative.metadata().set_string("glaciological_units", "Gt m-2");
    grounded_basal_flux_2D_cumulative.write_in_glaciological_units = true;
  }

  if (set_contains(extras, "floating_basal_flux_cumulative")) {
    floating_basal_flux_2D_cumulative.create(m_grid, "floating_basal_flux_cumulative", WITHOUT_GHOSTS);
    floating_basal_flux_2D_cumulative.set_attrs("diagnostic",
                                                "cumulative basal flux into the ice "
                                                "in floating areas (positive means ice gain)",
                                                "kg m-2", "");
    floating_basal_flux_2D_cumulative.set_time_independent(false);
    floating_basal_flux_2D_cumulative.metadata().set_string("glaciological_units", "Gt m-2");
    floating_basal_flux_2D_cumulative.write_in_glaciological_units = true;
  }

  if (set_contains(extras, "nonneg_flux_cumulative")) {
    nonneg_flux_2D_cumulative.create(m_grid, "nonneg_flux_cumulative", WITHOUT_GHOSTS);
    nonneg_flux_2D_cumulative.set_attrs("diagnostic",
                                        "cumulative nonnegative rule flux (positive means ice gain)",
                                        "kg m-2", "");
    nonneg_flux_2D_cumulative.set_time_independent(false);
    nonneg_flux_2D_cumulative.metadata().set_string("glaciological_units", "Gt m-2");
    nonneg_flux_2D_cumulative.write_in_glaciological_units = true;
  }

  if (set_contains(extras, "discharge_flux_cumulative")) {
    discharge_flux_2D_cumulative.create(m_grid, "discharge_flux_cumulative", WITHOUT_GHOSTS);
    discharge_flux_2D_cumulative.set_attrs("diagnostic",
                                           "cumulative discharge (calving) flux (positive means ice loss)",
                                           "kg m-2", "");
    discharge_flux_2D_cumulative.set_time_independent(false);
    discharge_flux_2D_cumulative.metadata().set_string("glaciological_units", "Gt m-2");
    discharge_flux_2D_cumulative.write_in_glaciological_units = true;
  }
}

//! The contents of the main PISM time-step.
/*!
During the time-step we perform the following actions:
 */
void IceModel::step(bool do_mass_continuity,
                              bool do_energy,
                              bool do_age,
           										bool do_tracer,
                              bool do_skip) {

  const Profiling &profiling = m_ctx->profiling();

  double current_time = m_time->current();

  //! \li call additionalAtStartTimestep() to let derived classes do more
  additionalAtStartTimestep();  // might set dt_force,maxdt_temporary

  //! \li update the velocity field; in some cases the whole three-dimensional
  //! field is updated and in some cases just the vertically-averaged
  //! horizontal velocity is updated

  // always "update" ice velocity (possibly trivially); only update
  // SSA and only update velocities at depth if suggested by temp and age
  // stability criterion; note *lots* of communication is avoided by skipping
  // SSA (and temp/age)

  bool updateAtDepth = (skipCountDown == 0),
    do_energy_step = updateAtDepth && do_energy;

  //! \li update the yield stress for the plastic till model (if appropriate)
  if (updateAtDepth && basal_yield_stress_model) {
    profiling.begin("basal yield stress");
    basal_yield_stress_model->update(current_time, dt);
    profiling.end("basal yield stress");
    basal_yield_stress.copy_from(basal_yield_stress_model->basal_material_yield_stress());
    stdout_flags += "y";
  } else {
    stdout_flags += "$";
  }

  // Update the fractional grounded/floating mask (used by the SSA
  // stress balance and the energy code)
  if (m_config->get_boolean("sub_groundingline")) {
    updateSurfaceElevationAndMask(); // update h and mask
    update_floatation_mask();
  }

  double sea_level = ocean->sea_level_elevation();

  IceModelVec2S &melange_back_pressure = vWork2d[0];

  ocean->melange_back_pressure_fraction(melange_back_pressure);

  try {
    profiling.begin("stress balance");
    stress_balance->update(updateAtDepth == false,
                           sea_level,
                           melange_back_pressure);
    profiling.end("stress balance");
  } catch (RuntimeError &e) {
    options::String output_file("-o", "output file name",
                                "output.nc", options::DONT_ALLOW_EMPTY);

    std::string o_file = pism_filename_add_suffix(output_file,
                                                  "_stressbalance_failed", "");
    dumpToFile(o_file);

    throw RuntimeError::formatted("stress balance computation failed. Model state was saved to '%s'",
                                  o_file.c_str());
  }

  stdout_flags += stress_balance->stdout_report();

  stdout_flags += (updateAtDepth ? "v" : "V");

  //! \li determine the time step according to a variety of stability criteria;
  //!  see determineTimeStep()
  max_timestep(dt, skipCountDown);

  //! \li Update surface and ocean models.
  profiling.begin("surface");
  surface->update(current_time, dt);
  profiling.end("surface");

  profiling.begin("ocean");
  ocean->update(current_time,   dt);
  profiling.end("ocean");

  dt_TempAge += dt;
  // IceModel::dt,dtTempAge are now set correctly according to
  // mass-continuity-eqn-diffusivity criteria, horizontal CFL criteria, and
  // other criteria from derived class additionalAtStartTimestep(), and from
  // "-skip" mechanism

  //! \li update the age of the ice (if appropriate)
  if (do_age && updateAtDepth) {
    profiling.begin("age");
    ageStep();
    profiling.end("age");
    stdout_flags += "a";
  } else {
    stdout_flags += "$";
  }

  //! \li update the tracers in the ice (if appropriate)
  if (do_tracer && updateAtDepth) {
    profiling.begin("tracer");
    tracerStep();
    profiling.end("tracer");
    stdout_flags += "t";
  } else {
    stdout_flags += "$";
  }

  //! \li update the enthalpy (or temperature) field according to the conservation of
  //!  energy model based (especially) on the new velocity field; see
  //!  energyStep()
  if (do_energy_step) { // do the energy step
    profiling.begin("energy");
    energyStep();
    profiling.end("energy");
    stdout_flags += "E";
  } else {
    stdout_flags += "$";
  }

  // Combine basal melt rate in grounded (computed during the energy
  // step) and floating (provided by an ocean model) areas.
  combine_basal_melt_rate();

  //! \li update the state variables in the subglacial hydrology model (typically
  //!  water thickness and sometimes pressure)
  profiling.begin("basal hydrology");
  subglacial_hydrology->update(current_time, dt);
  profiling.end("basal hydrology");

  //! \li update the fracture density field; see calculateFractureDensity()
  if (m_config->get_boolean("do_fracture_density")) {
    profiling.begin("fracture density");
    calculateFractureDensity();
    profiling.end("fracture density");
  }

  //! \li update the thickness of the ice according to the mass conservation
  //!  model; see massContExplicitStep()
  if (do_mass_continuity) {
    profiling.begin("mass transport");
    massContExplicitStep();
    updateSurfaceElevationAndMask(); // update h and mask
    profiling.end("mass transport");

    // Note that there are three adaptive time-stepping criteria. Two of them
    // (using max. diffusion and 2D CFL) are limiting the mass-continuity
    // time-step and the third (3D CFL) limits the energy and age time-steps.

    // The mass-continuity time-step is usually smaller, and the skipping
    // mechanism lets us do several mass-continuity steps for each energy step.

    // When -no_mass is set, mass-continuity-related time-step restrictions are
    // disabled, making "skipping" unnecessary.

    // This is why the following two lines appear here and are executed only
    // if do_mass_continuity == true.
    if (do_skip == true && skipCountDown > 0) {
      skipCountDown--;
    }
    stdout_flags += "h";
  } else {
    stdout_flags += "$";
  }

  profiling.begin("calving");
  do_calving();
  profiling.end("calving");

  Href_cleanup();

  //! \li compute the bed deformation, which only depends on current thickness
  //! and bed elevation
  if (beddef) {
    const IceModelVec2S &bed_topography = beddef->bed_elevation();
    int topg_state_counter = bed_topography.get_state_counter();

    profiling.begin("bed deformation");
    beddef->update(current_time, dt);
    profiling.end("bed deformation");

    if (bed_topography.get_state_counter() != topg_state_counter) {
      stdout_flags += "b";
      updateSurfaceElevationAndMask();
    } else {
      stdout_flags += " ";
    }
  }

  //! \li call additionalAtEndTimestep() to let derived classes do more
  additionalAtEndTimestep();

  // Done with the step; now adopt the new time.
  m_time->step(dt);

  if (do_energy_step) {
    t_TempAge = m_time->current();
    dt_TempAge = 0.0;
  }

  // end the flag line
  stdout_flags += " " + m_adaptive_timestep_reason;
}


/**
 * Run the time-stepping loop from the current model time to `time`.
 *
 * This should be called by the coupler controlling PISM when it is
 * running alongside a GCM.
 *
 * @param run_end model time (in seconds) to run to
 *
 * @return 0 on success
 */
void IceModel::run_to(double run_end) {

  m_time->set_end(run_end);

  run();
}


/**
 * Run the time-stepping loop from the current time until the time
 * specified by the IceModel::grid::time object.
 *
 * This is the method used by PISM in the "standalone" mode.
 *
 * @return 0 on success
 */
void IceModel::run() {
  const Profiling &profiling = m_ctx->profiling();

  bool do_mass_conserve = m_config->get_boolean("do_mass_conserve");
  bool do_energy = m_config->get_boolean("do_energy");
  bool do_age = m_config->get_boolean("do_age");
  bool do_tracer = m_config->get_boolean("do_tracer");
  bool do_skip = m_config->get_boolean("do_skip");

  int stepcount = m_config->get_boolean("count_time_steps") ? 0 : -1;

  updateSurfaceElevationAndMask();

  // update diagnostics at the beginning of the run:
  write_timeseries();
  write_extras();

  m_log->message(2, "running forward ...\n");

  stdout_flags.erase(); // clear it out
  summaryPrintLine(true, do_energy, 0.0, 0.0, 0.0, 0.0, 0.0);
  m_adaptive_timestep_reason = '$'; // no reason for no timestep
  summary(do_energy);  // report starting state

  t_TempAge = m_time->current();
  dt_TempAge = 0.0;

  // main loop for time evolution
  // IceModel::step calls Time::step(dt), ensuring that this while loop
  // will terminate
  profiling.stage_begin("time-stepping loop");
  while (m_time->current() < m_time->end()) {

    stdout_flags.erase();  // clear it out
    dt_force = -1.0;
    maxdt_temporary = -1.0;

    step(do_mass_conserve, do_energy, do_age, do_tracer, do_skip);

    // report a summary for major steps or the last one
    bool updateAtDepth = skipCountDown == 0;
    bool tempAgeStep = updateAtDepth && (do_energy || do_age || do_tracer);

    const bool show_step = tempAgeStep || m_adaptive_timestep_reason == "end of the run";
    summary(show_step);

    // writing these fields here ensures that we do it after the last time-step
    profiling.begin("I/O during run");
    write_snapshot();
    write_timeseries();
    write_extras();
    write_backup();
    profiling.end("I/O during run");

    update_viewers();

    if (stepcount >= 0) {
      stepcount++;
    }
    if (endOfTimeStepHook() != 0) {
      break;
    }
  } // end of the time-stepping loop

  profiling.stage_end("time-stepping loop");

  options::Integer pause_time("-pause", "Pause after the run, seconds", 0);
  if (pause_time > 0) {
    m_log->message(2, "pausing for %d secs ...\n", pause_time.value());
    PetscErrorCode ierr = PetscSleep(pause_time); PISM_CHK(ierr, "PetscSleep");
  }

  if (stepcount >= 0) {
    m_log->message(1,
               "count_time_steps:  run() took %d steps\n"
               "average dt = %.6f years\n",
               stepcount,
               units::convert(m_sys, m_time->end() - m_time->start(), "seconds", "years")/(double)stepcount);
  }
}

//! Manage the initialization of the IceModel object.
/*!
Please see the documenting comments of the functions called below to find
explanations of their intended uses.
 */
void IceModel::init() {
  const Profiling &profiling = m_ctx->profiling();

  profiling.begin("initialization");

  // Build PISM with -DPISM_WAIT_FOR_GDB=1 and run with -wait_for_gdb to
  // make it wait for a connection.
#ifdef PISM_WAIT_FOR_GDB
  bool wait_for_gdb = options::Bool("-wait_for_gdb", "wait for GDB to attach");
  if (wait_for_gdb.is_set()) {
    pism_wait_for_gdb(grid.com, 0);
  }
#endif
  //! The IceModel initialization sequence is this:

  //! 1) Initialize model time:
  time_setup();

  //! 2) Process the options:
  setFromOptions();

  //! 3) Memory allocation:
  createVecs();

  //! 4) Allocate PISM components modeling some physical processes.
  allocate_submodels();

  //! 5) Allocate work vectors:
  allocate_internal_objects();

  //! 6) Initialize coupler models and fill the model state variables
  //! (from a PISM output file, from a bootstrapping file using some
  //! modeling choices or using formulas). Calls IceModel::regrid()
  model_state_setup();

  //! 7) Report grid parameters:
  m_grid->report_parameters();

  //! 8) Miscellaneous stuff: set up the bed deformation model, initialize the
  //! basal till model, initialize snapshots. This has to happen *after*
  //! regridding.
  misc_setup();

  // Get the start time in seconds and ensure that it is consistent
  // across all processors.
  start_time = GlobalMax(m_grid->com, GetTime());

  profiling.end("initialization");
}

// FIXME: THIS IS BAD! (Provides unguarded access to IceModel's internals.)
IceModelVec2S* IceModel::get_geothermal_flux() {
  return &this->geothermal_flux;
}

// FIXME: THIS IS BAD! (Provides unguarded access to IceModel's internals.)
stressbalance::StressBalance* IceModel::get_stress_balance() {
  return this->stress_balance;
}

} // end of namespace pism

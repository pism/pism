// Copyright (C) 2004-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "base/util/pism_utilities.hh"

namespace pism {

IceModel::IceModel(IceGrid::Ptr g, Context::Ptr context)
  : m_grid(g),
    m_config(context->config()),
    m_ctx(context),
    m_sys(context->unit_system()),
    m_log(context->log()),
    m_time(context->time()),
    m_output_global_attributes("PISM_GLOBAL", m_sys),
    mapping("mapping", m_sys),
    run_stats("run_stats", m_sys),
    m_extra_bounds("time_bounds", m_config->get_string("time_dimension_name"), m_sys),
    m_timestamp("timestamp", m_config->get_string("time_dimension_name"), m_sys) {

  m_extra_bounds.set_string("units", m_time->units_string());

  m_timestamp.set_string("units", "hours");
  m_timestamp.set_string("long_name", "wall-clock time since the beginning of the run");

  pism_signal = 0;
  signal(SIGTERM, pism_signal_handler);
  signal(SIGUSR1, pism_signal_handler);
  signal(SIGUSR2, pism_signal_handler);

  subglacial_hydrology = NULL;
  basal_yield_stress_model = NULL;

  m_stress_balance = NULL;

  m_external_surface_model = false;
  m_external_ocean_model   = false;

  m_surface = NULL;
  m_ocean   = NULL;
  m_beddef  = NULL;

  btu = NULL;

  iceberg_remover             = NULL;
  ocean_kill_calving          = NULL;
  float_kill_calving          = NULL;
  thickness_threshold_calving = NULL;
  eigen_calving               = NULL;

  // initialize maximum |u|,|v|,|w| in ice
  m_max_u_speed = 0;
  m_max_v_speed = 0;
  m_max_w_speed = 0;

  // set default locations of the column used by -view_system
  m_id = (m_grid->Mx() - 1)/2;
  m_jd = (m_grid->My() - 1)/2;

  m_output_global_attributes.set_string("Conventions", "CF-1.5");
  m_output_global_attributes.set_string("source", std::string("PISM ") + PISM_Revision);

  // Do not save snapshots by default:
  m_save_snapshots = false;
  // Do not save time-series by default:
  m_save_ts        = false;
  m_save_extra     = false;

  reset_counters();
}

void IceModel::reset_counters() {
  CFLmaxdt     = 0.0;
  CFLmaxdt2D   = 0.0;
  CFLviolcount = 0;
  dt_TempAge   = 0.0;
  dt_from_cfl  = 0.0;

  m_max_u_speed = 0.0;
  m_max_v_speed = 0.0;
  m_max_w_speed = 0.0;

  maxdt_temporary = 0.0;
  m_dt              = 0.0;
  skipCountDown   = 0;

  m_timestep_hit_multiples_last_time = m_time->current();

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

  delete m_stress_balance;

  if (not m_external_ocean_model) {
    delete m_ocean;
  }

  if (not m_external_surface_model) {
    delete m_surface;
  }

  delete m_beddef;

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

  m_ice_enthalpy.create(m_grid, "enthalpy", WITH_GHOSTS, WIDE_STENCIL);
  // POSSIBLE standard name = land_ice_enthalpy
  m_ice_enthalpy.set_attrs("model_state",
                  "ice enthalpy (includes sensible heat, latent heat, pressure)",
                  "J kg-1", "");
  m_grid->variables().add(m_ice_enthalpy);

  if (m_config->get_boolean("do_cold_ice_methods")) {
    // ice temperature
    m_ice_temperature.create(m_grid, "temp", WITH_GHOSTS);
    m_ice_temperature.set_attrs("model_state",
                 "ice temperature", "K", "land_ice_temperature");
    m_ice_temperature.metadata().set_double("valid_min", 0.0);
    m_grid->variables().add(m_ice_temperature);

    if (m_config->get_boolean("do_energy")) {
      m_ice_enthalpy.metadata().set_string("pism_intent", "diagnostic");
    } else {
      m_ice_temperature.metadata().set_string("pism_intent", "diagnostic");
    }
  }

  // age of ice but only if age will be computed
  if (m_config->get_boolean("do_age")) {
    m_ice_age.create(m_grid, "age", WITH_GHOSTS, WIDE_STENCIL);
    // PROPOSED standard_name = land_ice_age
    m_ice_age.set_attrs("model_state", "age of ice",
                   "s", "");
    m_ice_age.metadata().set_string("glaciological_units", "years");
    m_ice_age.write_in_glaciological_units = true;
    m_ice_age.metadata().set_double("valid_min", 0.0);
    m_grid->variables().add(m_ice_age);
  }

  // ice upper surface elevation
  m_ice_surface_elevation.create(m_grid, "usurf", WITH_GHOSTS, WIDE_STENCIL);
  m_ice_surface_elevation.set_attrs("diagnostic", "ice upper surface elevation",
                                  "m", "surface_altitude");
  m_grid->variables().add(m_ice_surface_elevation);

  // land ice thickness
  m_ice_thickness.create(m_grid, "thk", WITH_GHOSTS, WIDE_STENCIL);
  m_ice_thickness.set_attrs("model_state", "land ice thickness",
                          "m", "land_ice_thickness");
  m_ice_thickness.metadata().set_double("valid_min", 0.0);
  m_grid->variables().add(m_ice_thickness);

  if (m_config->get_boolean("sub_groundingline")) {
    m_gl_mask.create(m_grid, "gl_mask", WITHOUT_GHOSTS);
    m_gl_mask.set_attrs("internal",
                      "fractional grounded/floating mask (floating=0, grounded=1)",
                      "", "");
    m_grid->variables().add(m_gl_mask);
  }

  // grounded_dragging_floating integer mask
  m_cell_type.create(m_grid, "mask", WITH_GHOSTS, WIDE_STENCIL);
  m_cell_type.set_attrs("diagnostic", "ice-type (ice-free/grounded/floating/ocean) integer mask",
                  "", "");
  std::vector<double> mask_values(4);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_GROUNDED;
  mask_values[2] = MASK_FLOATING;
  mask_values[3] = MASK_ICE_FREE_OCEAN;
  m_cell_type.metadata().set_doubles("flag_values", mask_values);
  m_cell_type.metadata().set_string("flag_meanings",
                              "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
  m_grid->variables().add(m_cell_type);

  // upward geothermal flux at bedrock surface
  m_geothermal_flux.create(m_grid, "bheatflx", WITHOUT_GHOSTS);
  // PROPOSED standard_name = lithosphere_upward_heat_flux
  m_geothermal_flux.set_attrs("climate_steady", "upward geothermal flux at bedrock surface",
                            "W m-2", "");
  m_geothermal_flux.metadata().set_string("glaciological_units", "mW m-2");
  m_geothermal_flux.write_in_glaciological_units = true;
  m_geothermal_flux.set_time_independent(true);
  m_grid->variables().add(m_geothermal_flux);

  // temperature seen by top of bedrock thermal layer
  m_bedtoptemp.create(m_grid, "bedtoptemp", WITHOUT_GHOSTS);
  m_bedtoptemp.set_attrs("internal",
                       "temperature of top of bedrock thermal layer",
                       "K", "");
  m_bedtoptemp.metadata().set_string("glaciological_units", "K");
  m_grid->variables().add(m_bedtoptemp);

  // yield stress for basal till (plastic or pseudo-plastic model)
  {
    m_basal_yield_stress.create(m_grid, "tauc", WITH_GHOSTS, WIDE_STENCIL);
    // PROPOSED standard_name = land_ice_basal_material_yield_stress
    m_basal_yield_stress.set_attrs("diagnostic",
                                 "yield stress for basal till (plastic or pseudo-plastic model)",
                                 "Pa", "");
    m_grid->variables().add(m_basal_yield_stress);
  }

  // basal melt rate
  m_basal_melt_rate.create(m_grid, "bmelt", WITHOUT_GHOSTS);
  // ghosted to allow the "redundant" computation of tauc
  m_basal_melt_rate.set_attrs("model_state",
                            "ice basal melt rate from energy conservation and subshelf melt, in ice thickness per time",
                            "m s-1", "land_ice_basal_melt_rate");
  m_basal_melt_rate.metadata().set_string("glaciological_units", "m year-1");
  m_basal_melt_rate.write_in_glaciological_units = true;
  m_basal_melt_rate.metadata().set_string("comment", "positive basal melt rate corresponds to ice loss");
  m_grid->variables().add(m_basal_melt_rate);

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

  if (m_config->get_boolean("part_grid")) {
    // Href
    vHref.create(m_grid, "Href", WITH_GHOSTS);
    vHref.set_attrs("model_state", "temporary ice thickness at calving front boundary",
                    "m", "");
    m_grid->variables().add(vHref);
  }

  if (m_config->get_string("calving_methods").find("eigen_calving") != std::string::npos or
      m_config->get_boolean("do_fracture_density")) {

    m_strain_rates.create(m_grid, "strain_rates", WITH_GHOSTS,
                        2, // stencil width, has to match or exceed the "offset" in eigenCalving
                        2);

    m_strain_rates.metadata(0).set_name("eigen1");
    m_strain_rates.set_attrs("internal",
                           "major principal component of horizontal strain-rate",
                           "second-1", "", 0);

    m_strain_rates.metadata(1).set_name("eigen2");
    m_strain_rates.set_attrs("internal",
                           "minor principal component of horizontal strain-rate",
                           "second-1", "", 1);
  }

  if (m_config->get_boolean("do_fracture_density")) {

    m_deviatoric_stresses.create(m_grid, "sigma", WITH_GHOSTS,
                               2, // stencil width
                               3);

    m_deviatoric_stresses.metadata(0).set_name("sigma_xx");
    m_deviatoric_stresses.set_attrs("internal",
                                  "deviatoric stress in x direction",
                                  "Pa", "", 0);

    m_deviatoric_stresses.metadata(1).set_name("sigma_yy");
    m_deviatoric_stresses.set_attrs("internal",
                                  "deviatoric stress in y direction",
                                  "Pa", "", 1);

    m_deviatoric_stresses.metadata(2).set_name("sigma_xy");
    m_deviatoric_stresses.set_attrs("internal",
                                  "deviatoric shear stress",
                                  "Pa", "", 2);
  }

  if (m_config->get_boolean("ssa_dirichlet_bc")) {
    // bc_locations
    m_ssa_dirichlet_bc_mask.create(m_grid, "bc_mask", WITH_GHOSTS, WIDE_STENCIL);
    m_ssa_dirichlet_bc_mask.set_attrs("model_state", "Dirichlet boundary mask",
                      "", "");
    std::vector<double> bc_mask_values(2);
    bc_mask_values[0] = 0;
    bc_mask_values[1] = 1;
    m_ssa_dirichlet_bc_mask.metadata().set_doubles("flag_values", bc_mask_values);
    m_ssa_dirichlet_bc_mask.metadata().set_string("flag_meanings", "no_data bc_condition");
    m_grid->variables().add(m_ssa_dirichlet_bc_mask);

    double fill_value = units::convert(m_sys, m_config->get_double("fill_value"),
                                       "m year-1", "m second-1");
    // vel_bc
    m_ssa_dirichlet_bc_values.create(m_grid, "_ssa_bc", WITH_GHOSTS, WIDE_STENCIL); // u_ssa_bc and v_ssa_bc
    m_ssa_dirichlet_bc_values.set_attrs("model_state",
                     "X-component of the SSA velocity boundary conditions",
                     "m s-1", "", 0);
    m_ssa_dirichlet_bc_values.set_attrs("model_state",
                     "Y-component of the SSA velocity boundary conditions",
                     "m s-1", "", 1);
    for (int j = 0; j < 2; ++j) {
      m_ssa_dirichlet_bc_values.metadata(j).set_string("glaciological_units", "m year-1");
      m_ssa_dirichlet_bc_values.metadata(j).set_double("valid_min",  units::convert(m_sys, -1e6, "m year-1", "m second-1"));
      m_ssa_dirichlet_bc_values.metadata(j).set_double("valid_max",  units::convert(m_sys,  1e6, "m year-1", "m second-1"));
      m_ssa_dirichlet_bc_values.metadata(j).set_double("_FillValue", fill_value);
    }
    //just for diagnostics...
    m_grid->variables().add(m_ssa_dirichlet_bc_values);
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
      vFG.set_attrs("model_state", "fracture growth rate",       "second-1", "");
      vFG.metadata().set_double("valid_min", 0.0);
      m_grid->variables().add(vFG);

      vFH.create(m_grid, "fracture_healing_rate", WITH_GHOSTS, WIDE_STENCIL);
      vFH.set_attrs("model_state", "fracture healing rate",      "second-1", "");
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
  m_cell_area.create(m_grid, "cell_area", WITHOUT_GHOSTS);
  m_cell_area.set_attrs("diagnostic", "cell areas", "m2", "");
  m_cell_area.metadata().set_string("comment",
                                  "values are equal to dx*dy "
                                  "if projection parameters are not available; "
                                  "otherwise WGS84 ellipsoid is used");
  m_cell_area.set_time_independent(true);
  m_cell_area.metadata().set_string("glaciological_units", "km2");
  m_cell_area.write_in_glaciological_units = true;
  m_grid->variables().add(m_cell_area);

  // fields owned by IceModel but filled by SurfaceModel *surface:
  // mean annual net ice equivalent surface mass balance rate
  m_climatic_mass_balance.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  m_climatic_mass_balance.set_attrs("climate_from_SurfaceModel",  // FIXME: can we do better?
                                  "ice-equivalent surface mass balance (accumulation/ablation) rate",
                                  "kg m-2 s-1",
                                  "land_ice_surface_specific_mass_balance_flux");
  m_climatic_mass_balance.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_climatic_mass_balance.write_in_glaciological_units = true;
  m_climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // annual mean air temperature at "ice surface", at level below all firn
  //   processes (e.g. "10 m" or ice temperatures)
  m_ice_surface_temp.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
  m_ice_surface_temp.set_attrs("climate_from_SurfaceModel",  // FIXME: can we do better?
                             "annual average ice surface temperature, below firn processes",
                             "K",
                             "");

  m_liqfrac_surface.create(m_grid, "liqfrac_surface", WITHOUT_GHOSTS);
  m_liqfrac_surface.set_attrs("climate_from_SurfaceModel",
                            "liquid water fraction at the top surface of the ice",
                            "1", "");
  // grid.variables().add(liqfrac_surface);

  // ice mass balance rate at the base of the ice shelf; sign convention for
  //   vshelfbasemass matches standard sign convention for basal melt rate of
  //   grounded ice
  m_shelfbmassflux.create(m_grid, "shelfbmassflux", WITHOUT_GHOSTS); // no ghosts; NO HOR. DIFF.!
  m_shelfbmassflux.set_attrs("climate_state", "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                           "m s-1", "");
  // PROPOSED standard name = ice_shelf_basal_specific_mass_balance
  // rescales from m second-1 to m year-1 when writing to NetCDF and std out:
  m_shelfbmassflux.write_in_glaciological_units = true;
  m_shelfbmassflux.metadata().set_string("glaciological_units", "m year-1");
  // do not add; boundary models are in charge here
  // grid.variables().add(shelfbmassflux);

  // ice boundary tempature at the base of the ice shelf
  m_shelfbtemp.create(m_grid, "shelfbtemp", WITHOUT_GHOSTS); // no ghosts; NO HOR. DIFF.!
  m_shelfbtemp.set_attrs("climate_state", "absolute temperature at ice shelf base",
                       "K", "");
  // PROPOSED standard name = ice_shelf_basal_temperature
  // do not add; boundary models are in charge here
  // grid.variables().add(shelfbtemp);

  {
    m_climatic_mass_balance_cumulative.create(m_grid,
                                              "climatic_mass_balance_cumulative",
                                              WITHOUT_GHOSTS);
    m_climatic_mass_balance_cumulative.set_attrs("diagnostic",
                                                 "cumulative surface mass balance",
                                                 "kg m-2", "");
  }

  {
    m_flux_divergence.create(m_grid, "flux_divergence", WITHOUT_GHOSTS);
    m_flux_divergence.set_attrs("diagnostic",
                                "flux divergence",
                                "m s-1", "");
    m_flux_divergence.metadata().set_string("glaciological_units", "m year-1");
    m_flux_divergence.write_in_glaciological_units = true;
  }

  {
    m_grounded_basal_flux_2D_cumulative.create(m_grid, "grounded_basal_flux_cumulative", WITHOUT_GHOSTS);
    m_grounded_basal_flux_2D_cumulative.set_attrs("diagnostic",
                                                  "cumulative basal flux into the ice "
                                                  "in grounded areas (positive means ice gain)",
                                                  "kg m-2", "");
    m_grounded_basal_flux_2D_cumulative.set_time_independent(false);
    m_grounded_basal_flux_2D_cumulative.metadata().set_string("glaciological_units", "Gt m-2");
    m_grounded_basal_flux_2D_cumulative.write_in_glaciological_units = true;
  }

  {
    m_floating_basal_flux_2D_cumulative.create(m_grid, "floating_basal_flux_cumulative", WITHOUT_GHOSTS);
    m_floating_basal_flux_2D_cumulative.set_attrs("diagnostic",
                                                  "cumulative basal flux into the ice "
                                                  "in floating areas (positive means ice gain)",
                                                  "kg m-2", "");
    m_floating_basal_flux_2D_cumulative.set_time_independent(false);
    m_floating_basal_flux_2D_cumulative.metadata().set_string("glaciological_units", "Gt m-2");
    m_floating_basal_flux_2D_cumulative.write_in_glaciological_units = true;
  }

  {
    m_nonneg_flux_2D_cumulative.create(m_grid, "nonneg_flux_cumulative", WITHOUT_GHOSTS);
    m_nonneg_flux_2D_cumulative.set_attrs("diagnostic",
                                          "cumulative nonnegative rule flux (positive means ice gain)",
                                          "kg m-2", "");
    m_nonneg_flux_2D_cumulative.set_time_independent(false);
    m_nonneg_flux_2D_cumulative.metadata().set_string("glaciological_units", "Gt m-2");
    m_nonneg_flux_2D_cumulative.write_in_glaciological_units = true;
  }

  {
    m_discharge_flux_2D_cumulative.create(m_grid, "discharge_flux_cumulative", WITHOUT_GHOSTS);
    m_discharge_flux_2D_cumulative.set_attrs("diagnostic",
                                             "cumulative discharge (calving) flux (positive means ice loss)",
                                             "kg m-2", "");
    m_discharge_flux_2D_cumulative.set_time_independent(false);
    m_discharge_flux_2D_cumulative.metadata().set_string("glaciological_units", "Gt m-2");
    m_discharge_flux_2D_cumulative.write_in_glaciological_units = true;
  }
}

//! The contents of the main PISM time-step.
/*!
During the time-step we perform the following actions:
 */
void IceModel::step(bool do_mass_continuity,
                              bool do_energy,
                              bool do_age,
                              bool do_skip) {

  const Profiling &profiling = m_ctx->profiling();

  double current_time = m_time->current();

  //! \li call additionalAtStartTimestep() to let derived classes do more
  additionalAtStartTimestep();  // might set maxdt_temporary

  //! \li update the velocity field; in some cases the whole three-dimensional
  //! field is updated and in some cases just the vertically-averaged
  //! horizontal velocity is updated

  // always "update" ice velocity (possibly trivially); only update
  // SSA and only update velocities at depth if suggested by temp and age
  // stability criterion; note *lots* of communication is avoided by skipping
  // SSA (and temp/age)

  bool updateAtDepth = (skipCountDown == 0),
    do_energy_step = updateAtDepth and do_energy;

  //! \li update the yield stress for the plastic till model (if appropriate)
  if (updateAtDepth and basal_yield_stress_model) {
    profiling.begin("basal yield stress");
    basal_yield_stress_model->update(current_time, m_dt);
    profiling.end("basal yield stress");
    m_basal_yield_stress.copy_from(basal_yield_stress_model->basal_material_yield_stress());
    stdout_flags += "y";
  } else {
    stdout_flags += "$";
  }

  // Update the fractional grounded/floating mask (used by the SSA
  // stress balance and the energy code)
  if (m_config->get_boolean("sub_groundingline")) {
    updateSurfaceElevationAndMask(); // update h and mask
    update_grounded_cell_fraction();
  }

  double sea_level = m_ocean->sea_level_elevation();

  IceModelVec2S &melange_back_pressure = vWork2d[0];

  m_ocean->melange_back_pressure_fraction(melange_back_pressure);

  try {
    profiling.begin("stress balance");
    m_stress_balance->update(not updateAtDepth,
                           sea_level,
                           melange_back_pressure);
    profiling.end("stress balance");
  } catch (RuntimeError &e) {
    options::String output_file("-o", "output file name",
                                "output.nc", options::DONT_ALLOW_EMPTY);

    std::string o_file = pism_filename_add_suffix(output_file,
                                                  "_stressbalance_failed", "");
    dumpToFile(o_file);

    e.add_context("performing a time step. (Note: Model state was saved to '%s'.)",
                  o_file.c_str());
    throw;
  }

  stdout_flags += m_stress_balance->stdout_report();

  stdout_flags += (updateAtDepth ? "v" : "V");

  //! \li determine the time step according to a variety of stability criteria;
  //!  see determineTimeStep()
  max_timestep(m_dt, skipCountDown);

  //! \li Update surface and ocean models.
  profiling.begin("surface");
  m_surface->update(current_time, m_dt);
  profiling.end("surface");

  profiling.begin("ocean");
  m_ocean->update(current_time, m_dt);
  profiling.end("ocean");

  dt_TempAge += m_dt;
  // IceModel::dt,dtTempAge are now set correctly according to
  // mass-continuity-eqn-diffusivity criteria, horizontal CFL criteria, and
  // other criteria from derived class additionalAtStartTimestep(), and from
  // "-skip" mechanism

  //! \li update the age of the ice (if appropriate)
  if (do_age and updateAtDepth) {
    profiling.begin("age");
    ageStep();
    profiling.end("age");
    stdout_flags += "a";
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
  subglacial_hydrology->update(current_time, m_dt);
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
    // if do_mass_continuity is true.
    if (do_skip and skipCountDown > 0) {
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
  if (m_beddef) {
    const IceModelVec2S &bed_topography = m_beddef->bed_elevation();
    int topg_state_counter = bed_topography.get_state_counter();

    profiling.begin("bed deformation");
    m_beddef->update(current_time, m_dt);
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
  m_time->step(m_dt);

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
    maxdt_temporary = -1.0;

    step(do_mass_conserve, do_energy, do_age, do_skip);

    // report a summary for major steps or the last one
    bool updateAtDepth = skipCountDown == 0;
    bool tempAgeStep = updateAtDepth and (do_energy or do_age);

    const bool show_step = tempAgeStep or m_adaptive_timestep_reason == "end of the run";
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

  //! The following flow-chart illustrates the process.
  //!
  //! \dotfile initialization-sequence.dot "IceModel initialization sequence"

  // Get the start time in seconds and ensure that it is consistent
  // across all processors.
  m_start_time = GlobalMax(m_grid->com, GetTime());

  profiling.end("initialization");
}

// FIXME: THIS IS BAD! (Provides unguarded access to IceModel's internals.)
IceModelVec2S* IceModel::get_geothermal_flux() {
  return &this->m_geothermal_flux;
}

// FIXME: THIS IS BAD! (Provides unguarded access to IceModel's internals.)
stressbalance::StressBalance* IceModel::get_stress_balance() {
  return this->m_stress_balance;
}

} // end of namespace pism

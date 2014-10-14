// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <petscdmda.h>

#include "iceModel.hh"
#include "pism_signal.h"
#include "PISMStressBalance.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "PISMBedDef.hh"
#include "bedrockThermalUnit.hh"
#include "PISMHydrology.hh"
#include "PISMYieldStress.hh"
#include "basal_resistance.hh"
#include "enthalpyConverter.hh"
#include "pism_options.hh"
#include "IceGrid.hh"
#include "PISMDiagnostic.hh"
#include "PISMIcebergRemover.hh"
#include "PISMOceanKill.hh"
#include "PISMFloatKill.hh"
#include "PISMCalvingAtThickness.hh"
#include "PISMEigenCalving.hh"
#include "Mask.hh"

namespace pism {

IceModel::IceModel(IceGrid &g, Config &conf, Config &conf_overrides)
  : grid(g),
    config(conf),
    overrides(conf_overrides),
    global_attributes("PISM_GLOBAL", g.get_unit_system()),
    mapping("mapping", g.get_unit_system()),
    run_stats("run_stats", g.get_unit_system()),
    extra_bounds("time_bounds", config.get_string("time_dimension_name"), g.get_unit_system()),
    timestamp("timestamp", config.get_string("time_dimension_name"), g.get_unit_system()) {

  extra_bounds.set_units(grid.time->units_string());

  timestamp.set_units("hours");
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

  EC  = NULL;
  btu = NULL;

  iceberg_remover             = NULL;
  ocean_kill_calving          = NULL;
  float_kill_calving          = NULL;
  thickness_threshold_calving = NULL;
  eigen_calving               = NULL;

  executable_short_name = "pism"; // drivers typically override this

  // initializr maximum |u|,|v|,|w| in ice
  gmaxu = 0;
  gmaxv = 0;
  gmaxw = 0;

  // set default locations of the column used by -view_system
  id = (grid.Mx - 1)/2;
  jd = (grid.My - 1)/2;

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

  timestep_hit_multiples_last_time = grid.time->current();

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
  while (i != ts_diagnostics.end()) delete (i++)->second;

  // de-allocate diagnostics
  std::map<std::string,Diagnostic*>::iterator j = diagnostics.begin();
  while (j != diagnostics.end()) delete (j++)->second;


  // de-allocate viewers
  std::map<std::string,PetscViewer>::iterator k = viewers.begin();
  while (k != viewers.end()) {
    if ((*k).second != NULL) {
      PetscViewerDestroy(&(*k).second);
      ++k;
    }
  }

  delete stress_balance;

  if (external_ocean_model == false)
    delete ocean;

  if (external_surface_model == false)
    delete surface;

  delete beddef;

  delete subglacial_hydrology;
  delete basal_yield_stress_model;
  delete EC;
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
PetscErrorCode IceModel::createVecs() {
  PetscErrorCode ierr;

  const unsigned int WIDE_STENCIL = config.get("grid_max_stencil_width");

  ierr = verbPrintf(3, grid.com,
                    "Allocating memory...\n"); CHKERRQ(ierr);

  // get the list of selected calving methods:
  std::istringstream calving_methods_list(config.get_string("calving_methods"));
  std::string calving_method_name;
  std::set<std::string> calving_methods;

  while (getline(calving_methods_list, calving_method_name, ','))
    calving_methods.insert(calving_method_name);

  // The following code creates (and documents -- to some extent) the
  // variables. The main (and only) principle here is using standard names from
  // the CF conventions; see
  // http://cf-pcmdi.llnl.gov/documents/cf-standard-names

  ierr = Enth3.create(grid, "enthalpy", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
  // POSSIBLE standard name = land_ice_enthalpy
  ierr = Enth3.set_attrs("model_state",
                         "ice enthalpy (includes sensible heat, latent heat, pressure)",
                         "J kg-1", ""); CHKERRQ(ierr);
  ierr = variables.add(Enth3); CHKERRQ(ierr);

  if (config.get_flag("do_cold_ice_methods")) {
    // ice temperature
    ierr = T3.create(grid, "temp", WITH_GHOSTS); CHKERRQ(ierr);
    ierr = T3.set_attrs("model_state",
                        "ice temperature", "K", "land_ice_temperature"); CHKERRQ(ierr);
    T3.metadata().set_double("valid_min", 0.0);
    ierr = variables.add(T3); CHKERRQ(ierr);

    if (config.get_flag("do_energy") == true) {
      Enth3.metadata().set_string("pism_intent", "diagnostic");
    } else {
      T3.metadata().set_string("pism_intent", "diagnostic");
    }
  }

  // age of ice but only if age will be computed
  if (config.get_flag("do_age")) {
    ierr = tau3.create(grid, "age", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
    // PROPOSED standard_name = land_ice_age
    ierr = tau3.set_attrs("model_state", "age of ice",
                          "s", ""); CHKERRQ(ierr);
    ierr = tau3.set_glaciological_units("years");
    tau3.write_in_glaciological_units = true;
    tau3.metadata().set_double("valid_min", 0.0);
    ierr = variables.add(tau3); CHKERRQ(ierr);
  }

  // ice upper surface elevation
  ierr = ice_surface_elevation.create(grid, "usurf", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = ice_surface_elevation.set_attrs("diagnostic", "ice upper surface elevation",
                      "m", "surface_altitude"); CHKERRQ(ierr);
  ierr = variables.add(ice_surface_elevation); CHKERRQ(ierr);

  // land ice thickness
  ierr = ice_thickness.create(grid, "thk", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = ice_thickness.set_attrs("model_state", "land ice thickness",
                                 "m", "land_ice_thickness"); CHKERRQ(ierr);
  ice_thickness.metadata().set_double("valid_min", 0.0);
  ierr = variables.add(ice_thickness); CHKERRQ(ierr);

  // bedrock surface elevation
  ierr = bed_topography.create(grid, "topg", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = bed_topography.set_attrs("model_state", "bedrock surface elevation",
                                  "m", "bedrock_altitude"); CHKERRQ(ierr);
  ierr = variables.add(bed_topography); CHKERRQ(ierr);

  if (config.get_flag("sub_groundingline")) {
    ierr = gl_mask.create(grid, "gl_mask", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = gl_mask.set_attrs("internal",
                             "fractional grounded/floating mask (floating=0, grounded=1)",
                             "", ""); CHKERRQ(ierr);
    ierr = variables.add(gl_mask); CHKERRQ(ierr);

    ierr = gl_mask_x.create(grid, "gl_mask_x", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = gl_mask_x.set_attrs("internal",
                               "fractional grounded/floating mask in x-direction (floating=0, grounded=1)",
                               "", ""); CHKERRQ(ierr);
    ierr = variables.add(gl_mask_x); CHKERRQ(ierr);

    ierr = gl_mask_y.create(grid, "gl_mask_y", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = gl_mask_y.set_attrs("internal",
                               "fractional grounded/floating mask in y-direction (floating=0, grounded=1)",
                               "", ""); CHKERRQ(ierr);
    ierr = variables.add(gl_mask_y); CHKERRQ(ierr);
  }

  // grounded_dragging_floating integer mask
  ierr = vMask.create(grid, "mask", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = vMask.set_attrs("diagnostic", "ice-type (ice-free/grounded/floating/ocean) integer mask",
                         "", ""); CHKERRQ(ierr);
  std::vector<double> mask_values(4);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_GROUNDED;
  mask_values[2] = MASK_FLOATING;
  mask_values[3] = MASK_ICE_FREE_OCEAN;
  vMask.metadata().set_doubles("flag_values", mask_values);
  vMask.metadata().set_string("flag_meanings",
                              "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
  ierr = variables.add(vMask); CHKERRQ(ierr);

  // upward geothermal flux at bedrock surface
  ierr = geothermal_flux.create(grid, "bheatflx", WITHOUT_GHOSTS); CHKERRQ(ierr);
  // PROPOSED standard_name = lithosphere_upward_heat_flux
  ierr = geothermal_flux.set_attrs("climate_steady", "upward geothermal flux at bedrock surface",
                        "W m-2", ""); CHKERRQ(ierr);
  ierr = geothermal_flux.set_glaciological_units("mW m-2");
  geothermal_flux.write_in_glaciological_units = true;
  geothermal_flux.set_time_independent(true);
  ierr = variables.add(geothermal_flux); CHKERRQ(ierr);

  // temperature seen by top of bedrock thermal layer
  ierr = bedtoptemp.create(grid, "bedtoptemp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = bedtoptemp.set_attrs("internal",
                              "temperature of top of bedrock thermal layer",
                              "K", ""); CHKERRQ(ierr);
  ierr = bedtoptemp.set_glaciological_units("K");
  ierr = variables.add(bedtoptemp); CHKERRQ(ierr);

  // yield stress for basal till (plastic or pseudo-plastic model)
  {
    ierr = basal_yield_stress.create(grid, "tauc", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
    // PROPOSED standard_name = land_ice_basal_material_yield_stress
    ierr = basal_yield_stress.set_attrs("diagnostic",
                                        "yield stress for basal till (plastic or pseudo-plastic model)",
                                        "Pa", ""); CHKERRQ(ierr);
    ierr = variables.add(basal_yield_stress); CHKERRQ(ierr);
  }

  // bedrock uplift rate
  ierr = bed_uplift_rate.create(grid, "dbdt", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = bed_uplift_rate.set_attrs("model_state", "bedrock uplift rate",
                                   "m s-1", "tendency_of_bedrock_altitude"); CHKERRQ(ierr);
  ierr = bed_uplift_rate.set_glaciological_units("m year-1");
  bed_uplift_rate.write_in_glaciological_units = true;
  ierr = variables.add(bed_uplift_rate); CHKERRQ(ierr);

  // basal melt rate
  ierr = basal_melt_rate.create(grid, "bmelt", WITHOUT_GHOSTS); CHKERRQ(ierr);
  // ghosted to allow the "redundant" computation of tauc
  ierr = basal_melt_rate.set_attrs("model_state",
                                   "ice basal melt rate from energy conservation and subshelf melt, in ice thickness per time",
                                   "m s-1", "land_ice_basal_melt_rate"); CHKERRQ(ierr);
  ierr = basal_melt_rate.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  basal_melt_rate.write_in_glaciological_units = true;
  basal_melt_rate.metadata().set_string("comment", "positive basal melt rate corresponds to ice loss");
  ierr = variables.add(basal_melt_rate); CHKERRQ(ierr);

  // longitude
  ierr = vLongitude.create(grid, "lon", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = vLongitude.set_attrs("mapping", "longitude", "degree_east", "longitude"); CHKERRQ(ierr);
  vLongitude.set_time_independent(true);
  vLongitude.metadata().set_string("coordinates", "");
  vLongitude.metadata().set_string("grid_mapping", "");
  vLongitude.metadata().set_double("valid_min", -180.0);
  vLongitude.metadata().set_double("valid_max",  180.0);
  ierr = variables.add(vLongitude); CHKERRQ(ierr);

  // latitude
  ierr = vLatitude.create(grid, "lat", WITH_GHOSTS); CHKERRQ(ierr); // has ghosts so that we can compute cell areas
  ierr = vLatitude.set_attrs("mapping", "latitude", "degree_north", "latitude"); CHKERRQ(ierr);
  vLatitude.set_time_independent(true);
  vLatitude.metadata().set_string("coordinates", "");
  vLatitude.metadata().set_string("grid_mapping", "");
  vLatitude.metadata().set_double("valid_min", -90.0);
  vLatitude.metadata().set_double("valid_max",  90.0);
  ierr = variables.add(vLatitude); CHKERRQ(ierr);

  if (config.get_flag("part_grid") == true) {
    // Href
    ierr = vHref.create(grid, "Href", WITH_GHOSTS); CHKERRQ(ierr);
    ierr = vHref.set_attrs("model_state", "temporary ice thickness at calving front boundary",
                           "m", ""); CHKERRQ(ierr);
    ierr = variables.add(vHref); CHKERRQ(ierr);
  }

  if (config.get_string("calving_methods").find("eigen_calving") != std::string::npos ||
      config.get_flag("do_fracture_density") == true) {

    ierr = strain_rates.create(grid, "edot", WITH_GHOSTS,
                               2, // stencil width, has to match or exceed the "offset" in eigenCalving
                               2); CHKERRQ(ierr);

    ierr = strain_rates.set_name("edot_1", 0); CHKERRQ(ierr);
    ierr = strain_rates.set_attrs("internal",
                                  "major principal component of horizontal strain-rate",
                                  "1/s", "", 0); CHKERRQ(ierr);

    ierr = strain_rates.set_name("edot_2", 1); CHKERRQ(ierr);
    ierr = strain_rates.set_attrs("internal",
                                  "minor principal component of horizontal strain-rate",
                                  "1/s", "", 1); CHKERRQ(ierr);
  }

  if (config.get_flag("do_fracture_density") == true) {
    
    ierr = deviatoric_stresses.create(grid, "sigma", WITH_GHOSTS,
                                      2, // stencil width
                                      3); CHKERRQ(ierr);
    
    ierr = deviatoric_stresses.set_name("sigma_xx", 0); CHKERRQ(ierr);
    ierr = deviatoric_stresses.set_attrs("internal",
                                         "deviatoric stress in x direction",
                                         "Pa", "", 0); CHKERRQ(ierr);
                                  
    ierr = deviatoric_stresses.set_name("sigma_yy", 1); CHKERRQ(ierr);
    ierr = deviatoric_stresses.set_attrs("internal",
                                         "deviatoric stress in y direction",
                                         "Pa", "", 1); CHKERRQ(ierr);   
                                         
    ierr = deviatoric_stresses.set_name("sigma_xy", 2); CHKERRQ(ierr);
    ierr = deviatoric_stresses.set_attrs("internal",
                                         "deviatoric shear stress",
                                         "Pa", "", 2); CHKERRQ(ierr);
  }

  if (config.get_flag("ssa_dirichlet_bc") == true) {
    // bc_locations
    ierr = vBCMask.create(grid, "bcflag", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vBCMask.set_attrs("model_state", "Dirichlet boundary mask",
                             "", ""); CHKERRQ(ierr);
    std::vector<double> bc_mask_values(2);
    bc_mask_values[0] = 0;
    bc_mask_values[1] = 1;
    vBCMask.metadata().set_doubles("flag_values", bc_mask_values);
    vBCMask.metadata().set_string("flag_meanings", "no_data bc_condition");
    ierr = variables.add(vBCMask); CHKERRQ(ierr);


    // vel_bc
    ierr = vBCvel.create(grid, "_ssa_bc", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr); // u_ssa_bc and v_ssa_bc
    ierr = vBCvel.set_attrs("model_state",
                            "X-component of the SSA velocity boundary conditions",
                            "m s-1", "", 0); CHKERRQ(ierr);
    ierr = vBCvel.set_attrs("model_state",
                            "Y-component of the SSA velocity boundary conditions",
                            "m s-1", "", 1); CHKERRQ(ierr);
    ierr = vBCvel.set_glaciological_units("m year-1"); CHKERRQ(ierr);
    for (int j = 0; j < 2; ++j) {
      vBCvel.metadata(j).set_double("valid_min",  grid.convert(-1e6, "m/year", "m/second"));
      vBCvel.metadata(j).set_double("valid_max",  grid.convert( 1e6, "m/year", "m/second"));
      vBCvel.metadata(j).set_double("_FillValue", config.get("fill_value", "m/year", "m/s"));
    }
    //just for diagnostics...
    ierr = variables.add(vBCvel); CHKERRQ(ierr);
  }

  // fracture density field
  if (config.get_flag("do_fracture_density")) {
    ierr = vFD.create(grid, "fracture_density", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr); 
    ierr = vFD.set_attrs("model_state", "fracture density in ice shelf", "", ""); CHKERRQ(ierr);
    vFD.metadata().set_double("valid_max", 1.0);
    vFD.metadata().set_double("valid_min", 0.0);
    ierr = variables.add(vFD); CHKERRQ(ierr);

    if (config.get_flag("write_fd_fields")) {
      ierr = vFG.create(grid, "fracture_growth_rate", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr); 
      ierr = vFG.set_attrs("model_state", "fracture growth rate",       "1/s", ""); CHKERRQ(ierr);
      vFG.metadata().set_double("valid_min", 0.0);
      ierr = variables.add(vFG); CHKERRQ(ierr);

      ierr = vFH.create(grid, "fracture_healing_rate", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr); 
      ierr = vFH.set_attrs("model_state", "fracture healing rate",      "1/s", ""); CHKERRQ(ierr);
      ierr = variables.add(vFH); CHKERRQ(ierr);

      ierr = vFE.create(grid, "fracture_flow_enhancement", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr); 
      ierr = vFE.set_attrs("model_state", "fracture-induced flow enhancement", "", ""); CHKERRQ(ierr);
      ierr = variables.add(vFE); CHKERRQ(ierr);

      ierr = vFA.create(grid, "fracture_age", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr); 
      ierr = vFA.set_attrs("model_state", "age since fracturing",       "years", ""); CHKERRQ(ierr);
      ierr = variables.add(vFA); CHKERRQ(ierr);
      
      ierr = vFT.create(grid, "fracture_toughness", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr); 
      ierr = vFT.set_attrs("model_state", "fracture toughness", "Pa", ""); CHKERRQ(ierr);
      ierr = variables.add(vFT); CHKERRQ(ierr);
    }
  }

  // cell areas
  ierr = cell_area.create(grid, "cell_area", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = cell_area.set_attrs("diagnostic", "cell areas", "m2", ""); CHKERRQ(ierr);
  cell_area.metadata().set_string("comment",
                                  "values are equal to dx*dy "
                                  "if projection parameters are not available; "
                                  "otherwise WGS84 ellipsoid is used");
  cell_area.set_time_independent(true);
  ierr = cell_area.set_glaciological_units("km2"); CHKERRQ(ierr);
  cell_area.write_in_glaciological_units = true;
  ierr = variables.add(cell_area); CHKERRQ(ierr);

  // fields owned by IceModel but filled by SurfaceModel *surface:
  // mean annual net ice equivalent surface mass balance rate
  ierr = climatic_mass_balance.create(grid, "climatic_mass_balance", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = climatic_mass_balance.set_attrs(
                        "climate_from_SurfaceModel",  // FIXME: can we do better?
                        "ice-equivalent surface mass balance (accumulation/ablation) rate",
                        "kg m-2 s-1",
                        "land_ice_surface_specific_mass_balance_flux");  // CF standard_name
  CHKERRQ(ierr);
  ierr = climatic_mass_balance.set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);
  climatic_mass_balance.write_in_glaciological_units = true;
  climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // annual mean air temperature at "ice surface", at level below all firn
  //   processes (e.g. "10 m" or ice temperatures)
  ierr = ice_surface_temp.create(grid, "ice_surface_temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = ice_surface_temp.set_attrs(
                                    "climate_from_SurfaceModel",  // FIXME: can we do better?
                                    "annual average ice surface temperature, below firn processes",
                                    "K",
                                    "");  // PROPOSED CF standard_name = land_ice_surface_temperature_below_firn
  CHKERRQ(ierr);

  ierr = liqfrac_surface.create(grid, "liqfrac_surface", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = liqfrac_surface.set_attrs("climate_from_SurfaceModel",
                                   "liquid water fraction at the top surface of the ice",
                                   "1", ""); CHKERRQ(ierr);
  // ierr = variables.add(liqfrac_surface); CHKERRQ(ierr);

  // ice mass balance rate at the base of the ice shelf; sign convention for
  //   vshelfbasemass matches standard sign convention for basal melt rate of
  //   grounded ice
  ierr = shelfbmassflux.create(grid, "shelfbmassflux", WITHOUT_GHOSTS); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = shelfbmassflux.set_attrs(
                                  "climate_state", "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                                  "m s-1", ""); CHKERRQ(ierr);
  // PROPOSED standard name = ice_shelf_basal_specific_mass_balance
  // rescales from m/s to m/year when writing to NetCDF and std out:
  shelfbmassflux.write_in_glaciological_units = true;
  ierr = shelfbmassflux.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  // do not add; boundary models are in charge here
  //ierr = variables.add(shelfbmassflux); CHKERRQ(ierr);

  // ice boundary tempature at the base of the ice shelf
  ierr = shelfbtemp.create(grid, "shelfbtemp", WITHOUT_GHOSTS); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = shelfbtemp.set_attrs(
                              "climate_state", "absolute temperature at ice shelf base",
                              "K", ""); CHKERRQ(ierr);
  // PROPOSED standard name = ice_shelf_basal_temperature
  // do not add; boundary models are in charge here
  // ierr = variables.add(shelfbtemp); CHKERRQ(ierr);

  // take care of 2D cumulative fluxes: we need to allocate special storage if
  // the user asked for climatic_mass_balance_cumulative or some others (below).

  std::string extra_vars_argument;
  bool extra_vars_set = false;
  ierr = OptionsString("-extra_vars", "", extra_vars_argument, extra_vars_set); CHKERRQ(ierr);
  std::set<std::string> ex_vars;
  if (extra_vars_set == true) {
    std::istringstream arg(extra_vars_argument);
    std::string var_name;

    while (getline(arg, var_name, ',')) {
      ex_vars.insert(var_name);
    }
  }

  if (set_contains(ex_vars, "climatic_mass_balance_cumulative")) {
    ierr = climatic_mass_balance_cumulative.create(grid,
                                                   "climatic_mass_balance_cumulative",
                                                   WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = climatic_mass_balance_cumulative.set_attrs("diagnostic",
                                                      "cumulative surface mass balance",
                                                      "kg m-2", ""); CHKERRQ(ierr);
  }

  std::string o_size = get_output_size("-o_size");

  if (set_contains(ex_vars, "flux_divergence") || o_size == "big") {
    ierr = flux_divergence.create(grid, "flux_divergence", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = flux_divergence.set_attrs("diagnostic",
                                     "flux divergence",
                                     "m s-1", ""); CHKERRQ(ierr);
    ierr = flux_divergence.set_glaciological_units("m year-1"); CHKERRQ(ierr);
    flux_divergence.write_in_glaciological_units = true;
  }

  if (set_contains(ex_vars, "grounded_basal_flux_cumulative")) {
    ierr = grounded_basal_flux_2D_cumulative.create(grid, "grounded_basal_flux_cumulative", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = grounded_basal_flux_2D_cumulative.set_attrs("diagnostic",
                                                       "cumulative basal flux into the ice "
                                                       "in grounded areas (positive means ice gain)",
                                                       "kg m-2", ""); CHKERRQ(ierr);
    grounded_basal_flux_2D_cumulative.set_time_independent(false);
    ierr = grounded_basal_flux_2D_cumulative.set_glaciological_units("Gt m-2"); CHKERRQ(ierr);
    grounded_basal_flux_2D_cumulative.write_in_glaciological_units = true;
  }

  if (set_contains(ex_vars, "floating_basal_flux_cumulative")) {
    ierr = floating_basal_flux_2D_cumulative.create(grid, "floating_basal_flux_cumulative", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = floating_basal_flux_2D_cumulative.set_attrs("diagnostic",
                                                       "cumulative basal flux into the ice "
                                                       "in floating areas (positive means ice gain)",
                                                       "kg m-2", ""); CHKERRQ(ierr);
    floating_basal_flux_2D_cumulative.set_time_independent(false);
    ierr = floating_basal_flux_2D_cumulative.set_glaciological_units("Gt m-2"); CHKERRQ(ierr);
    floating_basal_flux_2D_cumulative.write_in_glaciological_units = true;
  }

  if (set_contains(ex_vars, "nonneg_flux_cumulative")) {
    ierr = nonneg_flux_2D_cumulative.create(grid, "nonneg_flux_cumulative", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = nonneg_flux_2D_cumulative.set_attrs("diagnostic",
                                               "cumulative nonnegative rule flux (positive means ice gain)",
                                               "kg m-2", ""); CHKERRQ(ierr);
    nonneg_flux_2D_cumulative.set_time_independent(false);
    ierr = nonneg_flux_2D_cumulative.set_glaciological_units("Gt m-2"); CHKERRQ(ierr);
    nonneg_flux_2D_cumulative.write_in_glaciological_units = true;
  }

  if (set_contains(ex_vars, "discharge_flux_cumulative")) {
    ierr = discharge_flux_2D_cumulative.create(grid, "discharge_flux_cumulative", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = discharge_flux_2D_cumulative.set_attrs("diagnostic",
                                                  "cumulative discharge (calving) flux (positive means ice loss)",
                                                  "kg m-2", ""); CHKERRQ(ierr);
    discharge_flux_2D_cumulative.set_time_independent(false);
    ierr = discharge_flux_2D_cumulative.set_glaciological_units("Gt m-2"); CHKERRQ(ierr);
    discharge_flux_2D_cumulative.write_in_glaciological_units = true;
  }

  return 0;
}

PetscErrorCode IceModel::setExecName(const std::string &my_executable_short_name) {
  executable_short_name = my_executable_short_name;
  return 0;
}

//! The contents of the main PISM time-step.
/*!
During the time-step we perform the following actions:
 */
PetscErrorCode IceModel::step(bool do_mass_continuity,
                              bool do_energy,
                              bool do_age,
                              bool do_skip) {
  PetscErrorCode ierr;

  //! \li call additionalAtStartTimestep() to let derived classes do more
  ierr = additionalAtStartTimestep(); CHKERRQ(ierr);  // might set dt_force,maxdt_temporary

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
    grid.profiling.begin("basal yield stress");
    ierr = basal_yield_stress_model->update(grid.time->current(), dt); CHKERRQ(ierr);
    grid.profiling.end("basal yield stress");
    ierr = basal_yield_stress_model->basal_material_yield_stress(basal_yield_stress); CHKERRQ(ierr);
    stdout_flags += "y";
  } else stdout_flags += "$";

  // Update the fractional grounded/floating mask (used by the SSA
  // stress balance and the energy code)
  if (config.get_flag("sub_groundingline")) {
    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr); // update h and mask
    ierr = update_floatation_mask(); CHKERRQ(ierr);
  }

  double sea_level = 0;
  ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);

  IceModelVec2S &melange_back_pressure = vWork2d[0];

  ierr = ocean->melange_back_pressure_fraction(melange_back_pressure); CHKERRQ(ierr);

  grid.profiling.begin("stress balance");
  ierr = stress_balance->update(updateAtDepth == false,
                                sea_level,
                                melange_back_pressure);
  grid.profiling.end("stress balance");
  if (ierr != 0) {
    std::string o_file = "stressbalance_failed.nc";
    bool o_file_set;
    ierr = OptionsString("-o", "output file name",
                             o_file, o_file_set); CHKERRQ(ierr);

    o_file = pism_filename_add_suffix(o_file, "_stressbalance_failed", "");
    ierr = PetscPrintf(grid.com,
                       "PISM ERROR: stress balance computation failed. Saving model state to '%s'...\n",
                       o_file.c_str());

    ierr = dumpToFile(o_file); CHKERRQ(ierr);

    ierr = PetscPrintf(grid.com, "ending...\n");
    PISMEnd();
  }

  std::string sb_stdout;
  ierr = stress_balance->stdout_report(sb_stdout); CHKERRQ(ierr);

  stdout_flags += sb_stdout;

  stdout_flags += (updateAtDepth ? "v" : "V");

  //! \li determine the time step according to a variety of stability criteria;
  //!  see determineTimeStep()
  ierr = max_timestep(dt, skipCountDown); CHKERRQ(ierr);

  //! \li Update surface and ocean models.
  grid.profiling.begin("surface");
  ierr = surface->update(grid.time->current(), dt); CHKERRQ(ierr);
  grid.profiling.end("surface");

  grid.profiling.begin("ocean");
  ierr = ocean->update(grid.time->current(),   dt); CHKERRQ(ierr);
  grid.profiling.end("ocean");

  dt_TempAge += dt;
  // IceModel::dt,dtTempAge are now set correctly according to
  // mass-continuity-eqn-diffusivity criteria, horizontal CFL criteria, and
  // other criteria from derived class additionalAtStartTimestep(), and from
  // "-skip" mechanism

  //! \li update the age of the ice (if appropriate)
  if (do_age && updateAtDepth) {
    grid.profiling.begin("age");
    ierr = ageStep(); CHKERRQ(ierr);
    grid.profiling.end("age");
    stdout_flags += "a";
  } else {
    stdout_flags += "$";
  }

  //! \li update the enthalpy (or temperature) field according to the conservation of
  //!  energy model based (especially) on the new velocity field; see
  //!  energyStep()
  if (do_energy_step) { // do the energy step
    grid.profiling.begin("energy");
    ierr = energyStep(); CHKERRQ(ierr);
    grid.profiling.end("energy");
    stdout_flags += "E";
  } else {
    stdout_flags += "$";
  }

  // Combine basal melt rate in grounded (computed during the energy
  // step) and floating (provided by an ocean model) areas.
  ierr = combine_basal_melt_rate(); CHKERRQ(ierr);

  //! \li update the state variables in the subglacial hydrology model (typically
  //!  water thickness and sometimes pressure)
  grid.profiling.begin("basal hydrology");
  ierr = subglacial_hydrology->update(grid.time->current(), dt); CHKERRQ(ierr);
  grid.profiling.end("basal hydrology");

  //! \li update the fracture density field; see calculateFractureDensity()
  if (config.get_flag("do_fracture_density")) {
    grid.profiling.begin("fracture density");
    ierr = calculateFractureDensity(); CHKERRQ(ierr);
    grid.profiling.end("fracture density");
  }

  //! \li update the thickness of the ice according to the mass conservation
  //!  model; see massContExplicitStep()
  if (do_mass_continuity) {
    grid.profiling.begin("mass transport");
    ierr = massContExplicitStep(); CHKERRQ(ierr);
    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr); // update h and mask
    grid.profiling.end("mass transport");

    // Note that there are three adaptive time-stepping criteria. Two of them
    // (using max. diffusion and 2D CFL) are limiting the mass-continuity
    // time-step and the third (3D CFL) limits the energy and age time-steps.

    // The mass-continuity time-step is usually smaller, and the skipping
    // mechanism lets us do several mass-continuity steps for each energy step.

    // When -no_mass is set, mass-continuity-related time-step restrictions are
    // disabled, making "skipping" unnecessary.

    // This is why the following two lines appear here and are executed only
    // if do_mass_continuity == true.
    if (do_skip == true && skipCountDown > 0)
      skipCountDown--;
    stdout_flags += "h";
  } else {
    stdout_flags += "$";
  }

  grid.profiling.begin("calving");
  ierr = do_calving(); CHKERRQ(ierr);
  grid.profiling.end("calving");

  ierr = Href_cleanup(); CHKERRQ(ierr);

  //! \li compute the bed deformation, which only depends on current thickness
  //! and bed elevation
  if (beddef) {
    int topg_state_counter = bed_topography.get_state_counter();

    grid.profiling.begin("bed deformation");
    ierr = beddef->update(grid.time->current(), dt); CHKERRQ(ierr);
    grid.profiling.end("bed deformation");

    if (bed_topography.get_state_counter() != topg_state_counter) {
      stdout_flags += "b";
      ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);
    } else
      stdout_flags += " ";
  }

  //! \li call additionalAtEndTimestep() to let derived classes do more
  ierr = additionalAtEndTimestep(); CHKERRQ(ierr);

  // Done with the step; now adopt the new time.
  grid.time->step(dt);

  if (do_energy_step) {
    t_TempAge = grid.time->current();
    dt_TempAge = 0.0;
  }

  // end the flag line
  stdout_flags += " " + m_adaptive_timestep_reason;

#if (PISM_DEBUG==1)
  ierr = variables.check_for_nan(); CHKERRQ(ierr);
#endif

  return 0;
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
PetscErrorCode IceModel::run_to(double run_end) {
  PetscErrorCode  ierr;

  grid.time->set_end(run_end);

  ierr = run(); CHKERRQ(ierr);

  return 0;
}


/**
 * Run the time-stepping loop from the current time until the time
 * specified by the IceModel::grid::time object.
 *
 * This is the method used by PISM in the "standalone" mode.
 *
 * @return 0 on success
 */
PetscErrorCode IceModel::run() {
  PetscErrorCode  ierr;

  bool do_mass_conserve = config.get_flag("do_mass_conserve");
  bool do_energy = config.get_flag("do_energy");
  bool do_age = config.get_flag("do_age");
  bool do_skip = config.get_flag("do_skip");

  int stepcount = config.get_flag("count_time_steps") ? 0 : -1;

  ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);

  // update diagnostics at the beginning of the run:
  ierr = write_timeseries(); CHKERRQ(ierr);
  ierr = write_extras(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "running forward ...\n"); CHKERRQ(ierr);

  stdout_flags.erase(); // clear it out
  ierr = summaryPrintLine(PETSC_TRUE, do_energy, 0.0, 0.0, 0.0, 0.0, 0.0); CHKERRQ(ierr);
  m_adaptive_timestep_reason = '$'; // no reason for no timestep
  ierr = summary(do_energy); CHKERRQ(ierr);  // report starting state

  t_TempAge = grid.time->current();
  dt_TempAge = 0.0;

  // main loop for time evolution
  // IceModel::step calls grid.time->step(dt), ensuring that this while loop
  // will terminate
  grid.profiling.stage_begin("time-stepping loop");
  while (grid.time->current() < grid.time->end()) {

    stdout_flags.erase();  // clear it out
    dt_force = -1.0;
    maxdt_temporary = -1.0;

    ierr = step(do_mass_conserve, do_energy, do_age, do_skip); CHKERRQ(ierr);

    // report a summary for major steps or the last one
    bool updateAtDepth = skipCountDown == 0;
    bool tempAgeStep = updateAtDepth && (do_energy || do_age);

    const bool show_step = tempAgeStep || m_adaptive_timestep_reason == "end of the run";
    ierr = summary(show_step); CHKERRQ(ierr);

    // writing these fields here ensures that we do it after the last time-step
    grid.profiling.begin("I/O during run");
    ierr = write_snapshot(); CHKERRQ(ierr);
    ierr = write_timeseries(); CHKERRQ(ierr);
    ierr = write_extras(); CHKERRQ(ierr);
    ierr = write_backup(); CHKERRQ(ierr);
    grid.profiling.end("I/O during run");

    ierr = update_viewers(); CHKERRQ(ierr);

    if (stepcount >= 0) stepcount++;
    if (endOfTimeStepHook() != 0) break;
  } // end of the time-stepping loop

  grid.profiling.stage_end("time-stepping loop");

  bool flag;
  int pause_time = 0;
  ierr = OptionsInt("-pause", "Pause after the run, seconds",
                        pause_time, flag); CHKERRQ(ierr);
  if (pause_time > 0) {
    ierr = verbPrintf(2,grid.com,"pausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);
  }

  if (stepcount >= 0) {
    ierr = verbPrintf(1,grid.com,
                      "count_time_steps:  run() took %d steps\n"
                      "average dt = %.6f years\n",
                      stepcount,
                      grid.convert(grid.time->end() - grid.time->start(), "seconds", "years")/(double)stepcount); CHKERRQ(ierr);
  }

  return 0;
}

//! Manage the initialization of the IceModel object.
/*!
Please see the documenting comments of the functions called below to find
explanations of their intended uses.
 */
PetscErrorCode IceModel::init() {
  PetscErrorCode ierr;

  ierr = PetscOptionsBegin(grid.com, "", "PISM options", ""); CHKERRQ(ierr);

  // Build PISM with -DPISM_WAIT_FOR_GDB=1 and run with -wait_for_gdb to
  // make it wait for a connection.
#ifdef PISM_WAIT_FOR_GDB
  bool wait_for_gdb = false;
  ierr = OptionsIsSet("-wait_for_gdb", wait_for_gdb); CHKERRQ(ierr);
  if (wait_for_gdb) {
    ierr = pism_wait_for_gdb(grid.com, 0); CHKERRQ(ierr);
  }
#endif
  //! The IceModel initialization sequence is this:

  //! 1) Initialize the computational grid:
  ierr = grid_setup(); CHKERRQ(ierr);

  //! 2) Process the options:
  ierr = setFromOptions(); CHKERRQ(ierr);

  //! 3) Memory allocation:
  ierr = createVecs(); CHKERRQ(ierr);

  //! 4) Allocate PISM components modeling some physical processes.
  ierr = allocate_submodels(); CHKERRQ(ierr);

  //! 5) Allocate work vectors:
  ierr = allocate_internal_objects(); CHKERRQ(ierr);

  //! 6) Initialize coupler models and fill the model state variables
  //! (from a PISM output file, from a bootstrapping file using some
  //! modeling choices or using formulas). Calls IceModel::regrid()
  ierr = model_state_setup(); CHKERRQ(ierr);

  //! 7) Report grid parameters:
  ierr = grid.report_parameters(); CHKERRQ(ierr);

  //! 8) Miscellaneous stuff: set up the bed deformation model, initialize the
  //! basal till model, initialize snapshots. This has to happen *after*
  //! regridding.
  ierr = misc_setup();

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  //! The following flow-chart illustrates the process.
  //!
  //! \dotfile initialization-sequence.dot "IceModel initialization sequence"

  // Get the start time in seconds and ensure that it is consistent
  // across all processors.
  {
    MPI_Datatype mpi_type;
    double my_start_time;
    ierr = PetscDataTypeToMPIDataType(PETSC_DOUBLE, &mpi_type); CHKERRQ(ierr);

    ierr = GetTime(&my_start_time); CHKERRQ(ierr);
    MPI_Allreduce(&my_start_time, &start_time, 1, mpi_type, MPI_MAX, grid.com);

  }
  return 0;
}

// FIXME: THIS IS BAD! (Provides unguarded access to IceModel's internals.)
IceModelVec2S* IceModel::get_geothermal_flux() {
  return &this->geothermal_flux;
}

// FIXME: THIS IS BAD! (Provides unguarded access to IceModel's internals.)
StressBalance* IceModel::get_stress_balance() {
  return this->stress_balance;
}

} // end of namespace pism

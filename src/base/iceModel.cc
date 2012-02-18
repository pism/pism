// Copyright (C) 2004-2012 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include <cstring>
#include <petscdmda.h>

#include "iceModel.hh"
#include "pism_signal.h"
#include "PISMStressBalance.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "PISMBedDef.hh"
#include "bedrockThermalUnit.hh"
#include "PISMYieldStress.hh"
#include "basal_resistance.hh"
#include "enthalpyConverter.hh"
#include "PISMProf.hh"
#include "pism_options.hh"


IceModel::IceModel(IceGrid &g, NCConfigVariable &conf, NCConfigVariable &conf_overrides)
  : grid(g), config(conf), overrides(conf_overrides), ice(NULL) {

  if (utIsInit() == 0) {
    if (utInit(NULL) != 0) {
      PetscPrintf(grid.com, "PISM ERROR: UDUNITS initialization failed.\n");
      PISMEnd();
    }
  }

  mapping.init("mapping", grid.com, grid.rank);
  global_attributes.init("global_attributes", grid.com, grid.rank);

  pism_signal = 0;
  signal(SIGTERM, pism_signal_handler);
  signal(SIGUSR1, pism_signal_handler);
  signal(SIGUSR2, pism_signal_handler);

  basal_yield_stress = NULL;
  basal = NULL;

  stress_balance = NULL;

  surface = NULL;
  ocean   = NULL;
  beddef  = NULL;

  EC = NULL;
  btu = NULL;

  executable_short_name = "pism"; // drivers typically override this

  shelvesDragToo = PETSC_FALSE;

  // initializr maximum |u|,|v|,|w| in ice
  gmaxu = gmaxv = gmaxw = 0;

  // set default locations of soundings and slices
  id = (grid.Mx - 1)/2;
  jd = (grid.My - 1)/2;

  // frequently used physical constants and parameters:
  standard_gravity = config.get("standard_gravity");

  global_attributes.set_string("Conventions", "CF-1.4");
  global_attributes.set_string("source", string("PISM ") + PISM_Revision);

  // Do not save snapshots by default:
  save_snapshots = false;
  // Do not save time-series by default:
  save_ts = false;
  save_extra = false;

  reset_counters();

  allowAboveMelting = PETSC_FALSE;  // only IceCompModel ever sets it to true
}

void IceModel::reset_counters() {
  CFLmaxdt = CFLmaxdt2D = 0.0;
  CFLviolcount = 0;
  dt_TempAge = 0.0;
  dt_from_diffus = dt_from_cfl = 0.0;
  gmaxu = gmaxv = gmaxw = 0;
  maxdt_temporary = dt = dt_force = 0.0;
  skipCountDown = 0;

  cumulative_basal_ice_flux = 0;
  cumulative_float_kill_flux = 0;
  cumulative_discharge_flux = 0;
  cumulative_nonneg_rule_flux = 0;
  cumulative_ocean_kill_flux = 0;
  cumulative_sub_shelf_ice_flux = 0;
  cumulative_surface_ice_flux = 0;
}


IceModel::~IceModel() {

  deallocate_internal_objects();

  // de-allocate time-series diagnostics
  map<string,PISMTSDiagnostic*>::iterator i = ts_diagnostics.begin();
  while (i != ts_diagnostics.end()) delete (i++)->second;

  // de-allocate diagnostics
  map<string,PISMDiagnostic*>::iterator j = diagnostics.begin();
  while (j != diagnostics.end()) delete (j++)->second;


  // de-allocate viewers
  map<string,PetscViewer>::iterator k = viewers.begin();
  while (k != viewers.end()) {
    if ((*k).second != PETSC_NULL) {
      PetscViewerDestroy(&(*k).second);
      ++k;
    }
  }

  delete stress_balance;

  delete ocean;
  delete surface;
  delete beddef;

  delete basal;
  delete EC;
  delete btu;
  delete ice;

  utTerm(); // Clean up after UDUNITS
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
  PetscInt WIDE_STENCIL = grid.max_stencil_width;

  ierr = verbPrintf(3, grid.com,
		    "Allocating memory...\n"); CHKERRQ(ierr);

  // The following code creates (and documents -- to some extent) the
  // variables. The main (and only) principle here is using standard names from
  // the CF conventions; see
  // http://cf-pcmdi.llnl.gov/documents/cf-standard-names

  ierr = Enth3.create(grid, "enthalpy", true, WIDE_STENCIL); CHKERRQ(ierr);
  // POSSIBLE standard name = land_ice_enthalpy
  ierr = Enth3.set_attrs(
                         "model_state",
                         "ice enthalpy (includes sensible heat, latent heat, pressure)",
                         "J kg-1", ""); CHKERRQ(ierr);
  ierr = variables.add(Enth3); CHKERRQ(ierr);

  if (config.get_flag("do_cold_ice_methods")) {
    // ice temperature
    ierr = T3.create(grid, "temp", true); CHKERRQ(ierr);
    ierr = T3.set_attrs("model_state", "ice temperature", "K", "land_ice_temperature"); CHKERRQ(ierr);
    ierr = T3.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    ierr = variables.add(T3); CHKERRQ(ierr);

    ierr = Enth3.set_attr("pism_intent", "diagnostic"); CHKERRQ(ierr); 
  }

  // age of ice but only if age will be computed
  if (config.get_flag("do_age")) {
    ierr = tau3.create(grid, "age", true, WIDE_STENCIL); CHKERRQ(ierr);
    // PROPOSED standard_name = land_ice_age
    ierr = tau3.set_attrs("model_state", "age of ice",
                          "s", ""); CHKERRQ(ierr);
    ierr = tau3.set_glaciological_units("years");
    tau3.write_in_glaciological_units = true;
    ierr = tau3.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    ierr = variables.add(tau3); CHKERRQ(ierr);
  }

  // ice upper surface elevation
  ierr = vh.create(grid, "usurf", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = vh.set_attrs("diagnostic", "ice upper surface elevation",
		      "m", "surface_altitude"); CHKERRQ(ierr);
  ierr = variables.add(vh); CHKERRQ(ierr);

  // land ice thickness
  ierr = vH.create(grid, "thk", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = vH.set_attrs("model_state", "land ice thickness",
		      "m", "land_ice_thickness"); CHKERRQ(ierr);
  ierr = vH.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = variables.add(vH); CHKERRQ(ierr);

  // bedrock surface elevation
  ierr = vbed.create(grid, "topg", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = vbed.set_attrs("model_state", "bedrock surface elevation",
			"m", "bedrock_altitude"); CHKERRQ(ierr);
  ierr = variables.add(vbed); CHKERRQ(ierr);

  if (config.get_flag("ocean_kill")) {
    ierr = ocean_kill_mask.create(grid, "ocean_kill_mask", false); CHKERRQ(ierr);
    ierr = ocean_kill_mask.set_attrs("internal",
                                     "mask specifying fixed calving front locations",
                                     "", ""); CHKERRQ(ierr);
    ocean_kill_mask.time_independent = true;
    ierr = variables.add(ocean_kill_mask); CHKERRQ(ierr);

  }

  // grounded_dragging_floating integer mask
  if(config.get_flag("do_eigen_calving")) {
    ierr = vMask.create(grid, "mask", true, 3); CHKERRQ(ierr); 
    // The wider stencil is needed for parallel calculation in iMcalving.cc when asking for mask values at the front (offset+1)
  } else {
    ierr = vMask.create(grid, "mask", true, WIDE_STENCIL); CHKERRQ(ierr);
  }
  ierr = vMask.set_attrs("diagnostic", "grounded_dragging_floating integer mask",
			 "", ""); CHKERRQ(ierr);
  vector<double> mask_values(4);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_GROUNDED;
  mask_values[2] = MASK_FLOATING;
  mask_values[3] = MASK_ICE_FREE_OCEAN;
  ierr = vMask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
  ierr = vMask.set_attr("flag_meanings",
			"ice_free_bedrock grounded_ice floating_ice ice_free_ocean"); CHKERRQ(ierr);
  vMask.output_data_type = NC_BYTE;
  ierr = variables.add(vMask); CHKERRQ(ierr);

  // iceberg identifying integer mask
  if (config.get_flag("kill_icebergs")) {
    ierr = vIcebergMask.create(grid, "IcebergMask", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vIcebergMask.set_attrs("internal", 
                                  "iceberg-identifying integer mask",
                                  "", ""); CHKERRQ(ierr);
    vector<double> icebergmask_values(5);
    icebergmask_values[0] = ICEBERGMASK_NO_ICEBERG;
    icebergmask_values[1] = ICEBERGMASK_NOT_SET;
    icebergmask_values[2] = ICEBERGMASK_ICEBERG_CAND;
    icebergmask_values[3] = ICEBERGMASK_STOP_OCEAN;
    icebergmask_values[4] = ICEBERGMASK_STOP_ATTACHED;
    //more values to identify lakes?

    ierr = vIcebergMask.set_attr("flag_values", icebergmask_values); CHKERRQ(ierr);
    ierr = vIcebergMask.set_attr("flag_meanings",
                                 "no_iceberg not_set iceberg_candidate ocean_boundary grounded_boundary"); CHKERRQ(ierr);
    vIcebergMask.output_data_type = NC_BYTE;
    ierr = variables.add(vIcebergMask); CHKERRQ(ierr);
  }

  // upward geothermal flux at bedrock surface
  ierr = vGhf.create(grid, "bheatflx", true, WIDE_STENCIL); CHKERRQ(ierr); // never differentiated
  // PROPOSED standard_name = lithosphere_upward_heat_flux
  ierr = vGhf.set_attrs("climate_steady", "upward geothermal flux at bedrock surface",
			"W m-2", ""); CHKERRQ(ierr);
  ierr = vGhf.set_glaciological_units("mW m-2");
  vGhf.write_in_glaciological_units = true;
  vGhf.time_independent = true;
  ierr = variables.add(vGhf); CHKERRQ(ierr);

  // temperature seen by top of bedrock thermal layer
  ierr = bedtoptemp.create(grid, "bedtoptemp", false); CHKERRQ(ierr); // never differentiated
  ierr = bedtoptemp.set_attrs("internal",
                              "temperature of top of bedrock thermal layer",
                              "K", ""); CHKERRQ(ierr);
  ierr = bedtoptemp.set_glaciological_units("K");
  ierr = variables.add(bedtoptemp); CHKERRQ(ierr);

  // effective thickness of subglacial melt water
  ierr = vbwat.create(grid, "bwat", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = vbwat.set_attrs("model_state", "effective thickness of subglacial melt water",
                         "m", ""); CHKERRQ(ierr);
  ierr = vbwat.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = vbwat.set_attr("valid_max", config.get("bwat_max")); CHKERRQ(ierr);
  ierr = variables.add(vbwat); CHKERRQ(ierr);

  if (config.get_flag("use_ssa_velocity") || config.get_flag("do_blatter")) {
    // yield stress for basal till (plastic or pseudo-plastic model)
    ierr = vtauc.create(grid, "tauc", true, WIDE_STENCIL); CHKERRQ(ierr);
    // PROPOSED standard_name = land_ice_basal_material_yield_stress
    ierr = vtauc.set_attrs("diagnostic", 
                           "yield stress for basal till (plastic or pseudo-plastic model)",
                           "Pa", ""); CHKERRQ(ierr);
    ierr = variables.add(vtauc); CHKERRQ(ierr);
  }

  // bedrock uplift rate
  ierr = vuplift.create(grid, "dbdt", false); CHKERRQ(ierr);
  ierr = vuplift.set_attrs("model_state", "bedrock uplift rate",
			   "m s-1", "tendency_of_bedrock_altitude"); CHKERRQ(ierr);
  ierr = vuplift.set_glaciological_units("m year-1");
  vuplift.write_in_glaciological_units = true;
  ierr = variables.add(vuplift); CHKERRQ(ierr);

  // basal melt rate
  ierr = vbmr.create(grid, "bmelt", true, WIDE_STENCIL); CHKERRQ(ierr);
  // ghosted to allow the "redundant" computation of tauc
  ierr = vbmr.set_attrs("model_state",
                        "ice basal melt rate in ice thickness per time",
                        "m s-1", "land_ice_basal_melt_rate"); CHKERRQ(ierr);
  ierr = vbmr.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vbmr.write_in_glaciological_units = true;
  vbmr.set_attr("comment", "positive basal melt rate corresponds to ice loss");
  ierr = variables.add(vbmr); CHKERRQ(ierr);

  // longitude
  ierr = vLongitude.create(grid, "lon", true); CHKERRQ(ierr);
  ierr = vLongitude.set_attrs("mapping", "longitude", "degree_east", "longitude"); CHKERRQ(ierr);
  vLongitude.time_independent = true;
  ierr = vLongitude.set_attr("coordinates", ""); CHKERRQ(ierr);
  ierr = vLongitude.set_attr("grid_mapping", ""); CHKERRQ(ierr);
  ierr = vLongitude.set_attr("valid_min", -180.0); CHKERRQ(ierr);
  ierr = vLongitude.set_attr("valid_max",  180.0); CHKERRQ(ierr);
  ierr = variables.add(vLongitude); CHKERRQ(ierr);

  // latitude
  ierr = vLatitude.create(grid, "lat", true); CHKERRQ(ierr); // has ghosts so that we can compute cell areas
  ierr = vLatitude.set_attrs("mapping", "latitude", "degree_north", "latitude"); CHKERRQ(ierr);
  vLatitude.time_independent = true;
  ierr = vLatitude.set_attr("coordinates", ""); CHKERRQ(ierr);
  ierr = vLatitude.set_attr("grid_mapping", ""); CHKERRQ(ierr);
  ierr = vLatitude.set_attr("valid_min", -90.0); CHKERRQ(ierr);
  ierr = vLatitude.set_attr("valid_max",  90.0); CHKERRQ(ierr);
  ierr = variables.add(vLatitude); CHKERRQ(ierr);

  if (config.get_flag("part_grid") == true) {
    // Href
    ierr = vHref.create(grid, "Href", true); CHKERRQ(ierr);
    ierr = vHref.set_attrs("model_state", "temporary ice thickness at calving front boundary",
                           "m", ""); CHKERRQ(ierr);
    ierr = variables.add(vHref); CHKERRQ(ierr);

    if (config.get_flag("part_redist") == true){
      // Hav
      ierr = vHresidual.create(grid, "Hresidual", true); CHKERRQ(ierr);
      ierr = vHresidual.set_attrs("diagnostic", "residual ice thickness in recently filled boundary grid cell",
                                  "m", ""); CHKERRQ(ierr);
      ierr = variables.add(vHresidual); CHKERRQ(ierr);
    }
  }

  if (config.get_flag("do_eigen_calving") == true) {
    ierr = vPrinStrain1.create(grid, "edot_1", true); CHKERRQ(ierr);
    ierr = vPrinStrain1.set_attrs("internal", 
                                  "major principal component of horizontal strain-rate",
                                  "1/s", ""); CHKERRQ(ierr);
    ierr = variables.add(vPrinStrain1); CHKERRQ(ierr);
    ierr = vPrinStrain2.create(grid, "edot_2", true); CHKERRQ(ierr);
    ierr = vPrinStrain2.set_attrs("internal",
                                  "minor principal component of horizontal strain-rate",
                                  "1/s", ""); CHKERRQ(ierr);
    ierr = variables.add(vPrinStrain2); CHKERRQ(ierr);
  }

  if (config.get_flag("ssa_dirichlet_bc") == true) {
    // bc_locations
    ierr = vBCMask.create(grid, "bcflag", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vBCMask.set_attrs("model_state", "Dirichlet boundary mask",
                             "", ""); CHKERRQ(ierr);
    vector<double> bc_mask_values(2);
    bc_mask_values[0] = 0;
    bc_mask_values[1] = 1;
    ierr = vBCMask.set_attr("flag_values", bc_mask_values); CHKERRQ(ierr);
    ierr = vBCMask.set_attr("flag_meanings", "no_data bc_condition"); CHKERRQ(ierr);
    vBCMask.output_data_type = NC_BYTE;
    ierr = variables.add(vBCMask); CHKERRQ(ierr);


    // vel_bc
    ierr = vBCvel.create(grid, "_ssa_bc", true, WIDE_STENCIL); CHKERRQ(ierr); // u_ssa_bc and v_ssa_bc
    ierr = vBCvel.set_attrs("model_state",
                            "X-component of the SSA velocity boundary conditions",
                            "m s-1", "", 0); CHKERRQ(ierr);
    ierr = vBCvel.set_attrs("model_state",
                            "Y-component of the SSA velocity boundary conditions",
                            "m s-1", "", 1); CHKERRQ(ierr);
    ierr = vBCvel.set_glaciological_units("m year-1"); CHKERRQ(ierr);
    for (int j = 0; j < 2; ++j) {
      ierr = vBCvel.set_attr("valid_min",  convert(-1e6, "m/year", "m/second"), j); CHKERRQ(ierr);
      ierr = vBCvel.set_attr("valid_max",  convert(1e6, "m/year", "m/second"), j); CHKERRQ(ierr);
      ierr = vBCvel.set_attr("_FillValue", convert(2e6, "m/year", "m/second"), j); CHKERRQ(ierr);
    }
    //just for diagnostics...
    ierr = variables.add(vBCvel); CHKERRQ(ierr);
  }


  // cell areas
  ierr = cell_area.create(grid, "cell_area", false); CHKERRQ(ierr);
  ierr = cell_area.set_attrs("diagnostic", "cell areas", "m2", ""); CHKERRQ(ierr);
  ierr = cell_area.set_attr("comment",
                            "values are equal to dx*dy "
                            "if projection parameters are not available; "
                            "otherwise WGS84 ellipsoid is used"); CHKERRQ(ierr); 
  cell_area.time_independent = true;
  ierr = cell_area.set_glaciological_units("km2"); CHKERRQ(ierr);
  cell_area.write_in_glaciological_units = true;
  ierr = variables.add(cell_area); CHKERRQ(ierr);

  // fields owned by IceModel but filled by PISMSurfaceModel *surface:
  // mean annual net ice equivalent surface mass balance rate
  ierr = acab.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = acab.set_attrs(
                        "climate_from_PISMSurfaceModel",  // FIXME: can we do better?
                        "ice-equivalent surface mass balance (accumulation/ablation) rate",
                        "m s-1",  // m *ice-equivalent* per second
                        "land_ice_surface_specific_mass_balance");  // CF standard_name
  CHKERRQ(ierr);
  ierr = acab.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  acab.write_in_glaciological_units = true;
  acab.set_attr("comment", "positive values correspond to ice gain");

  if (config.get_flag("compute_cumulative_acab")) {
    ierr = acab_cumulative.create(grid, "acab_cumulative", false); CHKERRQ(ierr);
    ierr = acab_cumulative.set_attrs("diagnostic",
                                     "cumulative ice-equivalent surface mass balance",
                                     "m", ""); CHKERRQ(ierr);
    ierr = variables.add(acab_cumulative); CHKERRQ(ierr);
  }

  // annual mean air temperature at "ice surface", at level below all firn
  //   processes (e.g. "10 m" or ice temperatures)
  ierr = artm.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = artm.set_attrs(
                        "climate_from_PISMSurfaceModel",  // FIXME: can we do better?
                        "annual average ice surface temperature, below firn processes",
                        "K", 
                        "");  // PROPOSED CF standard_name = land_ice_surface_temperature_below_firn
  CHKERRQ(ierr);

  ierr = liqfrac_surface.create(grid, "liqfrac_surface", false); CHKERRQ(ierr);
  ierr = liqfrac_surface.set_attrs("climate_from_PISMSurfaceModel",
                                   "liquid water fraction at the top surface of the ice",
                                   "1", ""); CHKERRQ(ierr);
  // ierr = variables.add(liqfrac_surface); CHKERRQ(ierr);

  // ice mass balance rate at the base of the ice shelf; sign convention for
  //   vshelfbasemass matches standard sign convention for basal melt rate of
  //   grounded ice
  ierr = shelfbmassflux.create(grid, "shelfbmassflux", false); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = shelfbmassflux.set_attrs(
                                  "climate_state", "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                                  "m s-1", ""); CHKERRQ(ierr); 
  // PROPOSED standard name = ice_shelf_basal_specific_mass_balance
  // rescales from m/s to m/a when writing to NetCDF and std out:
  shelfbmassflux.write_in_glaciological_units = true;
  ierr = shelfbmassflux.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  // do not add; boundary models are in charge here
  //ierr = variables.add(shelfbmassflux); CHKERRQ(ierr);

  // ice boundary tempature at the base of the ice shelf
  ierr = shelfbtemp.create(grid, "shelfbtemp", false); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = shelfbtemp.set_attrs(
                              "climate_state", "absolute temperature at ice shelf base",
                              "K", ""); CHKERRQ(ierr);
  // PROPOSED standard name = ice_shelf_basal_temperature
  // do not add; boundary models are in charge here
  // ierr = variables.add(shelfbtemp); CHKERRQ(ierr);

  return 0;
}


//! De-allocate internal objects.
/*! This includes Vecs that are not in an IceModelVec, SSA tools and the bed
  deformation model.
 */
PetscErrorCode IceModel::deallocate_internal_objects() {
  return 0;
}

PetscErrorCode IceModel::setExecName(string my_executable_short_name) {
  executable_short_name = my_executable_short_name;
  return 0;
}

//! Do the contents of the main PISM time-step.
/*!
During the time-step we perform the following actions:
 */
PetscErrorCode IceModel::step(bool do_mass_continuity,
			      bool do_energy,
			      bool do_age,
			      bool do_skip) {
  PetscErrorCode ierr;

  grid.profiler->begin(event_step);

  //! \li call additionalAtStartTimestep() to let derived classes do more
  ierr = additionalAtStartTimestep(); CHKERRQ(ierr);  // might set dt_force,maxdt_temporary

  //! \li determine the maximum time-step boundary models can take
  double apcc_dt;
  bool restrict_dt;
  ierr = surface->max_timestep(grid.time->current(), apcc_dt, restrict_dt); CHKERRQ(ierr);
  if (restrict_dt) {
    if (maxdt_temporary > 0)
      maxdt_temporary = PetscMin(apcc_dt, maxdt_temporary);
    else
      maxdt_temporary = apcc_dt;
  }

  double opcc_dt;
  ierr = ocean->max_timestep(grid.time->current(), opcc_dt, restrict_dt); CHKERRQ(ierr);
  if (restrict_dt) {
    if (maxdt_temporary > 0)
      maxdt_temporary = PetscMin(opcc_dt, maxdt_temporary);
    else
      maxdt_temporary = opcc_dt;
  }

  double ts_dt;
  ierr = ts_max_timestep(grid.time->current(), ts_dt, restrict_dt); CHKERRQ(ierr);
  if (restrict_dt) {
    if (maxdt_temporary > 0)
      maxdt_temporary = PetscMin(ts_dt, maxdt_temporary);
    else
      maxdt_temporary = ts_dt;
  }

  //! \li apply the time-step restriction from the -extra_{times,file,vars}
  //! mechanism
  double extras_dt;
  ierr = extras_max_timestep(grid.time->current(), extras_dt, restrict_dt); CHKERRQ(ierr);
  if (restrict_dt) {
    if (maxdt_temporary > 0)
      maxdt_temporary = PetscMin(extras_dt, maxdt_temporary);
    else
      maxdt_temporary = extras_dt;
  }

  //! \li update the velocity field; in some cases the whole three-dimensional
  //! field is updated and in some cases just the vertically-averaged
  //! horizontal velocity is updated; see velocity()

  // always do SIA velocity calculation (unless -no_sia is given); only update
  // SSA and only update velocities at depth if suggested by temp and age
  // stability criterion; note *lots* of communication is avoided by skipping
  // SSA (and temp/age)

  bool updateAtDepth = (skipCountDown == 0),
    do_energy_step = updateAtDepth && do_energy;

  //! \li update the yield stress for the plastic till model (if appropriate)
  if (updateAtDepth && basal_yield_stress) {
    ierr = basal_yield_stress->update(grid.time->current(), dt); CHKERRQ(ierr);
    ierr = basal_yield_stress->basal_material_yield_stress(vtauc); CHKERRQ(ierr);
    stdout_flags += "y";
  } else stdout_flags += "$";

  grid.profiler->begin(event_velocity);

  ierr = stress_balance->update(updateAtDepth == false); CHKERRQ(ierr); 

  grid.profiler->end(event_velocity);

  string sb_stdout;
  ierr = stress_balance->stdout_report(sb_stdout); CHKERRQ(ierr);

  stdout_flags += sb_stdout;

  stdout_flags += (updateAtDepth ? "v" : "V");

  // communication here for global max; sets CFLmaxdt2D
  ierr = computeMax2DSlidingSpeed(); CHKERRQ(ierr);   

  if (updateAtDepth) {
    // communication here for global max; sets CFLmaxdt
    ierr = computeMax3DVelocities(); CHKERRQ(ierr); 
  }

  //! \li determine the time step according to a variety of stability criteria;
  //!  see determineTimeStep()
  ierr = determineTimeStep(do_energy); CHKERRQ(ierr);

  //! \li Update surface and ocean models.
  ierr = surface->update(grid.time->current(), dt); CHKERRQ(ierr);
  ierr = ocean->update(grid.time->current(),   dt); CHKERRQ(ierr);

  dt_TempAge += dt;
  // IceModel::dt,dtTempAge are now set correctly according to
  // mass-continuity-eqn-diffusivity criteria, horizontal CFL criteria, and
  // other criteria from derived class additionalAtStartTimestep(), and from
  // "-skip" mechanism

  grid.profiler->begin(event_age);

  //! \li update the age of the ice (if appropriate)
  if (do_age && updateAtDepth) {
    ierr = ageStep(); CHKERRQ(ierr);
    stdout_flags += "a";
  } else {
    stdout_flags += "$";
  }

  grid.profiler->end(event_age);

  grid.profiler->begin(event_energy);

  //! \li update the enthalpy (or temperature) field according to the conservation of
  //!  energy model based (especially) on the new velocity field; see
  //!  energyStep()
  if (do_energy_step) { // do the energy step
    ierr = energyStep(); CHKERRQ(ierr);
    stdout_flags += "E";
  } else {
    stdout_flags += "$";
  }

  grid.profiler->end(event_energy);

  // finally, diffuse the stored basal water once per energy step, if it is requested
  if (do_energy_step && config.get_flag("do_diffuse_bwat")) {
    ierr = diffuse_bwat(); CHKERRQ(ierr);
  }

  grid.profiler->begin(event_mass);

  //! \li update the thickness of the ice according to the mass conservation
  //!  model; see massContExplicitStep()
  if (do_mass_continuity) {
    ierr = massContExplicitStep(); CHKERRQ(ierr); // update H
    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr); // update h and mask

    // Note that there are three adaptive time-stepping criteria. Two of them
    // (using max. diffusion and 2D CFL) are limiting the mass-continuity
    // time-step and the third (3D CFL) limits the energy and age time-steps.

    // The mass-continuity time-step is usually smaller, and the skipping
    // mechanism lets us do several mass-continuity steps for each energy step.

    // When -no_mass is set, mass-continuity-related time-step restrictions are
    // disabled, making "skipping" unnecessary.

    // This is why the following two lines appear here and are executed only
    // if do_mass_continuity == true.
    if (do_skip == PETSC_TRUE && skipCountDown > 0)
      skipCountDown--;
    stdout_flags += "h";
  } else {
    stdout_flags += "$";
  }

  grid.profiler->end(event_mass);

  //! \li compute the bed deformation, which only depends on current thickness
  //! and bed elevation
  if (beddef) {
    grid.profiler->begin(event_beddef);
    int topg_state_counter = vbed.get_state_counter();

    ierr = beddef->update(grid.time->current(), dt); CHKERRQ(ierr);

    if (vbed.get_state_counter() != topg_state_counter) {
      stdout_flags += "b";
      ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);
    } else
      stdout_flags += " ";
    grid.profiler->end(event_beddef);
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
  char tempstr[5];  snprintf(tempstr,5," %c", adaptReasonFlag);
  stdout_flags += tempstr;

#if (PISM_DEBUG==1)
  ierr = variables.check_for_nan(); CHKERRQ(ierr);
#endif

  grid.profiler->end(event_step);

  return 0;
}


//! Do the time-stepping for an evolution run.
/*! 
This procedure is the main time-stepping loop.
 */
PetscErrorCode IceModel::run() {
  PetscErrorCode  ierr;

  bool do_mass_conserve = config.get_flag("do_mass_conserve"),
    do_energy = config.get_flag("do_energy"),
    do_age = config.get_flag("do_age"),
    do_skip = config.get_flag("do_skip");
  int stepcount = (config.get_flag("count_time_steps")) ? 0 : -1;

  // do a one-step diagnostic run:
  ierr = verbPrintf(2,grid.com,
      "doing preliminary step to fill diagnostic quantities ...\n"); CHKERRQ(ierr);

  // set verbosity to 1 to suppress reporting
  PetscInt tmp_verbosity = getVerbosityLevel(); 
  ierr = setVerbosityLevel(1); CHKERRQ(ierr);

  dt_force = -1.0;
  maxdt_temporary = -1.0;
  skipCountDown = 0;
  dt_TempAge = 0.0;
  dt = 0.0;
  PetscReal run_end = grid.time->end();

  // FIXME:  In the case of derived class diagnostic time series this fixed
  //         step-length can be problematic.  The fix may have to be in the derived class.
  //         The problem is that unless the derived class fully reinitializes its
  //         time series then there can be a request for an interpolation on [A,B]
  //         where A>B.  See IcePSTexModel.
  grid.time->set_end(grid.time->start() + 1); // run for 1 second

  ierr = step(do_mass_conserve, do_energy, do_age, do_skip); CHKERRQ(ierr);

  // print verbose messages according to user-set verbosity
  if (tmp_verbosity > 2) {
    ierr = PetscPrintf(grid.com,
      " done; reached time %.4f a\n", grid.time->year()); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com,
      "  re-setting model state as initialized ...\n"); CHKERRQ(ierr);
  }

  // re-initialize the model:
  global_attributes.set_string("history", "");
  grid.time->set(grid.time->start());
  t_TempAge = grid.time->start();
  grid.time->set_end(run_end);
  ierr = model_state_setup(); CHKERRQ(ierr);

  // restore verbosity:
  ierr = setVerbosityLevel(tmp_verbosity); CHKERRQ(ierr);

  // Write snapshots and time-series at the beginning of the run.
  ierr = write_snapshot(); CHKERRQ(ierr);
  ierr = write_timeseries(); CHKERRQ(ierr);
  ierr = write_extras(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "running forward ...\n"); CHKERRQ(ierr);

  stdout_flags.erase(); // clear it out
  ierr = summaryPrintLine(PETSC_TRUE,do_energy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); CHKERRQ(ierr);
  adaptReasonFlag = '$'; // no reason for no timestep
  reset_counters();
  ierr = summary(do_energy); CHKERRQ(ierr);  // report starting state

  // main loop for time evolution
  // IceModel::step calls grid.time->step(dt), ensuring that this while loop
  // will terminate
  while (grid.time->current() < grid.time->end()) {

    stdout_flags.erase();  // clear it out
    dt_force = -1.0;
    maxdt_temporary = -1.0;

    ierr = step(do_mass_conserve, do_energy, do_age, do_skip); CHKERRQ(ierr);

    // report a summary for major steps or the last one
    bool updateAtDepth = (skipCountDown == 0);
    bool tempAgeStep = ( updateAtDepth && ((do_energy) || (do_age)) );

    const bool show_step = tempAgeStep || (adaptReasonFlag == 'e');
    ierr = summary(show_step); CHKERRQ(ierr);

    // writing these fields here ensures that we do it after the last time-step
    ierr = write_snapshot(); CHKERRQ(ierr);
    ierr = write_timeseries(); CHKERRQ(ierr);
    ierr = write_extras(); CHKERRQ(ierr);
    ierr = write_backup(); CHKERRQ(ierr);

    ierr = update_viewers(); CHKERRQ(ierr);

    if (stepcount >= 0) stepcount++;
    if (endOfTimeStepHook() != 0) break;
  } // end of the time-stepping loop

  bool flag;
  PetscInt pause_time = 0;
  ierr = PISMOptionsInt("-pause", "Pause after the run, seconds",
			pause_time, flag); CHKERRQ(ierr);
  if (pause_time > 0) {
    ierr = verbPrintf(2,grid.com,"pausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);
  }

  if (stepcount >= 0) {
    ierr = verbPrintf(1,grid.com,
                      "count_time_steps:  run() took %d steps\n"
                      "average dt = %.6f years\n",
                      stepcount, grid.time->run_length_years()/(double)stepcount); CHKERRQ(ierr);
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
  ierr = PISMOptionsIsSet("-wait_for_gdb", wait_for_gdb); CHKERRQ(ierr);
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

  //! 5) Initialize atmosphere and ocean couplers:
  ierr = init_couplers(); CHKERRQ(ierr);

  //! 6) Allocate work vectors:
  ierr = allocate_internal_objects(); CHKERRQ(ierr);

  //! 7) Fill the model state variables (from a PISM output file, from a
  //! bootstrapping file using some modeling choices or using formulas). Calls
  //! IceModel::regrid()
  ierr = model_state_setup(); CHKERRQ(ierr);

  //! 8) Report grid parameters:
  ierr = grid.report_parameters(); CHKERRQ(ierr);

  //! 9) Miscellaneous stuff: set up the bed deformation model, initialize the
  //! basal till model, initialize snapshots. This has to happen *after*
  //! regridding.
  ierr = misc_setup();

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  //! The following flow-chart illustrates the process.
  //!
  //! \dotfile initialization-sequence.dot "IceModel initialization sequence"

  ierr = MPI_Barrier(grid.com); CHKERRQ(ierr);
  ierr = PetscGetTime(&start_time); CHKERRQ(ierr);

  return 0; 
}

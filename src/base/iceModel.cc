// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <petscda.h>

#include "iceModel.hh"
#include "pism_signal.h"

IceModel::IceModel(IceGrid &g, NCConfigVariable &conf, NCConfigVariable &conf_overrides)
  : grid(g), config(conf), overrides(conf_overrides), iceFactory(grid.com,NULL, conf), ice(NULL) {
  PetscErrorCode ierr;

  if (utIsInit() == 0) {
    if (utInit(NULL) != 0) {
      PetscPrintf(grid.com, "PISM ERROR: UDUNITS initialization failed.\n");
      PetscEnd();
    }
  }

  mapping.init("mapping", grid.com, grid.rank);
  global_attributes.init("global_attributes", grid.com, grid.rank);

  pism_signal = 0;
  signal(SIGTERM, pism_signal_handler);
  signal(SIGUSR1, pism_signal_handler);
  signal(SIGUSR2, pism_signal_handler);

  doAdaptTimeStep = PETSC_TRUE;
  basal = NULL;
  CFLviolcount = 0;

  surface = NULL;
  ocean   = NULL;
  beddef  = NULL;

  EC = NULL;
  sia_bed_smoother = NULL;

  ierr = setDefaults();  // lots of parameters and flags set here, including by reading from a config file
  if (ierr != 0) {
    verbPrintf(1,grid.com, "Error setting defaults.\n");
    PetscEnd();
  }

  // Special diagnostic viewers are off by default:
  view_diffusivity = false;
  view_nuH = false;
  view_log_nuH = false;

  // Do not save snapshots by default:
  save_snapshots = false;
  // Do not save time-series by default:
  save_ts = false;
  save_extra = false;

  dvoldt = gdHdtav = 0;
  total_surface_ice_flux = 0;
  total_basal_ice_flux = 0;
  total_sub_shelf_ice_flux = 0;

  have_ssa_velocities = 0;	// no SSA velocities at the start of the run

  allowAboveMelting = PETSC_FALSE;  // only IceCompModel ever sets it to true

  // Default ice type:
  iceFactory.setType(ICE_PB);

  prof = NULL;
}


IceModel::~IceModel() {

  deallocate_internal_objects();

  // write (and deallocate) time-series
  vector<DiagnosticTimeseries*>::iterator i;
  for (i = timeseries.begin(); i < timeseries.end(); ++i)
    delete (*i);

  delete ocean;
  delete surface;
  delete beddef;

  delete basal;
  delete EC;
  delete ice;
  delete prof;

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
  PetscInt WIDE_STENCIL = 2;

  ierr = verbPrintf(3, grid.com,
		    "Allocating memory...\n"); CHKERRQ(ierr);

  // The following code creates (and documents -- to some extent) the
  // variables. The main (and only) principle here is using standard names from
  // the CF conventions; see
  // http://cf-pcmdi.llnl.gov/documents/cf-standard-names

  ierr =     u3.create(grid, "uvel", true); CHKERRQ(ierr);
  ierr =     u3.set_attrs("diagnostic", "horizontal velocity of ice in the X direction",
			  "m s-1", "land_ice_x_velocity"); CHKERRQ(ierr);
  ierr =     u3.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  u3.write_in_glaciological_units = true;
  ierr = variables.add(u3); CHKERRQ(ierr);

  ierr =     v3.create(grid, "vvel", true); CHKERRQ(ierr);
  ierr =     v3.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
			  "m s-1", "land_ice_y_velocity"); CHKERRQ(ierr);
  ierr =     v3.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  v3.write_in_glaciological_units = true;
  ierr = variables.add(v3); CHKERRQ(ierr);

  ierr =     w3.create(grid, "wvel_rel", false); CHKERRQ(ierr); // never diff'ed in hor dirs
  ierr =     w3.set_attrs(
                "diagnostic", 
                "vertical velocity of ice, relative to base of ice directly below point",
		"m s-1", ""); CHKERRQ(ierr);
  ierr =     w3.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  w3.write_in_glaciological_units = true;
  ierr = variables.add(w3); CHKERRQ(ierr);

  ierr = Sigma3.create(grid, "strainheat", false); CHKERRQ(ierr); // never diff'ed in hor dirs
  ierr = Sigma3.set_attrs("internal",
                          "rate of strain heating in ice (dissipation heating)",
	        	  "W m-3", ""); CHKERRQ(ierr);
  ierr = Sigma3.set_glaciological_units("mW m-3"); CHKERRQ(ierr);
  ierr = variables.add(Sigma3); CHKERRQ(ierr);

  if (config.get_flag("do_cold_ice_methods")) {
    // ice temperature
    ierr = T3.create(grid, "temp", true); CHKERRQ(ierr);
    ierr = T3.set_attrs("model_state","ice temperature",
			"K", "land_ice_temperature"); CHKERRQ(ierr);
    ierr = T3.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    ierr = variables.add(T3); CHKERRQ(ierr);
  }

  ierr = Enth3.create(grid, "enthalpy", true, WIDE_STENCIL); CHKERRQ(ierr);
  // POSSIBLE standard name = land_ice_enthalpy
  ierr = Enth3.set_attrs(
     "model_state",
     "ice enthalpy (includes sensible heat, latent heat, pressure)",
     "J kg-1", ""); CHKERRQ(ierr);
  ierr = variables.add(Enth3); CHKERRQ(ierr);

  // age of ice but only if age will be computed
  if (config.get_flag("do_age")) {
    ierr = tau3.create(grid, "age", true, WIDE_STENCIL); CHKERRQ(ierr);
    // PROPOSED standard_name = land_ice_age
    ierr = tau3.set_attrs("diagnostic", "age of ice",
                          "s", ""); CHKERRQ(ierr);
    ierr = tau3.set_glaciological_units("years");
    tau3.write_in_glaciological_units = true;
    ierr = tau3.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    ierr = variables.add(tau3); CHKERRQ(ierr);
  }

  // bedrock temperature
  ierr = Tb3.create(grid,"litho_temp", false); CHKERRQ(ierr);
  // PROPOSED standard_name = lithosphere_temperature
  ierr = Tb3.set_attrs("model_state", "lithosphere (bedrock) temperature",
		       "K", ""); CHKERRQ(ierr);
  ierr = Tb3.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = variables.add(Tb3); CHKERRQ(ierr);

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

  // grounded_dragging_floating integer mask
  ierr = vMask.create(grid, "mask", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = vMask.set_attrs("model_state", "grounded_dragging_floating integer mask",
			 "", ""); CHKERRQ(ierr);
  vector<double> mask_values(6);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_SHEET;
  mask_values[2] = MASK_DRAGGING_SHEET;
  mask_values[3] = MASK_FLOATING;
  mask_values[4] = MASK_ICE_FREE_OCEAN;
  mask_values[5] = MASK_OCEAN_AT_TIME_0;
  ierr = vMask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
  ierr = vMask.set_attr("flag_meanings",
			"ice_free_bedrock sheet dragging_sheet floating ice_free_ocean ocean_at_time_zero"); CHKERRQ(ierr);
  vMask.output_data_type = NC_BYTE;
  ierr = variables.add(vMask); CHKERRQ(ierr);

  // upward geothermal flux at bedrock surface
  ierr = vGhf.create(grid, "bheatflx", false); CHKERRQ(ierr); // never differentiated
  // PROPOSED standard_name = lithosphere_upward_heat_flux
  ierr = vGhf.set_attrs("climate_steady", "upward geothermal flux at bedrock surface",
			"W m-2", ""); CHKERRQ(ierr);
  ierr = vGhf.set_glaciological_units("mW m-2");
  vGhf.time_independent = true;
  ierr = variables.add(vGhf); CHKERRQ(ierr);

  // u bar and v bar
  ierr = vel_bar.create(grid, "bar", false); CHKERRQ(ierr); // components are ubar and vbar
  ierr = vel_bar.set_attrs("diagnostic", 
			   "vertical mean of horizontal ice velocity in the X direction",
			   "m s-1", "land_ice_vertical_mean_x_velocity", 0); CHKERRQ(ierr);
  ierr = vel_bar.set_attrs("diagnostic", 
			   "vertical mean of horizontal ice velocity in the Y direction",
			   "m s-1", "land_ice_vertical_mean_y_velocity", 1); CHKERRQ(ierr);
  ierr = vel_bar.set_glaciological_units("m year-1");
  vel_bar.write_in_glaciological_units = true;
  ierr = variables.add(vel_bar, "uvbar"); CHKERRQ(ierr);

  // basal velocities on standard grid
  ierr = vel_basal.create(grid, "basal", true); CHKERRQ(ierr); // components are ubasal and vbasal
  ierr = vel_basal.set_attrs("diagnostic", "basal ice velocity in the X direction",
		       "m s-1", "land_ice_basal_x_velocity", 0); CHKERRQ(ierr);
  ierr = vel_basal.set_attrs("diagnostic", "basal ice velocity in the Y direction",
		       "m s-1", "land_ice_basal_y_velocity", 1); CHKERRQ(ierr);
  ierr = vel_basal.set_glaciological_units("m year-1");
  vel_basal.write_in_glaciological_units = true;
  ierr = variables.add(vel_basal, "uvbasal"); CHKERRQ(ierr);

  // basal frictional heating on regular grid
  ierr = vRb.create(grid, "bfrict", false); CHKERRQ(ierr);
  // PROPOSED standard_name = land_ice_basal_frictional_heating
  ierr = vRb.set_attrs("diagnostic",
                       "basal frictional heating from ice sliding (= till dissipation)",
		       "W m-2", ""); CHKERRQ(ierr);
  ierr = vRb.set_glaciological_units("mW m-2");
  vRb.write_in_glaciological_units = true;
  ierr = vRb.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = variables.add(vRb); CHKERRQ(ierr);

  // effective thickness of subglacial melt water
  ierr = vHmelt.create(grid, "bwat", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = vHmelt.set_attrs("model_state", "effective thickness of subglacial melt water",
			  "m", ""); CHKERRQ(ierr);
  // NB! Effective thickness of subglacial melt water *does* vary from 0 to hmelt_max meters only.
  ierr = vHmelt.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = vHmelt.set_attr("valid_max", config.get("hmelt_max")); CHKERRQ(ierr);
  ierr = variables.add(vHmelt); CHKERRQ(ierr);

  // rate of change of ice thickness
  ierr = vdHdt.create(grid, "dHdt", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = vdHdt.set_attrs("diagnostic", "rate of change of ice thickness",
			 "m s-1", "tendency_of_land_ice_thickness"); CHKERRQ(ierr);
  ierr = vdHdt.set_glaciological_units("m year-1");
  vdHdt.write_in_glaciological_units = true;
  const PetscScalar  huge_dHdt = 1.0e6;      // million m a-1 is out-of-range
  ierr = vdHdt.set_attr("valid_min", -huge_dHdt / secpera); CHKERRQ(ierr);
  ierr = vdHdt.set_attr("valid_max", huge_dHdt / secpera); CHKERRQ(ierr);
  ierr = variables.add(vdHdt); CHKERRQ(ierr);

  // yield stress for basal till (plastic or pseudo-plastic model)
  ierr = vtauc.create(grid, "tauc", false); CHKERRQ(ierr);
  // PROPOSED standard_name = land_ice_basal_material_yield_stress
  ierr = vtauc.set_attrs("diagnostic", 
             "yield stress for basal till (plastic or pseudo-plastic model)",
	     "Pa", ""); CHKERRQ(ierr);
  ierr = variables.add(vtauc); CHKERRQ(ierr);

  // bedrock uplift rate
  ierr = vuplift.create(grid, "dbdt", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = vuplift.set_attrs("model_state", "bedrock uplift rate",
			   "m s-1", "tendency_of_bedrock_altitude"); CHKERRQ(ierr);
  ierr = vuplift.set_glaciological_units("m year-1");
  vuplift.write_in_glaciological_units = true;
  ierr = variables.add(vuplift); CHKERRQ(ierr);

  // basal melt rate
  ierr = vbmr.create(grid, "bmelt", false); CHKERRQ(ierr);
  ierr = vbmr.set_attrs("model_state",
                                  "ice basal melt rate in ice thickness per time",
				  "m s-1", "land_ice_basal_melt_rate"); CHKERRQ(ierr);
  ierr = vbmr.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vbmr.write_in_glaciological_units = true;
  vbmr.set_attr("comment", "positive basal melt rate corresponds to ice loss");
  ierr = variables.add(vbmr); CHKERRQ(ierr);

  // friction angle for till under grounded ice sheet
  ierr = vtillphi.create(grid, "tillphi", false); // never differentiated
  // PROPOSED standard_name = land_ice_basal_material_friction_angle
  ierr = vtillphi.set_attrs("climate_steady", "friction angle for till under grounded ice sheet",
			    "degrees", ""); CHKERRQ(ierr);
  vtillphi.time_independent = true;
  ierr = variables.add(vtillphi); CHKERRQ(ierr);

  // longitude
  ierr = vLongitude.create(grid, "lon", true); CHKERRQ(ierr);
  ierr = vLongitude.set_attrs("mapping", "longitude", "degree_east", "longitude"); CHKERRQ(ierr);
  vLongitude.time_independent = true;
  ierr = vLongitude.set_attr("coordinates", ""); CHKERRQ(ierr);
  ierr = vLongitude.set_attr("grid_mapping", ""); CHKERRQ(ierr);
  ierr = variables.add(vLongitude); CHKERRQ(ierr);

  // latitude
  ierr = vLatitude.create(grid, "lat", true); CHKERRQ(ierr); // has ghosts so that we can compute cell areas
  ierr = vLatitude.set_attrs("mapping", "latitude", "degree_north", "latitude"); CHKERRQ(ierr);
  vLatitude.time_independent = true;
  ierr = vLatitude.set_attr("coordinates", ""); CHKERRQ(ierr);
  ierr = vLatitude.set_attr("grid_mapping", ""); CHKERRQ(ierr);
  ierr = variables.add(vLatitude); CHKERRQ(ierr);

  // cell areas
  ierr = cell_area.create(grid, "cell_area", false); CHKERRQ(ierr);
  ierr = cell_area.set_attrs("diagnostic", "cell areas", "m2", ""); CHKERRQ(ierr);
  ierr = cell_area.set_attr("comment",
                            "values are equal to dx*dy "
                            "if latitude and longitude fields are not available; "
                            "otherwise WGS84 ellipsoid is used"); CHKERRQ(ierr); 
  cell_area.time_independent = true;
  ierr = cell_area.set_glaciological_units("km2"); CHKERRQ(ierr);
  cell_area.write_in_glaciological_units = true;
  ierr = variables.add(cell_area); CHKERRQ(ierr);

  // u bar and v bar on staggered grid
  ierr = uvbar.create(grid, "uvbar", true); CHKERRQ(ierr);
  ierr = uvbar.set_attrs("internal", 
			 "vertically averaged ice velocity, on staggered grid offset in X direction,"
			 " from SIA, in the X direction",
			 "m s-1", "", 0); CHKERRQ(ierr);
  ierr = uvbar.set_attrs("internal", 
			 "vertically averaged ice velocity, on staggered grid offset in Y direction,"
			 " from SIA, in the Y direction",
			 "m s-1", "", 1); CHKERRQ(ierr);

  // initial guesses of SSA velocities
  ierr = vel_ssa.create(grid, "bar_ssa", true, WIDE_STENCIL); // components are ubar_ssa and vbar_ssa
  ierr = vel_ssa.set_attrs("internal_restart", "SSA model ice velocity in the X direction",
                            "m s-1", "", 0); CHKERRQ(ierr);
  ierr = vel_ssa.set_attrs("internal_restart", "SSA model ice velocity in the Y direction",
                            "m s-1", "", 1); CHKERRQ(ierr);
  ierr = vel_ssa.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  ierr = variables.add(vel_ssa, "uvbar_ssa"); CHKERRQ(ierr);

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
  ierr = variables.add(acab); CHKERRQ(ierr);

  // annual mean air temperature at "ice surface", at level below all firn
  //   processes (e.g. "10 m" or ice temperatures)
  ierr = artm.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = artm.set_attrs(
            "climate_from_PISMSurfaceModel",  // FIXME: can we do better?
            "annual average ice surface temperature, below firn processes",
            "K", 
            "");  // PROPOSED CF standard_name = land_ice_surface_temperature_below_firn
            CHKERRQ(ierr);
  ierr = variables.add(artm); CHKERRQ(ierr);

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
  ierr = variables.add(shelfbmassflux); CHKERRQ(ierr);

  // ice boundary tempature at the base of the ice shelf
  ierr = shelfbtemp.create(grid, "shelfbtemp", false); CHKERRQ(ierr); // no ghosts; NO HOR. DIFF.!
  ierr = shelfbtemp.set_attrs(
           "climate_state", "absolute temperature at ice shelf base",
	   "K", ""); CHKERRQ(ierr);
  // PROPOSED standard name = ice_shelf_basal_temperature
  ierr = variables.add(shelfbtemp); CHKERRQ(ierr);

  return 0;
}


//! De-allocate internal objects.
/*! This includes Vecs that are not in an IceModelVec, SSA tools and the bed
  deformation model.
 */
PetscErrorCode IceModel::deallocate_internal_objects() {
  PetscErrorCode ierr;

  // since SSA tools are currently part of IceModel, de-allocate them here
  ierr = destroySSAobjects(); CHKERRQ(ierr);

  // SIA has Schoof (2003)-type smoother, allocate it here
  delete sia_bed_smoother;

  return 0;
}

void IceModel::setConstantNuHForSSA(PetscScalar nuH) {
  config.set_flag("use_constant_nuh_for_ssa", true);
  ssaStrengthExtend.set_notional_strength(nuH);
}


PetscErrorCode IceModel::setExecName(const char *my_executable_short_name) {
  executable_short_name = my_executable_short_name;
  return 0;
}

//! The contents of the time-step.
/*!
During the time-step we perform the following actions:
 */
PetscErrorCode IceModel::step(bool do_mass_continuity,
			      bool do_energy,
			      bool do_age,
			      bool do_skip,
			      bool use_ssa_when_grounded) {
  PetscErrorCode ierr;

  //! \li call additionalAtStartTimestep() to let derived classes do more
  ierr = additionalAtStartTimestep(); CHKERRQ(ierr);  // might set dt_force,maxdt_temporary

  //! \li determine the maximum time-step boundary models can take
  double apcc_dt;
  ierr = surface->max_timestep(grid.year, apcc_dt); CHKERRQ(ierr);
  apcc_dt *= secpera;
  if (apcc_dt > 0.0) {
    if (maxdt_temporary > 0)
      maxdt_temporary = PetscMin(apcc_dt, maxdt_temporary);
    else
      maxdt_temporary = apcc_dt;
  }

  double opcc_dt;
  ierr = ocean->max_timestep(grid.year, opcc_dt); CHKERRQ(ierr);
  opcc_dt *= secpera;
  if (opcc_dt > 0.0) {
    if (maxdt_temporary > 0)
      maxdt_temporary = PetscMin(opcc_dt, maxdt_temporary);
    else
      maxdt_temporary = opcc_dt;
  }

  //! \li apply the time-step restriction from the -extra_{times,file,vars}
  //! mechanism
  double extras_dt;
  ierr = extras_max_timestep(grid.year, extras_dt); CHKERRQ(ierr);
  extras_dt *= secpera;
  if (extras_dt > 0.0) {
    if (maxdt_temporary > 0)
      maxdt_temporary = PetscMin(extras_dt, maxdt_temporary);
    else
      maxdt_temporary = extras_dt;
  }

  PetscLogEventBegin(beddefEVENT,0,0,0,0);

  prof->begin(event_beddef);
  //! \li compute the bed deformation, which only depends on current thickness
  //! and bed elevation
  if (beddef) {
    bool bed_changed;
    ierr = bed_def_step(bed_changed); CHKERRQ(ierr);
    if (bed_changed) {
      stdout_flags += "b";
      // since bed changed we need to update info in bed smoother
      ierr = sia_bed_smoother->preprocess_bed(vbed,
               config.get("Glen_exponent"), config.get("bed_smoother_range") );
               CHKERRQ(ierr);
    } else stdout_flags += "$";
  } else stdout_flags += " ";
  prof->end(event_beddef);

  PetscLogEventEnd(beddefEVENT,0,0,0,0);

  //! \li update the yield stress for the plastic till model (if appropriate)
  if (use_ssa_when_grounded) {
    ierr = updateYieldStressUsingBasalWater();  CHKERRQ(ierr);
    stdout_flags += "y";
  } else stdout_flags += "$";

  //! \li update the velocity field; in some cases the whole three-dimensional
  //! field is updated and in some cases just the vertically-averaged
  //! horizontal velocity is updated; see velocity()

  // always do SIA velocity calculation; only update SSA and 
  //   only update velocities at depth if suggested by temp and age
  //   stability criterion; note *lots* of communication is avoided by 
  //   skipping SSA (and temp/age)

  prof->begin(event_velocity);

  bool updateAtDepth = (skipCountDown == 0);
  ierr = velocity(updateAtDepth); CHKERRQ(ierr);  // event logging in here
  stdout_flags += (updateAtDepth ? "v" : "V");

  prof->end(event_velocity);

  //! \li determine the time step according to a variety of stability criteria;
  //!  see determineTimeStep()
  ierr = determineTimeStep(do_energy); CHKERRQ(ierr);

  dtTempAge += dt;
  grid.year += dt / secpera;  // adopt it
  // IceModel::dt,dtTempAge,grid.year are now set correctly according to
  //    mass-continuity-eqn-diffusivity criteria, horizontal CFL criteria, and
  //    other criteria from derived class additionalAtStartTimestep(), and from
  //    "-skip" mechanism

  PetscLogEventBegin(tempEVENT,0,0,0,0);

  prof->begin(event_age);
  //! \li update the age of the ice (if appropriate)
  if (do_age) {
    ierr = ageStep(); CHKERRQ(ierr);
    stdout_flags += "a";
  } else {
    stdout_flags += "$";
  }
  prof->end(event_age);

  prof->begin(event_energy);

  //! \li update the enthalpy (or temperature) field according to the conservation of
  //!  energy model based (especially) on the new velocity field; see
  //!  energyStep()
  if (updateAtDepth && do_energy) { // do the temperature step
    ierr = energyStep(); CHKERRQ(ierr);
    if (updateHmelt == PETSC_TRUE) {
      ierr = diffuseHmelt(); CHKERRQ(ierr);
    }
    stdout_flags += "t";
  } else {
    stdout_flags += "$";
  }

  dtTempAge = 0.0;

  prof->end(event_energy);

  PetscLogEventEnd(tempEVENT,0,0,0,0);

  //! \li compute fluxes through ice boundaries; this method frequently updates
  //! the surface process models pre-emptively, so that massContExplicitStep()
  //! does not do that work again
  ierr = ice_mass_bookkeeping(); CHKERRQ(ierr);

  PetscLogEventBegin(massbalEVENT,0,0,0,0);

  prof->begin(event_mass);
  //! \li update the thickness of the ice according to the mass conservation
  //!  model; see massContExplicitStep()
  if (do_mass_continuity) {
    ierr = massContExplicitStep(); CHKERRQ(ierr); // update H
    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr); // update h and mask
    if ((do_skip == PETSC_TRUE) && (skipCountDown > 0))
      skipCountDown--;
    stdout_flags += "h";
  } else {
    stdout_flags += "$";
    // if do_mass_conserve is false, then ice thickness does not change and
    // dH/dt = 0:
    ierr = vdHdt.set(0.0); CHKERRQ(ierr);
  }
  prof->end(event_mass);

  PetscLogEventEnd(massbalEVENT,0,0,0,0);

  //! \li call additionalAtEndTimestep() to let derived classes do more
  ierr = additionalAtEndTimestep(); CHKERRQ(ierr);

  // end the flag line
  char tempstr[5];  snprintf(tempstr,5," %c", adaptReasonFlag);
  stdout_flags += tempstr;

#ifdef PISM_DEBUG
  ierr = variables.check_for_nan(); CHKERRQ(ierr);
#endif

  return 0;
}


//! Do the time-stepping for an evolution run.
/*! 
This procedure is the main time-stepping loop.
 */
PetscErrorCode IceModel::run() {
  PetscErrorCode  ierr;

  bool do_mass_conserve = config.get_flag("do_mass_conserve"),
    do_temp = config.get_flag("do_temp"),
    do_age = config.get_flag("do_age"),
    do_skip = config.get_flag("do_skip"),
    use_ssa_when_grounded = config.get_flag("use_ssa_when_grounded");

#if PETSC_VERSION_MAJOR >= 3
# define PismLogEventRegister(name,cookie,event) PetscLogEventRegister((name),(cookie),(event))
#else
# define PismLogEventRegister(name,cookie,event) PetscLogEventRegister((event),(name),(cookie))
#endif
PismLogEventRegister("sia velocity", 0,&siaEVENT);
PismLogEventRegister("ssa velocity", 0,&ssaEVENT);
PismLogEventRegister("misc vel calc",0,&velmiscEVENT);
PismLogEventRegister("bed deform",   0,&beddefEVENT);
PismLogEventRegister("mass bal calc",0,&massbalEVENT);
PismLogEventRegister("temp age calc",0,&tempEVENT);

  // do a one-step diagnostic run:
  ierr = verbPrintf(2,grid.com,
      "doing preliminary step to fill diagnostic quantities ...\n"); CHKERRQ(ierr);

  // set verbosity to 1 to suppress reporting
  PetscInt tmp_verbosity = getVerbosityLevel(); 
  ierr = setVerbosityLevel(1); CHKERRQ(ierr);

  dt_force = -1.0;
  maxdt_temporary = -1.0;
  skipCountDown = 0;
  dtTempAge = 0.0;
  dt = 0.0;
  PetscReal end_year = grid.end_year;
  
  // FIXME:  In the case of derived class diagnostic time series this fixed
  //         step-length can be problematic.  The fix may have to be in the derived class.
  //         The problem is that unless the derived class fully reinitializes its
  //         time series then there can be a request for an interpolation on [A,B]
  //         where A>B.  See IcePSTexModel.
  //grid.end_year = grid.start_year + 1; // all what matters is that it is
  //				       // greater than start_year
  grid.end_year = grid.start_year + 0.01; // all what matters is that it is
				       // greater than start_year

  prof->begin(event_step);
  ierr = step(do_mass_conserve, do_temp, do_age,
	      do_skip, use_ssa_when_grounded); CHKERRQ(ierr);
  prof->end(event_step);

  // print verbose messages according to user-set verbosity
  if (tmp_verbosity > 2) {
    ierr = PetscPrintf(grid.com,
      " done; reached time %.4f a\n", grid.year); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com,
      "  re-setting model state as initialized ...\n"); CHKERRQ(ierr);
  }

  // re-initialize the model:
  global_attributes.set_string("history", "");
  grid.year = grid.start_year;
  grid.end_year = end_year;
  ierr = model_state_setup(); CHKERRQ(ierr);

  // restore verbosity:
  ierr = setVerbosityLevel(tmp_verbosity); CHKERRQ(ierr);

  // Write snapshots and time-series at the beginning of the run.
  ierr = write_snapshot(); CHKERRQ(ierr);
  ierr = write_timeseries(); CHKERRQ(ierr);
  ierr = write_extras(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, "running forward ...\n"); CHKERRQ(ierr);

  stdout_flags.erase(); // clear it out
  ierr = summaryPrintLine(PETSC_TRUE,do_temp, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); CHKERRQ(ierr);
  adaptReasonFlag = '$'; // no reason for no timestep
  skipCountDown = 0;
  dtTempAge = 0.0;
  maxdt_temporary = dt = dt_force = 0.0;
  dt_from_diffus = dt_from_cfl = 0.0;
  CFLmaxdt = CFLmaxdt2D = 0.0;
  gDmax = dvoldt = gdHdtav = 0;
  total_surface_ice_flux = 0;
  total_basal_ice_flux = 0;
  total_sub_shelf_ice_flux = 0;
  gmaxu = gmaxv = gmaxw = -1;
  ierr = summary(do_temp,reportPATemps); CHKERRQ(ierr);  // report starting state

  // main loop for time evolution
  for (PetscScalar year = grid.start_year; year < grid.end_year; year += dt/secpera) {
    
    stdout_flags.erase();  // clear it out
    dt_force = -1.0;
    maxdt_temporary = -1.0;

    prof->begin(event_step);
    ierr = step(do_mass_conserve, do_temp, do_age,
		do_skip, use_ssa_when_grounded); CHKERRQ(ierr);
    prof->end(event_step);

    // report a summary for major steps or the last one
    bool updateAtDepth = (skipCountDown == 0);
    bool tempAgeStep = ( updateAtDepth && ((do_temp) || (do_age)) );

    const bool show_step = tempAgeStep || (adaptReasonFlag == 'e');
    ierr = summary(show_step,reportPATemps); CHKERRQ(ierr);

    // writing these fields here ensures that we do it after the last time-step
    ierr = write_snapshot(); CHKERRQ(ierr);
    ierr = write_timeseries(); CHKERRQ(ierr);
    ierr = write_extras(); CHKERRQ(ierr);

    ierr = update_viewers(); CHKERRQ(ierr);

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

  //! 4) Initialize the IceFlowLaw and (possibly) other physics.
  ierr = init_physics(); CHKERRQ(ierr);

  //! 5) Initialize atmosphere and ocean couplers:
  ierr = init_couplers(); CHKERRQ(ierr);

  //! 6) Fill the model state variables (from a PISM output file, from a
  //! bootstrapping file using some modeling choices or using formulas). Calls
  //! IceModel::regrid()
  ierr = model_state_setup(); CHKERRQ(ierr);

  //! 7) Report grid parameters:
  ierr = report_grid_parameters(); CHKERRQ(ierr);

  //! 8) Allocate SSA tools and work vectors:
  ierr = allocate_internal_objects(); CHKERRQ(ierr);

  //! 9) Miscellaneous stuff: set up the bed deformation model, initialize the
  //! basal till model, initialize snapshots. This has to happen *after*
  //! regridding.
  ierr = misc_setup();

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  //! The following flow-chart illustrates the process.
  //! \dotfile initialization-sequence.dot IceModel initialization sequence

  ierr = PetscGetTime(&start_time); CHKERRQ(ierr);

  return 0; 
}

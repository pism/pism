// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "ssaJed/pismssa.hh"

IceModel::IceModel(IceGrid &g)
  : grid(g), iceFactory(grid.com,NULL), ice(NULL), shelfExtensionJed(grid.com,NULL) {
  PetscErrorCode ierr;

  if (utIsInit() == 0) {
    if (utInit(NULL) != 0) {
      PetscPrintf(grid.com, "PISM ERROR: UDUNITS initialization failed.\n");
      PetscEnd();
    }
  }

  bootstrapLIC = PETSC_NULL;
  
  history_size = TEMPORARY_STRING_LENGTH;
  history = new char[history_size];
  history[0] = 0;               // Initialize with empty string so that prepending works correctly

  have_ssa_velocities = false;

  pism_signal = 0;
  signal(SIGTERM, pism_signal_handler);
  signal(SIGUSR1, pism_signal_handler);

  createBasal_done = PETSC_FALSE;
  top0ctx = PETSC_NULL;
  g2natural = PETSC_NULL;
  createVecs_done = PETSC_FALSE;
  CFLviolcount = 0;

  for (PetscInt nn = 0; nn < tnN; nn++) {
    runtimeViewers[nn] = PETSC_NULL;
  }
  createViewers_done = PETSC_FALSE;

  atmosPCC = PETSC_NULL;
  info_atmoscoupler.lat = PETSC_NULL;
  info_atmoscoupler.lon = PETSC_NULL;  
  info_atmoscoupler.mask = PETSC_NULL;
  info_atmoscoupler.surfelev = PETSC_NULL;
  oceanPCC = PETSC_NULL;
  info_oceancoupler.lat = PETSC_NULL;
  info_oceancoupler.lon = PETSC_NULL;  
  info_oceancoupler.mask = PETSC_NULL;
  info_oceancoupler.thk = PETSC_NULL;

  ierr = setDefaults();  // lots of parameters and flags set here
  if (ierr != 0) {
    verbPrintf(1,grid.com, "Error setting defaults.\n");
    PetscEnd();
  }

  // Do not save snapshots by default:
  save_snapshots = false;
  dvoldt = gdHdtav = 0;

  // Default ice type:
  iceFactory.setType(ICE_PB);
}


IceModel::~IceModel() {
  // actions to de-allocate Vecs first
  if (createVecs_done == PETSC_TRUE) {
    destroyVecs();
  }
  // other deallocations
  if (createViewers_done == PETSC_TRUE) {
    destroyViewers();
  }
  if (createBasal_done == PETSC_TRUE) {
    delete basal;
    delete basalSIA;
  }
  if (bootstrapLIC != PETSC_NULL) {
    delete bootstrapLIC;
    bootstrapLIC = PETSC_NULL;
  }
  if (ssa) SSADestroy(ssa);
  delete[] history;
  delete ice;
  utTerm(); // Clean up after UDUNITS
}


//! Allocate all Vecs defined in IceModel.
/*!

  This procedure allocates the memory used to store model state, diagnostic and
  work vectors and sets metadata.

  Default values should not be set here; please use set_vars_from_options().
  
*/
PetscErrorCode IceModel::createVecs() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com,
		    "Allocating memory...\n"); CHKERRQ(ierr);

  if (createVecs_done == PETSC_TRUE) {
    SETERRQ(1, "IceModel::createVecs(): IceModelVecs are created already.");
  }

  // The following code creates (and documents -- to some extent) the
  // variables. The main (and only) principle here is using standard names from
  // the CF conventions; see
  // http://cf-pcmdi.llnl.gov/documents/cf-standard-names

  ierr =     u3.create(grid, "uvel", true); CHKERRQ(ierr);
  ierr =     u3.set_attrs("diagnostic", "horizontal velocity of ice in the X direction",
			  "m s-1", "land_ice_x_velocity"); CHKERRQ(ierr);
  ierr =     u3.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  u3.write_in_glaciological_units = true;

  ierr =     v3.create(grid, "vvel", true); CHKERRQ(ierr);
  ierr =     v3.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
			  "m s-1", "land_ice_y_velocity"); CHKERRQ(ierr);
  ierr =     v3.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  v3.write_in_glaciological_units = true;

  ierr =     w3.create(grid, "wvel", false); CHKERRQ(ierr); // never diff'ed in hor dirs
  // PROPOSED standard name = land_ice_upward_velocity
  //   (compare "upward_air_velocity" and "upward_sea_water_velocity")
  ierr =     w3.set_attrs("diagnostic", "vertical velocity of ice",
			  "m s-1", NULL); CHKERRQ(ierr);
  ierr =     w3.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  w3.write_in_glaciological_units = true;

  ierr = Sigma3.create(grid, "Sigma", false); CHKERRQ(ierr); // never diff'ed in hor dirs
  ierr = Sigma3.set_attrs("internal","rate of strain heating",
	        	  "W m-3", NULL); CHKERRQ(ierr);

  // ice temperature
  ierr = T3.create(grid, "temp", true); CHKERRQ(ierr);
  ierr = T3.set_attrs("model_state","ice temperature",
		      "K", "land_ice_temperature"); CHKERRQ(ierr);
  ierr = T3.set_valid_min(0.0); CHKERRQ(ierr);

  // age of ice
  ierr = tau3.create(grid, "age", true); CHKERRQ(ierr);
  // PROPOSED standard_name = land_ice_age
  ierr = tau3.set_attrs("model_state", "age of ice",
			"s", NULL); CHKERRQ(ierr);
  ierr = tau3.set_glaciological_units("years");
  tau3.write_in_glaciological_units = true;
  ierr = tau3.set_valid_min(0.0); CHKERRQ(ierr);

  // bedrock temperature
  ierr = Tb3.create(grid,"litho_temp", false); CHKERRQ(ierr);
  // PROPOSED standard_name = lithosphere_temperature
  ierr = Tb3.set_attrs("model_state", "lithosphere (bedrock) temperature",
		       "K", NULL); CHKERRQ(ierr);
  ierr = Tb3.set_valid_min(0.0); CHKERRQ(ierr);

  // ice upper surface elevation
  ierr = vh.create(grid, "usurf", true); CHKERRQ(ierr);
  ierr = vh.set_attrs("diagnostic", "ice upper surface elevation",
		      "m", "surface_altitude"); CHKERRQ(ierr);

  // land ice thickness
  ierr = vH.create(grid, "thk", true); CHKERRQ(ierr);
  ierr = vH.set_attrs("model_state", "land ice thickness",
		      "m", "land_ice_thickness"); CHKERRQ(ierr);
  ierr = vH.set_valid_min(0.0); CHKERRQ(ierr);

  // bedrock surface elevation
  ierr = vbed.create(grid, "topg", true); CHKERRQ(ierr);
  ierr = vbed.set_attrs("model_state", "bedrock surface elevation",
			"m", "bedrock_altitude"); CHKERRQ(ierr);

  // grounded_dragging_floating integer mask
  ierr = vMask.create(grid, "mask",     true); CHKERRQ(ierr);
  ierr = vMask.set_attrs("model_state", "grounded_dragging_floating integer mask",
			 NULL, NULL); CHKERRQ(ierr);

  // upward geothermal flux at bedrock surface
  ierr = vGhf.create(grid, "bheatflx", true); CHKERRQ(ierr);
  // PROPOSED standard_name = lithosphere_upward_heat_flux
  ierr = vGhf.set_attrs("climate_steady",
                        "upward geothermal flux at bedrock surface",
			"W m-2",
			NULL); CHKERRQ(ierr);
  ierr = vGhf.set_glaciological_units("mW m-2");

  // u bar and v bar
  ierr = vubar.create(grid, "ubar", true); CHKERRQ(ierr);
  ierr = vubar.set_attrs("diagnostic", 
                         "vertical mean of horizontal ice velocity in the X direction",
			 "m s-1", "land_ice_vertical_mean_x_velocity"); CHKERRQ(ierr);
  ierr = vubar.set_glaciological_units("m year-1");
  vubar.write_in_glaciological_units = true;

  ierr = vvbar.create(grid, "vbar", true); CHKERRQ(ierr);
  ierr = vvbar.set_attrs("diagnostic", 
                         "vertical mean of horizontal ice velocity in the Y direction",
			 "m s-1", "land_ice_vertical_mean_y_velocity"); CHKERRQ(ierr);
  ierr = vvbar.set_glaciological_units("m year-1");
  vvbar.write_in_glaciological_units = true;

  // basal velocities on standard grid
  ierr = vub.create(grid, "ub", true); CHKERRQ(ierr);
  ierr = vub.set_attrs("diagnostic", "basal ice velocity in the X direction",
		       "m s-1", "land_ice_basal_x_velocity"); CHKERRQ(ierr);
  ierr = vub.set_glaciological_units("m year-1");
  vub.write_in_glaciological_units = true;
  
  ierr = vvb.create(grid, "vb", true); CHKERRQ(ierr);
  ierr = vvb.set_attrs("diagnostic", "basal ice velocity in the Y direction",
		       "m s-1", "land_ice_basal_y_velocity"); CHKERRQ(ierr);
  ierr = vvb.set_glaciological_units("m year-1");
  vvb.write_in_glaciological_units = true;

  // basal frictional heating on regular grid
  ierr = vRb.create(grid, "Rb", true); CHKERRQ(ierr);
  // PROPOSED standard_name = land_ice_basal_frictional_heating
  ierr = vRb.set_attrs("diagnostic", "basal frictional heating",
		       "W m-2", NULL); CHKERRQ(ierr);
  ierr = vRb.set_glaciological_units("mW m-2");
  vRb.write_in_glaciological_units = true;

  // effective thickness of subglacial melt water
  ierr = vHmelt.create(grid, "bwat", true); CHKERRQ(ierr);
  ierr = vHmelt.set_attrs("model_state", "effective thickness of subglacial melt water",
			  "m", NULL); CHKERRQ(ierr);
  // NB! Effective thickness of subglacial melt water *does* vary from 0 to 2 meters only.
  ierr = vHmelt.set_valid_range(0.0, 2.0); CHKERRQ(ierr);

  // rate of change of ice thickness
  ierr = vdHdt.create(grid, "dHdt", true); CHKERRQ(ierr);
  ierr = vdHdt.set_attrs("diagnostic", "rate of change of ice thickness",
			 "m s-1", "tendency_of_land_ice_thickness"); CHKERRQ(ierr);
  ierr = vdHdt.set_glaciological_units("m year-1");
  vdHdt.write_in_glaciological_units = true;

  // yield stress for basal till (plastic or pseudo-plastic model)
  ierr = vtauc.create(grid, "tauc", true); CHKERRQ(ierr);
  // PROPOSED standard_name = land_ice_basal_material_yield_stress
  ierr = vtauc.set_attrs("diagnostic", 
             "yield stress for basal till (plastic or pseudo-plastic model)",
	     "Pa", NULL); CHKERRQ(ierr);

  // bedrock uplift rate
  ierr = vuplift.create(grid, "dbdt", true); CHKERRQ(ierr);
  ierr = vuplift.set_attrs("model_state", "bedrock uplift rate",
			   "m s-1", "tendency_of_bedrock_altitude"); CHKERRQ(ierr);
  ierr = vuplift.set_glaciological_units("m year-1");
  vuplift.write_in_glaciological_units = true;

  // basal melt rate
  ierr = vbasalMeltRate.create(grid, "basal_melt_rate", true); CHKERRQ(ierr);
  ierr = vbasalMeltRate.set_attrs("diagnostic", "basal melt rate",
				  "m s-1", "land_ice_basal_melt_rate"); CHKERRQ(ierr);
  ierr = vbasalMeltRate.set(0.0); CHKERRQ(ierr);  // so vertical velocities do not use junk from 
                                                  //   uninitialized basal melt rate.
  ierr = vbasalMeltRate.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vbasalMeltRate.write_in_glaciological_units = true;

  // friction angle for till under grounded ice sheet
  ierr = vtillphi.create(grid, "tillphi", true);
  // PROPOSED standard_name = land_ice_basal_material_friction_angle
  ierr = vtillphi.set_attrs("climate_steady", "friction angle for till under grounded ice sheet",
			    "degrees", NULL); CHKERRQ(ierr);

  // longitude
  ierr = vLongitude.create(grid, "lon", true); CHKERRQ(ierr);
  ierr = vLongitude.set_attrs("mapping", "longitude", "degree_east", "longitude"); CHKERRQ(ierr);

  // latitude
  ierr = vLatitude.create(grid, "lat", true); CHKERRQ(ierr);
  ierr = vLatitude.set_attrs("mapping", "latitude", "degree_north", "latitude"); CHKERRQ(ierr);

  // u bar and v bar on staggered grid
  ierr = vuvbar[0].create(grid, "vuvbar[0]", true); CHKERRQ(ierr);
  ierr = vuvbar[0].set_attrs("internal", 
            "vertically averaged ice velocity, on staggered grid offset in X direction, from SIA, in the X direction",
	    "m s-1", NULL); CHKERRQ(ierr);
  ierr = vuvbar[1].create(grid, "vuvbar[1]", true); CHKERRQ(ierr);
  ierr = vuvbar[1].set_attrs("internal", 
            "vertically averaged ice velocity, on staggered grid offset in Y direction, from SIA, in the Y direction",
	    "m s-1", NULL); CHKERRQ(ierr);

  // initial guesses of SSA velocities
  ierr = vubarSSA.create(grid, "vubarSSA", true);
  ierr = vubarSSA.set_attrs("internal_restart", "SSA model ice velocity in the X direction",
                            "m s-1", NULL); CHKERRQ(ierr);
  ierr = vubarSSA.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  ierr = vvbarSSA.create(grid, "vvbarSSA", true);
  ierr = vvbarSSA.set_attrs("internal_restart", "SSA model ice velocity in the Y direction",
                            "m s-1", NULL); CHKERRQ(ierr);
  ierr = vvbarSSA.set_glaciological_units("m year-1"); CHKERRQ(ierr);


  createVecs_done = PETSC_TRUE;
  return 0;
}


//! De-allocate all Vecs defined in IceModel.
/*! 
  Undoes the actions of createVecs().
 */
PetscErrorCode IceModel::destroyVecs() {
  PetscErrorCode ierr;

  ierr = bedDefCleanup(); CHKERRQ(ierr);

  ierr = u3.destroy(); CHKERRQ(ierr);
  ierr = v3.destroy(); CHKERRQ(ierr);
  ierr = w3.destroy(); CHKERRQ(ierr);
  ierr = Sigma3.destroy(); CHKERRQ(ierr);
  ierr = T3.destroy(); CHKERRQ(ierr);
  ierr = tau3.destroy(); CHKERRQ(ierr);

  ierr = Tb3.destroy(); CHKERRQ(ierr);

  ierr = vh.destroy(); CHKERRQ(ierr);
  ierr = vH.destroy(); CHKERRQ(ierr);
  ierr = vbed.destroy(); CHKERRQ(ierr);
  ierr = vMask.destroy(); CHKERRQ(ierr);
  ierr = vGhf.destroy(); CHKERRQ(ierr);
  ierr = vubar.destroy(); CHKERRQ(ierr);
  ierr = vvbar.destroy(); CHKERRQ(ierr);
  ierr = vub.destroy(); CHKERRQ(ierr);
  ierr = vvb.destroy(); CHKERRQ(ierr);

  ierr = vRb.destroy(); CHKERRQ(ierr);

  ierr = vHmelt.destroy(); CHKERRQ(ierr);
  ierr = vbasalMeltRate.destroy(); CHKERRQ(ierr);
  ierr = vuplift.destroy(); CHKERRQ(ierr);
  ierr = vdHdt.destroy(); CHKERRQ(ierr);

  ierr = vtauc.destroy(); CHKERRQ(ierr);
  ierr = vtillphi.destroy(); CHKERRQ(ierr);

  ierr = vLongitude.destroy(); CHKERRQ(ierr);
  ierr = vLatitude.destroy(); CHKERRQ(ierr);

  ierr = vuvbar[0].destroy(); CHKERRQ(ierr);
  ierr = vuvbar[1].destroy(); CHKERRQ(ierr);

  for (int j = 0; j < nWork2d; j++) {
    ierr = vWork2d[j].destroy(); CHKERRQ(ierr);
  }

  ierr = vubarSSA.destroy(); CHKERRQ(ierr);
  ierr = vvbarSSA.destroy(); CHKERRQ(ierr);

  ierr = Tnew3.destroy(); CHKERRQ(ierr);
  ierr = taunew3.destroy(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].destroy(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].destroy(); CHKERRQ(ierr);
  ierr = Istag3[0].destroy(); CHKERRQ(ierr);
  ierr = Istag3[1].destroy(); CHKERRQ(ierr);

  ierr = VecDestroy(g2); CHKERRQ(ierr);

  ierr = KSPDestroy(SSAKSP); CHKERRQ(ierr);
  ierr = MatDestroy(SSAStiffnessMatrix); CHKERRQ(ierr);
  ierr = VecDestroy(SSAX); CHKERRQ(ierr);
  ierr = VecDestroy(SSARHS); CHKERRQ(ierr);
  ierr = VecDestroy(SSAXLocal); CHKERRQ(ierr);
  ierr = VecScatterDestroy(SSAScatterGlobalToLocal); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::attachAtmospherePCC(PISMAtmosphereCoupler &aPCC) {
  atmosPCC = &aPCC;
  return 0;
}


PetscErrorCode IceModel::attachOceanPCC(PISMOceanCoupler &oPCC) {
  oceanPCC = &oPCC;
  return 0;
}


void IceModel::setMaxTimeStepYears(PetscScalar y) {
  maxdt = y * secpera;
  doAdaptTimeStep = PETSC_TRUE;
}


void IceModel::setAdaptTimeStepRatio(PetscScalar c) {
  adaptTimeStepRatio = c;
}


PetscErrorCode IceModel::setStartYear(PetscScalar y0) {
  startYear = y0;
  return 0;
}


PetscErrorCode IceModel::setEndYear(PetscScalar ye) {    
  if (ye < startYear)   {
    PetscPrintf(grid.com,
		"PISM ERROR: ye < startYear. PISM cannot run backward in time.\n");
    PetscEnd();
  }
  endYear = ye;
  return 0;
}


void  IceModel::setInitialAgeYears(PetscScalar d) {
  tau3.set(d*secpera);
}


void IceModel::setAllGMaxVelocities(PetscScalar uvw_for_cfl) {
  gmaxu=uvw_for_cfl;
  gmaxv=uvw_for_cfl;
  gmaxw=uvw_for_cfl;
}


void IceModel::setConstantNuHForSSA(PetscScalar nuH) {
  useConstantNuHForSSA = PETSC_TRUE;
  ssaStrengthExtend.set_notional_strength(nuH);
}


PetscErrorCode IceModel::setExecName(const char *my_executable_short_name) {
  strcpy(executable_short_name, my_executable_short_name);
  return 0;
}

//! Do the time-stepping for an evolution run.
/*! 
This procedure is the main time-stepping loop.  The following actions are taken on each pass 
through the loop:
\li the yield stress for the plastic till model is updated (if appropriate)
\li the positive degree day model is invoked to compute the surface mass balance (if appropriate)
\li a step of the bed deformation model is taken (if appropriate)
\li the velocity field is updated; in some cases the whole three-dimensional field is updated 
    and in some cases just the vertically-averaged horizontal velocity is updated; see velocity()
\li the time step is determined according to a variety of stability criteria; 
    see determineTimeStep()
\li the temperature field is updated according to the conservation of energy model based 
    (especially) on the new velocity field; see temperatureAgeStep()
\li the thickness of the ice is updated according to the mass conservation model; see
    massContExplicitStep()
\li there is various reporting to the user on the current state; see summary() and updateViewers()

Note that at the beginning and ends of each pass throught the loop there is a chance for 
derived classes to do extra work.  See additionalAtStartTimestep() and additionalAtEndTimestep().
 */
PetscErrorCode IceModel::run() {
  PetscErrorCode  ierr;

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

  ierr = summaryPrintLine(PETSC_TRUE,doTemp, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); CHKERRQ(ierr);
  adaptReasonFlag = '$'; // no reason for no timestep
  skipCountDown = 0;
  ierr = summary(doTemp,reportHomolTemps); CHKERRQ(ierr);  // report starting state
  dtTempAge = 0.0;

  // main loop for time evolution
  for (PetscScalar year = startYear; year < endYear; year += dt/secpera) {
    write_snapshot();

    ierr = verbPrintf(2,grid.com, " "); CHKERRQ(ierr);
    dt_force = -1.0;
    maxdt_temporary = -1.0;
    ierr = additionalAtStartTimestep(); CHKERRQ(ierr);  // might set dt_force,maxdt_temporary
    
PetscLogEventBegin(beddefEVENT,0,0,0,0);

    // compute bed deformation, which only depends on current thickness and bed elevation
    if (doBedDef == PETSC_TRUE) {
      ierr = bedDefStepIfNeeded(); CHKERRQ(ierr); // prints "b" or "$" as appropriate
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }
    
PetscLogEventEnd(beddefEVENT,0,0,0,0);

    // update basal till yield stress if appropriate; will modify and communicate mask
    if (doPlasticTill == PETSC_TRUE) {
      ierr = updateYieldStressFromHmelt();  CHKERRQ(ierr);
      ierr = verbPrintf(2,grid.com, "y"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }

    // do vertical grid blow-out checking just before first 3D field gets
    // FIXME:  this version just ends; task #4218: expand grid upward instead of ending!
    ierr = thicknessTooLargeCheck(); CHKERRQ(ierr);

    // always do SIA velocity calculation; only update SSA and 
    //   only update velocities at depth if suggested by temp and age
    //   stability criterion; note *lots* of communication is avoided by 
    //   skipping SSA (and temp/age)
    bool updateAtDepth = (skipCountDown == 0);
    ierr = velocity(updateAtDepth); CHKERRQ(ierr);  // event logging in here
    ierr = verbPrintf(2,grid.com, updateAtDepth ? "v" : "V" ); CHKERRQ(ierr);
    
    // adapt time step using velocities and diffusivity, ..., just computed
    bool useCFLforTempAgeEqntoGetTimestep = (doTemp == PETSC_TRUE);
    ierr = determineTimeStep(useCFLforTempAgeEqntoGetTimestep); CHKERRQ(ierr);
    dtTempAge += dt;
    grid.year += dt / secpera;  // adopt it
    // IceModel::dt,dtTempAge,grid.year are now set correctly according to
    //    mass-continuity-eqn-diffusivity criteria, horizontal CFL criteria, and other 
    //    criteria from derived class additionalAtStartTimestep(), and from 
    //    "-skip" mechanism

    // ierr = PetscPrintf(PETSC_COMM_SELF,
    //           "\n[rank=%d, year=%f, dt=%f, startYear=%f, endYear=%f]",
    //           grid.rank, grid.year, dt/secpera, startYear, endYear);
    //        CHKERRQ(ierr);

PetscLogEventBegin(tempEVENT,0,0,0,0);
    
    bool tempAgeStep = (updateAtDepth && (doTemp == PETSC_TRUE));
    if (tempAgeStep) { // do temperature and age
      ierr = temperatureAgeStep(); CHKERRQ(ierr);
      dtTempAge = 0.0;
      if (updateHmelt == PETSC_TRUE) {
        ierr = diffuseHmelt(); CHKERRQ(ierr);
      }
      ierr = verbPrintf(2,grid.com, "at"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$$"); CHKERRQ(ierr);
    }

PetscLogEventEnd(tempEVENT,0,0,0,0);
PetscLogEventBegin(massbalEVENT,0,0,0,0);

    if (doMassConserve == PETSC_TRUE) {
      ierr = massContExplicitStep(); CHKERRQ(ierr); // update H
      ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr); // update h and mask
      if ((doSkip == PETSC_TRUE) && (skipCountDown > 0))
        skipCountDown--;
      ierr = verbPrintf(2,grid.com, "h"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }

PetscLogEventEnd(massbalEVENT,0,0,0,0);
    
    ierr = additionalAtEndTimestep(); CHKERRQ(ierr);

    // end the flag line and report a summary
    ierr = verbPrintf(2,grid.com, " %d%c  +%6.5f\n", skipCountDown, adaptReasonFlag,
                      dt / secpera); CHKERRQ(ierr);
    ierr = summary(tempAgeStep,reportHomolTemps); CHKERRQ(ierr);

    ierr = updateViewers(); CHKERRQ(ierr);

    if (endOfTimeStepHook() != 0) break;
  }

  return 0;
}


//! Calls the necessary routines to do a diagnostic calculation of velocity.
/*! 
This important routine can be replaced by derived classes; it is \c virtual.

This procedure has no loop but the following actions are taken:
\li the yield stress for the plastic till model is updated (if appropriate)
\li the velocity field is updated; in some cases the whole three-dimensional field is updated 
    and in some cases just the vertically-averaged horizontal velocity is updated; see velocity()
\li there is various reporting to the user on the current state; see summary() and updateViewers()
 */
PetscErrorCode IceModel::diagnosticRun() {
  PetscErrorCode  ierr;

  // print out some stats about input state
  ierr = summaryPrintLine(PETSC_TRUE,PETSC_TRUE, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
           CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep
  skipCountDown = 0;

  // update basal till yield stress if appropriate; will modify and communicate mask
  if (doPlasticTill == PETSC_TRUE) {
    ierr = updateYieldStressFromHmelt();  CHKERRQ(ierr);
  }

  ierr = velocity(true); CHKERRQ(ierr);  // compute velocities (at depth)

  ierr = summary(true,true); CHKERRQ(ierr);
  
  // update viewers and pause for a chance to view
  ierr = updateViewers(); CHKERRQ(ierr);
  PetscInt    pause_time = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, PETSC_NULL); CHKERRQ(ierr);
  if (pause_time > 0) {
    ierr = verbPrintf(2,grid.com,"pausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);
  }
  return 0;
}

//! Set default values of grid parameters.
/*!
  Derived classes (IceCompModel, for example) reimplement this to change the
  grid initialization when no -i option is set.
 */
PetscErrorCode IceModel::set_grid_defaults() {
  PetscErrorCode ierr;
  PetscTruth Mx_set, My_set, Mz_set, Lz_set, boot_from_set;
  char filename[PETSC_MAX_PATH_LEN];
  grid_info gi;

  // Get the bootstrapping file name:
  ierr = check_old_option_and_stop("-bif", "-boot_from"); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-boot_from",
			       filename, PETSC_MAX_PATH_LEN, &boot_from_set); CHKERRQ(ierr);

  if (!boot_from_set) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: Please specify an input file using -i or -boot_from.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  // Use a bootstrapping file to set some grid parameters (they can be
  // overridden later, in IceModel::set_grid_from_options()).

  // Determine the grid extent from a bootstrapping file:
  NCTool nc(&grid);
  bool x_dim_exists, y_dim_exists, t_exists;
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

  ierr = nc.find_dimension("x", NULL, x_dim_exists); CHKERRQ(ierr);
  ierr = nc.find_dimension("y", NULL, y_dim_exists); CHKERRQ(ierr);
  ierr = nc.find_variable("t", NULL, NULL, t_exists); CHKERRQ(ierr);
  ierr = nc.get_grid_info(gi);
  ierr = nc.close(); CHKERRQ(ierr);

  // if the horizontal dimensions are absent then we can not proceed
  if (!x_dim_exists) {
    ierr = PetscPrintf(grid.com,"bootstrapping file '%s' has no horizontal dimension 'x'\n",filename);
    CHKERRQ(ierr);
    PetscEnd();
  }
  if (!y_dim_exists) {
    ierr = PetscPrintf(grid.com,"bootstrapping file '%s' has no horizontal dimension 'y'\n",filename);
    CHKERRQ(ierr);
    PetscEnd();
  }

  // Set the grid center and horizontal extent:
  grid.x0 = gi.x0;
  grid.y0 = gi.y0;
  grid.Lx = gi.Lx;
  grid.Ly = gi.Ly;

  if (t_exists) {
    grid.year = gi.time / secpera; // set year from read-in time variable
    ierr = verbPrintf(2, grid.com, 
		      "  time t = %5.4f years found; setting current year\n",
		      grid.year); CHKERRQ(ierr);
  } else {
    grid.year = 0.0;
    ierr = verbPrintf(2, grid.com, 
		      "  time dimension was not found; setting current year to t = 0.0 years\n",
		      grid.year); CHKERRQ(ierr);
  }

  // Grid dimensions and its vertical extent should not be deduced from a
  // bootstrapping file, so we check if these options are set and stop if they
  // are not.
  // Note that here interpreting "-Mx 0" as "-Mx was not set" is OK.
  ierr = check_option("-Mx", Mx_set); CHKERRQ(ierr);
  ierr = check_option("-My", My_set); CHKERRQ(ierr);
  ierr = check_option("-Mz", Mz_set); CHKERRQ(ierr);
  ierr = check_option("-Lz", Lz_set); CHKERRQ(ierr);
  if ( !(Mx_set && My_set && Mz_set && Lz_set) ) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: All of -boot_from, -Mx, -My, -Mz, -Lz, are required for bootstrapping.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  return 0;
}

//! Initalizes the grid from options.
/*! 
  Reads all of -Mx, -My, -Mz, -Mbz, -Lx, -Ly, -Lz, -quadZ and -chebZ. Sets
  corresponding grid parameters.
 */
PetscErrorCode IceModel::set_grid_from_options() {
  PetscErrorCode ierr;
  PetscTruth Mx_set, My_set, Mz_set, Mbz_set, Lx_set, Ly_set, Lz_set,
    quadZ_set, chebZ_set;
  PetscScalar x_scale, y_scale, z_scale;
  int Mx, My, Mz, Mbz;

  // Process the options:
  ierr = PetscOptionsBegin(grid.com, PETSC_NULL,
			   "Options setting the computational grid extent and dimensions",
			   PETSC_NULL); CHKERRQ(ierr);

  // Read -Lx and -Ly. Note the transpose!
  ierr = PetscOptionsScalar("-Lx", "Half of the grid extent in the X direction, in km", "",
			    y_scale, &y_scale, &Ly_set); CHKERRQ(ierr);
  ierr = PetscOptionsScalar("-Ly", "Half of the grid extent in the Y direction, in km", "",
			    x_scale, &x_scale, &Lx_set); CHKERRQ(ierr);
  // Vertical extent (in the ice):
  ierr = PetscOptionsScalar("-Lz", "Grid extent in the Z (vertical) direction in the ice", "",
			    z_scale, &z_scale, &Lz_set); CHKERRQ(ierr);

  // Read -Mx, -My, -Mz and -Mbz. Note the transpose!
  ierr = PetscOptionsInt("-Mx", "Number of grid points in the X direction", "",
			 grid.My, &My, &My_set); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-My", "Number of grid points in the Y direction", "",
			 grid.Mx, &Mx, &Mx_set); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-Mz", "Number of grid points in the Z (vertical) direction in the ice", "",
			 grid.Mz, &Mz, &Mz_set); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-Mbz", "Number of grid points in the Z (vertical) direction in the bedrock", PETSC_NULL,
			 grid.Mbz, &Mbz, &Mbz_set); CHKERRQ(ierr);

  // Determine the vertical grid spacing in the ice:
  ierr = PetscOptionsName("-quadZ", "Chooses the quadratic vertical grid spacing",
			  PETSC_NULL, &quadZ_set); CHKERRQ(ierr);
  ierr = PetscOptionsName("-chebZ", "Chooses the Chebyshev vertical grid spacing",
			  PETSC_NULL, &chebZ_set); CHKERRQ(ierr);

  // Only one of -quadZ and -chebZ is allowed.
  if (quadZ_set && chebZ_set) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: at most one of -quadZ and -chebZ is allowed.\n"); CHKERRQ(ierr);
    PetscEnd();
  }

  // Done with the options.
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // Use the information obtained above:
  if (Lx_set)    grid.Lx = x_scale * 1000.0; // convert to meters
  if (Ly_set)    grid.Ly = y_scale * 1000.0; // convert to meters
  if (Lz_set)    grid.Lz = z_scale;	     // in meters already
  if (Mx_set)    grid.Mx = Mx;
  if (My_set)    grid.My = My;
  if (Mz_set)    grid.Mz = Mz;
  if (Mbz_set)   grid.Mbz = Mbz;
  if (quadZ_set) grid.vertical_spacing = QUADRATIC;
  if (chebZ_set) grid.vertical_spacing = CHEBYSHEV;

  ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
  ierr = grid.compute_vertical_levels();    CHKERRQ(ierr);

  // At this point all the fields except for da2, xs, xm, ys, ym should be
  // filled. We're ready to call grid.createDA().
  return 0;
}

//! Manage the initialization of the IceModel object.
/*!
The IceModel initialization sequence is this:
  
   1) Initialize the computational grid.

   2) Process the options.

   3) Memory allocation.

   4) Initialize IceType and (possibly) other physics.

   5) Initialize PDD and forcing.

   6) Fill the model state variables (from a PISM output file, from a
   bootstrapping file using some modeling choices or using formulas).

   7) Regrid.

   8) Report grid parameters.

   9) Allocate internal objects: SSA tools and work vectors. Some tasks in the
   next (tenth) item (bed deformation setup, for example) might need this.

   10) Miscellaneous stuff: update surface elevation and mask, set up the bed
   deformation model, initialize the basal till model, initialize snapshots.

Please see the documenting comments of the functions called below to find 
explanations of their intended uses.  See the flow-chart
doc/initialization_sequence.png for a graphical illustration of the process.
 */
PetscErrorCode IceModel::init() {
  PetscErrorCode ierr;

  // 1) Initialize the computational grid:
  ierr = grid_setup(); CHKERRQ(ierr);

  // 2) Process the options:
  ierr = setFromOptions(); CHKERRQ(ierr);

  // 3) Memory allocation:
  ierr = createVecs(); CHKERRQ(ierr);

  // 4) Initialize the IceType and (possibly) other physics.
  ierr = init_physics(); CHKERRQ(ierr);

  // 5) Initialize atmosphere and ocean couplers:
  ierr = init_couplers(); CHKERRQ(ierr);

  // 6) Fill the model state variables (from a PISM output file, from a
  // bootstrapping file using some modeling choices or using formulas).
  ierr = model_state_setup(); CHKERRQ(ierr);

  // 7) Regrid:
  ierr = regrid(); CHKERRQ(ierr);

  // 8) Report grid parameters:
  ierr = report_grid_parameters(); CHKERRQ(ierr);

  // 9) Allocate SSA tools and work vectors:
  ierr = allocate_internal_objects(); CHKERRQ(ierr);

  // 10) Miscellaneous stuff: update surface elevation and mask, set up the bed
  // deformation model, initialize the basal till model, initialize snapshots.
  // This has to happen *after* regridding.
  ierr = misc_setup();

  return 0; 
}

//! Sets up the computational grid.
/*!
  There are two cases here:

  1) Initializing from a PISM ouput file, in which case all the options
  influencing the grid (currently: -Mx, -My, -Mz, -Mbz, -Lx, -Ly, -Lz, -quadZ,
  -chebZ) are ignored.

  2) Initializing using defaults, command-line options and (possibly) a
  bootstrapping file. Derived classes requiring special grid setup should
  reimplement IceGrid::set_grid_from_options().

  No memory allocation should happen here.
 */
PetscErrorCode IceModel::grid_setup() {
  PetscErrorCode ierr;
  PetscTruth i_set;
  char filename[PETSC_MAX_PATH_LEN];

  ierr = verbPrintf(3, grid.com,
		    "Setting up the computational grid...\n"); CHKERRQ(ierr);

  // Check if we are initializing from a PISM output file:
  ierr = check_old_option_and_stop("-if", "-i"); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-i",
			       filename, PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);

  if (i_set) {
    NCTool nc(&grid);
    char *tmp;
    int length;

    // Get the 'source' global attribute to check if we are given a PISM output
    // file:
    ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
    ierr = nc.get_att_text(NC_GLOBAL, "source", &length, &tmp);
    ierr = nc.close();

    // If it's missing, print a warning
    if (tmp == NULL) {
      ierr = verbPrintf(1, grid.com,
			"PISM WARNING: file '%s' does not have the 'source' global attribute.\n"
			"     If '%s' is a PISM output file, please run the following to get rid of this warning:\n"
			"     ncatted -a source,global,c,c,PISM %s\n",
			filename, filename, filename); CHKERRQ(ierr);
    } else if (strstr(tmp, "PISM") == NULL) {
      // If the 'source' attribute does not contain the string "PISM", then print
      // a message and stop:
      ierr = verbPrintf(1, grid.com,
			"PISM WARNING: '%s' does not seem to be a PISM output file.\n"
			"     If it is, please make sure that the 'source' global attribute contains the string \"PISM\".\n",
			filename); CHKERRQ(ierr);
    }
    delete[] tmp;

    ierr = nc.get_grid(filename);   CHKERRQ(ierr);

    // These options are ignored because we're getting *all* the grid
    // parameters from a file.
    ierr = ignore_option("-Mx");    CHKERRQ(ierr);
    ierr = ignore_option("-My");    CHKERRQ(ierr);
    ierr = ignore_option("-Mz");    CHKERRQ(ierr);
    ierr = ignore_option("-Mbz");   CHKERRQ(ierr);
    ierr = ignore_option("-Lx");    CHKERRQ(ierr);
    ierr = ignore_option("-Ly");    CHKERRQ(ierr);
    ierr = ignore_option("-Lz");    CHKERRQ(ierr);
    ierr = ignore_option("-chebZ"); CHKERRQ(ierr);
    ierr = ignore_option("-quadZ"); CHKERRQ(ierr);
  } else {
    ierr = set_grid_defaults(); CHKERRQ(ierr);
    ierr = set_grid_from_options(); CHKERRQ(ierr);
  }

  ierr = grid.createDA(); CHKERRQ(ierr);
  
  return 0;
}

//! Sets the starting values of model state variables.
/*!
  There are two cases:
  
  1) Initializing from a PISM output file.

  2) Setting the values using command-line options only (verification and
  simplified geometry runs, for example) or from a bootstrapping file, using
  heuristics to fill in missing and 3D fields.

  This function is called after all the memory allocation is done and all the
  physical parameters are set.
 */
PetscErrorCode IceModel::model_state_setup() {
  PetscErrorCode ierr;
  PetscTruth i_set;
  char filename[PETSC_MAX_PATH_LEN];
  
  // Check if we are initializing from a PISM output file:
  ierr = PetscOptionsGetString(PETSC_NULL, "-i",
			       filename, PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);

  if (i_set) {
    ierr = initFromFile(filename); CHKERRQ(ierr);
  } else {
    ierr = set_vars_from_options(); CHKERRQ(ierr);
  }

  return 0;
}

//! Sets starting values of model state variables using command-line options.
/*!
  Sets starting values of model state variables using command-line options and
  (possibly) a bootstrapping file.

  In the base class there is only one case: bootstrapping.
 */
PetscErrorCode IceModel::set_vars_from_options() {
  PetscErrorCode ierr;
  PetscTruth boot_from_set;
  char filename[PETSC_MAX_PATH_LEN];

  ierr = verbPrintf(3, grid.com,
		    "Setting initial values of model state variables...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-boot_from",
			       filename, PETSC_MAX_PATH_LEN, &boot_from_set); CHKERRQ(ierr);
  
  if (boot_from_set) {
    ierr = bootstrapFromFile(filename); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(grid.com, "PISM ERROR: No input file specified.\n"); CHKERRQ(ierr);
    PetscEnd();
  }
  
  return 0;
}

//! Initialize some physical parameters.
/*!
  This is the place for all non-trivial initalization of physical parameters
  (non-trivial meaning requiring more than just setting a value of a
  parameter).

  This method is called after memory allocation but before filling any of
  IceModelVecs.

  Rationale: all the physical parameters should be initialized before setting
  up the coupling or filling model-state variables.

  In the base class we just initialize the IceType and the shelf extension.

  Also, this is the good place for setting parameters that a user should not be
  able to override using a command-line option.
 */
PetscErrorCode IceModel::init_physics() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com,
		    "Initializing IceType and shelfExtension...\n"); CHKERRQ(ierr);

  ierr = iceFactory.setFromOptions(); CHKERRQ(ierr);
  // Initialize the IceType object:
  if (ice == PETSC_NULL) {
    ierr = iceFactory.create(&ice); CHKERRQ(ierr);
    ierr = ice->setFromOptions(); CHKERRQ(ierr); // Set options specific to this particular ice type
  }

  return 0;
}

//! Miscellaneous initialization tasks plus tasks that need the fields that can come from regridding.
PetscErrorCode IceModel::misc_setup() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com, "Finishing initialization...\n"); CHKERRQ(ierr);

  ierr = init_snapshots_from_options(); CHKERRQ(ierr);
  ierr = stampHistoryCommand(); CHKERRQ(ierr);
  ierr = createViewers(); CHKERRQ(ierr);

  // consistency of geometry after initialization;
  ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);

  // allocate and setup bed deformation model
  ierr = bedDefSetup(); CHKERRQ(ierr);

  // init basal till model, possibly inverting for phi, if desired;
  //   reads options "-topg_to_phi phi_min,phi_max,phi_ocean,topg_min,topg_max"
  //   or "-surf_vel_to_phi foo.nc";
  //   initializes PlasticBasalType* basal; sets fields vtauc, vtillphi
  ierr = initBasalTillModel(); CHKERRQ(ierr);
  
  return 0;
}

//! Initializes atmosphere and ocean couplers.
PetscErrorCode IceModel::init_couplers() {
  PetscErrorCode ierr;

  // so that we can let atmosPCC, oceanPCC know about these fields in IceModel state
  info_atmoscoupler.lat = &vLatitude;
  info_atmoscoupler.lon = &vLongitude;  
  info_atmoscoupler.mask = &vMask;
  info_atmoscoupler.surfelev = &vh;

  info_oceancoupler.lat = &vLatitude;
  info_oceancoupler.lon = &vLongitude;  
  info_oceancoupler.mask = &vMask;
  info_oceancoupler.thk = &vH;

  ierr = verbPrintf(3, grid.com,
		    "Initializing atmosphere and ocean couplers...\n"); CHKERRQ(ierr);

  if (atmosPCC != PETSC_NULL) {
    ierr = atmosPCC->initFromOptions(&grid); CHKERRQ(ierr);
  } else {  SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");  }
 
 if (oceanPCC != PETSC_NULL) {
    if (isDrySimulation == PETSC_TRUE) {  oceanPCC->reportInitializationToStdOut = false;  }
    ierr = oceanPCC->initFromOptions(&grid); CHKERRQ(ierr);
  } else {  SETERRQ(2,"PISM ERROR: oceanPCC == PETSC_NULL");  }


  return 0;
}

//! Allocates SSA tools and work vectors.
PetscErrorCode IceModel::allocate_internal_objects() {
  PetscErrorCode ierr;

  // a global Vec is needed for things like viewers and comm to proc zero
  ierr = DACreateGlobalVector(grid.da2, &g2); CHKERRQ(ierr);

  // setup (classical) SSA tools
  const PetscInt M = 2 * grid.Mx * grid.My;
  ierr = MatCreateMPIAIJ(grid.com, PETSC_DECIDE, PETSC_DECIDE, M, M,
                         13, PETSC_NULL, 13, PETSC_NULL,
                         &SSAStiffnessMatrix); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com, PETSC_DECIDE, M, &SSAX); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &SSARHS); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, M, &SSAXLocal);
  ierr = VecScatterCreate(SSAX, PETSC_NULL, SSAXLocal, PETSC_NULL,
                          &SSAScatterGlobalToLocal); CHKERRQ(ierr);
  ierr = KSPCreate(grid.com, &SSAKSP); CHKERRQ(ierr);


  // various internal quantities
  // 2d work vectors
  for (int j = 0; j < nWork2d; j++) {
    ierr = vWork2d[j].create(grid, "a_work_vector", true); CHKERRQ(ierr);
  }

  // 3d dedicated work vectors
  ierr = Tnew3.createSameDA(T3,grid,"temp_new",false); CHKERRQ(ierr);
  ierr = Tnew3.set_attrs("internal", "ice temperature; temporary during update",
                         "K", NULL); CHKERRQ(ierr);
  ierr = taunew3.createSameDA(tau3,grid,"age_new",false); CHKERRQ(ierr);
  ierr = taunew3.set_attrs("internal", "age of ice; temporary during update",
                           "s", NULL); CHKERRQ(ierr);
  ierr = Sigmastag3[0].create(grid,"Sigma_stagx",true); CHKERRQ(ierr);
  ierr = Sigmastag3[0].set_attrs("internal",
             "rate of strain heating; on staggered grid offset in X direction",
	     "J s-1 m-3", NULL); CHKERRQ(ierr);
  ierr = Sigmastag3[1].create(grid,"Sigma_stagy",true); CHKERRQ(ierr);
  ierr = Sigmastag3[1].set_attrs("internal",
             "rate of strain heating; on staggered grid offset in Y direction",
	     "J s-1 m-3", NULL); CHKERRQ(ierr);
  ierr = Istag3[0].create(grid,"I_stagx",true); CHKERRQ(ierr);
  ierr = Istag3[0].set_attrs("internal",NULL,NULL,NULL); CHKERRQ(ierr);
  ierr = Istag3[1].create(grid,"I_stagy",true); CHKERRQ(ierr);
  ierr = Istag3[1].set_attrs("internal",NULL,NULL,NULL); CHKERRQ(ierr);

  return 0;
}

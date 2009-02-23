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
#include "ssa/pismssa.hh"

IceModel::IceModel(IceGrid &g): grid(g), iceFactory(grid.com,NULL), ice(NULL) {
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
  top0ctx_created = PETSC_FALSE;
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

  save_snapshots = false;
  dvoldt = gdHdtav = 0;

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
/*! Initialization of an IceModel is confusing.  Here is a description of the intended order:
	\li 1. The constructor for IceModel.  Note IceModel has a member "grid", of class IceGrid. 
	   The IceGrid constructor sets 
	   defaults for (grid.)Mx,My,Mz,Mbz,Lx,Ly,Lz,Lbz,dx,dy,dz,year.
	\li [1.5] derivedClass::setFromOptions() to get options special to derived class
	\li 2. setFromOptions() to get all options *including* Mx,My,Mz,Mbz
	\li [2.5] initFromFile() which reads Mx,My,Mz,Mbz from file and overwrites previous; if 
	   this represents a change the user is warned
	\li 3. createDA(), which uses only Mx,My,Mz,Mbz
	\li 4. createVecs() uses DA to create/allocate Vecs
	\li [4.5] derivedClass:createVecs() to create/allocate Vecs special to derived class
	\li 5. afterInitHook() which changes Lx,Ly,Lz if set by user

Note driver programs call only setFromOptions() and initFromOptions() (for IceModel 
or derived class).

Note IceModel::setFromOptions() should be called at the end of derivedClass:setFromOptions().

Note 2.5, 3, and 4 are called from initFromFile() in IceModel.

Note 3 and 4 are called from initFromOptions() in some derived classes (e.g. IceCompModel) 
in cases where initFromFile() is not called.

Note step 2.5 is skipped when bootstrapping (-boot_from and bootstrapFromFile()) or in
those derived classes which can start with no input files, e.g. IceCompModel and IceEISModel.
That is, 2.5 is only done when starting from a saved model state.
*/
PetscErrorCode IceModel::createVecs() {
  PetscErrorCode ierr;

  if (createVecs_done == PETSC_TRUE) {
    ierr = destroyVecs(); CHKERRQ(ierr);
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
	        	  "J s-1 m-3", NULL); CHKERRQ(ierr);

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
  ierr = tau3.set_glaciological_units("year");
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
  ierr = vGhf.set_attrs("climate_steady", "upward geothermal flux at bedrock surface",
			"W m-2", NULL); CHKERRQ(ierr);
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
  ierr = vtauc.set_attrs("diagnostic", 
             "yield stress for basal till (plastic or pseudo-plastic model)",
	     "Pa", NULL); CHKERRQ(ierr);

  // bedrock uplift rate
  ierr = vuplift.create(grid, "dbdt", true); CHKERRQ(ierr);
  ierr = vuplift.set_attrs("model_state", "bedrock uplift rate",
			   "m s-1", "tendency_of_bedrock_altitude"); CHKERRQ(ierr);
  ierr = vuplift.set_glaciological_units("m year-1");

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
  ierr = vtillphi.set_attrs("climate_steady", "friction angle for till under grounded ice sheet",
			    "degrees", NULL); CHKERRQ(ierr);

  // longitude
  ierr = vLongitude.create(grid, "lon", true); CHKERRQ(ierr);
  ierr = vLongitude.set_attrs("mapping", "longitude", "degrees_east", "longitude"); CHKERRQ(ierr);

  // latitude
  ierr = vLatitude.create(grid, "lat", true); CHKERRQ(ierr);
  ierr = vLatitude.set_attrs("mapping", "latitude", "degrees_north", "latitude"); CHKERRQ(ierr);

  // u bar and v bar on staggered grid
  ierr = vuvbar[0].create(grid, "vuvbar[0]", true); CHKERRQ(ierr);
  ierr = vuvbar[0].set_attrs("internal", 
            "vertically averaged ice velocity, on staggered grid offset in X direction, from SIA, in the X direction",
	    "m s-1", NULL); CHKERRQ(ierr);
  ierr = vuvbar[1].create(grid, "vuvbar[1]", true); CHKERRQ(ierr);
  ierr = vuvbar[1].set_attrs("internal", 
            "vertically averaged ice velocity, on staggered grid offset in Y direction, from SIA, in the Y direction",
	    "m s-1", NULL); CHKERRQ(ierr);

  // work vectors
  for (int j = 0; j < nWork2d; j++) {
    ierr = vWork2d[j].create(grid, "a_work_vector", true); CHKERRQ(ierr);
  }

  // initial guesses of SSA velocities
  ierr = vubarSSA.create(grid, "vubarSSA", true);
  ierr = vubarSSA.set_attrs("internal_restart", "SSA model ice velocity in the X direction",
                            "m s-1", NULL); CHKERRQ(ierr);
  ierr = vubarSSA.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  ierr = vvbarSSA.create(grid, "vvbarSSA", true);
  ierr = vvbarSSA.set_attrs("internal_restart", "SSA model ice velocity in the Y direction",
                            "m s-1", NULL); CHKERRQ(ierr);
  ierr = vvbarSSA.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  // various internal quantities
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

  ierr = DACreateGlobalVector(grid.da2, &g2); CHKERRQ(ierr);

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

  // so that we can let atmosPCC know about these fields in IceModel state
  info_atmoscoupler.lat = &vLatitude;
  info_atmoscoupler.lon = &vLongitude;  
  info_atmoscoupler.mask = &vMask;
  info_atmoscoupler.surfelev = &vh;

  // so that we can let oceanPCC know about these fields in IceModel state
  info_oceancoupler.lat = &vLatitude;
  info_oceancoupler.lon = &vLongitude;  
  info_oceancoupler.mask = &vMask;
  info_oceancoupler.thk = &vH;

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
  constantNuHForSSA = nuH;
}


PetscErrorCode IceModel::setExecName(const char *my_executable_short_name) {
  strcpy(executable_short_name, my_executable_short_name);
  return 0;
}


PetscTruth IceModel::isInitialized() const {
  return initialized_p;
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

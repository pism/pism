// Copyright (C) 2004-2008 Jed Brown, Ed Bueler and Constantine Khroulev
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

// following numerical values have some significance; see updateSurfaceElevationAndMask() below
const int IceModel::MASK_SHEET = 1;
const int IceModel::MASK_DRAGGING = 2;
const int IceModel::MASK_FLOATING = 3;
// (modMask(mask[i][j]) == MASK_FLOATING) is criteria for floating; ..._OCEAN0 only used if -ocean_kill 
const int IceModel::MASK_FLOATING_OCEAN0 = 7;


IceModel::IceModel(IceGrid &g, IceType *i): grid(g), ice(i) {
  PetscErrorCode ierr;

  bootstrapLIC = PETSC_NULL;
  
  history_size = TEMPORARY_STRING_LENGTH;
  history = new char[history_size];

  have_ssa_velocities = false;

  pism_signal = 0;
  signal(SIGTERM, pism_signal_handler);
  signal(SIGUSR1, pism_signal_handler);

  createBasal_done = PETSC_FALSE;
  top0ctx_created = PETSC_FALSE;
  createVecs_done = PETSC_FALSE;

  for (PetscInt nn = 0; nn < tnN; nn++) {
    runtimeViewers[nn] = PETSC_NULL;
  }
  createViewers_done = PETSC_FALSE;

  pddStuffCreated = PETSC_FALSE;
  pddRandStuffCreated = PETSC_FALSE;

  dTforcing = PETSC_NULL;
  dSLforcing = PETSC_NULL;

  vmonthlyTs = PETSC_NULL;
  
  ierr = setDefaults();  // lots of parameters and flags set here
  if (ierr != 0) {
    verbPrintf(1,grid.com, "Error setting defaults.\n");
    PetscEnd();
  }

  ierr = getFlowLawNumber(flowLawNumber, flowLawNumber); //CHKERRQ(ierr);
  if (flowLawNumber == 4)   flowLawUsesGrainSize = PETSC_TRUE;
  else                      flowLawUsesGrainSize = PETSC_FALSE;

  save_snapshots = false;
}


IceModel::~IceModel() {
  if (createVecs_done == PETSC_TRUE) {
    //verbPrintf(1,grid.com, "calling destroyVecs()\n");
    destroyVecs();
  }
  if (createViewers_done == PETSC_TRUE) {
    //verbPrintf(1,grid.com, "calling destroyViewers()\n");
    destroyViewers();
  }
  if (createBasal_done == PETSC_TRUE) {
    //verbPrintf(1,grid.com, "deleting basal()\n");
    delete basal;
    //verbPrintf(1,grid.com, "deleting basalSIA()\n");
    delete basalSIA;
  }
  if (bootstrapLIC != PETSC_NULL) {
    delete bootstrapLIC;
  }
  //verbPrintf(1,grid.com, "Cleaning up the history string.\n");
  delete[] history;
}


//! Allocate all Vecs defined in IceModel.
/*! Initialization of an IceModel is confusing.  Here is a description of the intended order:
	\li 1. The constructor for IceModel.  Note IceModel has a member "grid", of class IceGrid. 
	   The IceGrid constructor sets 
	   defaults for (grid.)Mx,My,Mz,Mbz,Lx,Ly,Lz,Lbz,dx,dy,dz,year.
	\li [1.5] derivedClass::setFromOptions() to get options special to derived class
	\li 2. setFromOptions() to get all options *including* Mx,My,Mz,Mbz
	\li [2.5] initFromFile_netCDF() which reads Mx,My,Mz,Mbz from file and overwrites previous; if 
	   this represents a change the user is warned
	\li 3. createDA(), which uses only Mx,My,Mz,Mbz
	\li 4. createVecs() uses DA to create/allocate Vecs
	\li [4.5] derivedClass:createVecs() to create/allocate Vecs special to derived class
	\li 5. afterInitHook() which changes Lx,Ly,Lz if set by user

Note driver programs call only setFromOptions() and initFromOptions() (for IceModel 
or derived class).

Note IceModel::setFromOptions() should be called at the end of derivedClass:setFromOptions().

Note 2.5, 3, and 4 are called from initFromFile_netCDF() in IceModel.

Note 3 and 4 are called from initFromOptions() in some derived classes (e.g. IceCompModel) 
in cases where initFromFile_netCDF() is not called.

Note step 2.5 is skipped when bootstrapping (-bif and bootstrapFromFile_netCDF()) or in
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
  ierr =     v3.create(grid, "vvel", true); CHKERRQ(ierr);
  ierr =     v3.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
			  "m s-1", "land_ice_y_velocity"); CHKERRQ(ierr);

  ierr =     w3.create(grid, "wvel", false); CHKERRQ(ierr); // never diff'ed in hor dirs
  // PROPOSED standard name = land_ice_upward_velocity
  //   (compare "upward_air_velocity" and "upward_sea_water_velocity")
  ierr =     w3.set_attrs("diagnostic", "vertical velocity of ice",
			  "m s-1", NULL); CHKERRQ(ierr);
  ierr = Sigma3.create(grid, "Sigma", false); CHKERRQ(ierr); // never diff'ed in hor dirs

  // ice temperature
  ierr = T3.create(grid, "temp", true); CHKERRQ(ierr);
  ierr = T3.set_attrs("model_state","ice temperature",
		      "K", "land_ice_temperature"); CHKERRQ(ierr);

  // age of ice
  ierr = tau3.create(grid, "age", true); CHKERRQ(ierr);
  // PROPOSED standard_name = land_ice_age
  ierr = tau3.set_attrs("model_state", "age of ice",
			"s", NULL); CHKERRQ(ierr);

  // bedrock temperature
  ierr = Tb3.create(grid,"litho_temp", false); CHKERRQ(ierr);
  // PROPOSED standard_name = lithosphere_temperature
  ierr = Tb3.set_attrs("model_state", "lithosphere (bedrock) temperature",
		       "K", NULL); CHKERRQ(ierr);

  // ice upper surface elevation
  ierr = vh.create(grid, "usurf", true); CHKERRQ(ierr);
  ierr = vh.set_attrs("diagnostic", "ice upper surface elevation",
		      "m", "surface_altitude"); CHKERRQ(ierr);

  // land ice thickness
  ierr = vH.create(grid, "thk", true); CHKERRQ(ierr);
  ierr = vH.set_attrs("model_state", "land ice thickness",
		      "m", "land_ice_thickness"); CHKERRQ(ierr);

  // bedrock surface elevation
  ierr = vbed.create(grid, "topg", true); CHKERRQ(ierr);
  ierr = vbed.set_attrs("model_state", "bedrock surface elevation",
			"m", "bedrock_altitude"); CHKERRQ(ierr);

  // mean annual net ice equivalent accumulation (ablation) rate
  ierr = vAccum.create(grid, "acab", true); CHKERRQ(ierr);
  ierr = vAccum.set_attrs("climate_steady", "mean annual net ice equivalent accumulation (ablation) rate",
			  "m s-1", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);

  // annual mean air temperature at "ice surface", i.e.
  //   at level below all firn processes
  // PROPOSED standard_name = land_ice_temperature_below_firn
  ierr = vTs.create(grid, "artm", true); CHKERRQ(ierr);
  ierr = vTs.set_attrs("climate_steady", "temperature at ice surface but below firn",
		       "K", NULL); CHKERRQ(ierr);

  // grounded_dragging_floating integer mask
  ierr = vMask.create(grid, "mask",     true); CHKERRQ(ierr);
  ierr = vMask.set_attrs("model_state", "grounded_dragging_floating integer mask",
			 "", NULL); CHKERRQ(ierr);

  // upward geothermal flux at bedrock surface
  ierr = vGhf.create(grid, "bheatflx", true); CHKERRQ(ierr);
  ierr = vGhf.set_attrs("climate_steady", "upward geothermal flux at bedrock surface",
			"W m-2", NULL); CHKERRQ(ierr);

  // u bar and v bar
  ierr = vubar.create(grid, "ubar", true); CHKERRQ(ierr);
  ierr = vubar.set_attrs(NULL, "vertical mean of horizontal ice velocity in the X direction",
			  "m s-1", "land_ice_vertical_mean_x_velocity"); CHKERRQ(ierr);
  ierr = vubar.set_glaciological_units("m year-1", secpera);

  ierr = vvbar.create(grid, "vbar", true); CHKERRQ(ierr);
  ierr = vvbar.set_attrs(NULL, "vertical mean of horizontal ice velocity in the Y direction",
			  "m s-1", "land_ice_vertical_mean_y_velocity"); CHKERRQ(ierr);
  ierr = vvbar.set_glaciological_units("m year-1", secpera);

  // basal velocities on standard grid
  ierr = vub.create(grid, "ub", true); CHKERRQ(ierr);
  ierr = vub.set_attrs(NULL, "basal ice velocity in the X direction",
		       "m s-1", "land_ice_basal_x_velocity"); CHKERRQ(ierr);
  ierr = vvb.create(grid, "vb", true); CHKERRQ(ierr);
  ierr = vvb.set_attrs(NULL, "basal ice velocity in the Y direction",
		       "m s-1", "land_ice_basal_y_velocity"); CHKERRQ(ierr);

  // basal frictional heating on regular grid
  ierr = vRb.create(grid, "Rb", true); CHKERRQ(ierr);
  ierr = vRb.set_attrs(NULL, "basal frictional heating",
		       NULL, NULL); CHKERRQ(ierr);

  // effective thickness of subglacial melt water
  ierr = vHmelt.create(grid, "bwat", true); CHKERRQ(ierr);
  ierr = vHmelt.set_attrs("model_state", "effective thickness of subglacial melt water",
			  "m", NULL); CHKERRQ(ierr);

  // rate of change of ice thickness
  ierr = vdHdt.create(grid, "dHdt", true); CHKERRQ(ierr);
  ierr = vdHdt.set_attrs("diagnostic", "rate of change of ice thickness",
			 "m s-1", "tendency_of_land_ice_thickness"); CHKERRQ(ierr);
  ierr = vdHdt.set_glaciological_units("m year-1", secpera);

  // yield stress for basal till (plastic or pseudo-plastic model)
  ierr = vtauc.create(grid, "tauc", true); CHKERRQ(ierr);
  ierr = vtauc.set_attrs("diagnostic", "yield stress for basal till (plastic or pseudo-plastic model)",
			 "Pa", NULL); CHKERRQ(ierr);

  // bedrock uplift rate
  ierr = vuplift.create(grid, "dbdt", true); CHKERRQ(ierr);
  ierr = vuplift.set_attrs("model_state", "bedrock uplift rate",
			   "m s-1", "tendency_of_bedrock_altitude"); CHKERRQ(ierr);

  // basal melt rate
  ierr = vbasalMeltRate.create(grid, "basal_melt_rate", true); CHKERRQ(ierr);
  ierr = vbasalMeltRate.set_attrs(NULL, "basal melt rate",
				  "m s-1", "land_ice_basal_melt_rate"); CHKERRQ(ierr);

  // friction angle for till under grounded ice sheet
  ierr = vtillphi.create(grid, "tillphi", true);
  ierr = vtillphi.set_attrs("climate_steady", "friction angle for till under grounded ice sheet",
			    "degrees", NULL); CHKERRQ(ierr);

  // Longitude
  ierr = vLongitude.create(grid, "lon", true); CHKERRQ(ierr);
  ierr = vLongitude.set_attrs("mapping", "longitude", "degrees_east", "longitude"); CHKERRQ(ierr);

  // Latitude
  ierr = vLatitude.create(grid, "lat", true); CHKERRQ(ierr);
  ierr = vLatitude.set_attrs("mapping", "latitude", "degrees_north", "latitude"); CHKERRQ(ierr);

  // u bar and v bar on staggered grid
  ierr = vuvbar[0].create(grid, "vuvbar[0]", true); CHKERRQ(ierr);
  ierr = vuvbar[1].create(grid, "vuvbar[1]", true); CHKERRQ(ierr);

  // work vectors
  for (int j = 0; j < nWork2d; j++) {
    ierr = vWork2d[j].create(grid, "a_work_vector", true); CHKERRQ(ierr);
  }

  // initial guesses of SSA velocities
  ierr = vubarSSA.create(grid, "vubarSSA", true);
  ierr = vvbarSSA.create(grid, "vvbarSSA", true);


  ierr = Tnew3.createSameDA(T3,grid,"temp_new",false); CHKERRQ(ierr);
  ierr = taunew3.createSameDA(tau3,grid,"age_new",false); CHKERRQ(ierr);
  ierr = Sigmastag3[0].create(grid,"Sigma_stagx",true); CHKERRQ(ierr);
  ierr = Sigmastag3[1].create(grid,"Sigma_stagy",true); CHKERRQ(ierr);
  ierr = Istag3[0].create(grid,"I_stagx",true); CHKERRQ(ierr);
  ierr = Istag3[1].create(grid,"I_stagy",true); CHKERRQ(ierr);

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
  ierr = PDDCleanup(); CHKERRQ(ierr);

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
  ierr = vAccum.destroy(); CHKERRQ(ierr);
  ierr = vTs.destroy(); CHKERRQ(ierr);
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
    SETERRQ(1, "ERROR: ye < startYear.  PISM cannot run backward in time.\n");
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

#if (PISM_LOG_EVENTS)
PetscLogEventRegister(&siaEVENT,    "sia velocity",0);
PetscLogEventRegister(&ssaEVENT,    "ssa velocity",0);
PetscLogEventRegister(&velmiscEVENT,"misc vel calc",0);
PetscLogEventRegister(&beddefEVENT, "bed deform",0);
PetscLogEventRegister(&pddEVENT,    "pos deg day",0);
PetscLogEventRegister(&massbalEVENT,"mass bal calc",0);
PetscLogEventRegister(&tempEVENT,   "temp age calc",0);
#endif

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
    
    // read in forcing data if present; (typically from ice/seabed core;
    //   modifies vTs and seaLevel)
    ierr = updateForcing(); CHKERRQ(ierr);
    
#if (PISM_LOG_EVENTS)
PetscLogEventBegin(beddefEVENT,0,0,0,0);
#endif

    // compute bed deformation, which only depends on current thickness and bed elevation
    if (doBedDef == PETSC_TRUE) {
      ierr = bedDefStepIfNeeded(); CHKERRQ(ierr); // prints "b" or "$" as appropriate
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }
    
#if (PISM_LOG_EVENTS)
PetscLogEventEnd(beddefEVENT,0,0,0,0);
#endif

    // update basal till yield stress if appropriate; will modify and communicate mask
    if (doPlasticTill == PETSC_TRUE) {
      ierr = updateYieldStressFromHmelt();  CHKERRQ(ierr);
      ierr = verbPrintf(2,grid.com, "y"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }

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

#if (PISM_LOG_EVENTS)
PetscLogEventBegin(tempEVENT,0,0,0,0);
#endif
    
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

#if (PISM_LOG_EVENTS)
PetscLogEventEnd(tempEVENT,0,0,0,0);
PetscLogEventBegin(pddEVENT,0,0,0,0);
#endif

    // compute PDD; generates surface mass balance, with appropriate ablation area,
    //   using snow accumulation
    if (doPDD == PETSC_TRUE) {
      ierr = updateSurfaceBalanceFromPDD();  CHKERRQ(ierr);
      ierr = verbPrintf(2,grid.com, "d"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }
    
#if (PISM_LOG_EVENTS)
PetscLogEventEnd(pddEVENT,0,0,0,0);
PetscLogEventBegin(massbalEVENT,0,0,0,0);
#endif

    if (doMassConserve == PETSC_TRUE) {
      ierr = massContExplicitStep(); CHKERRQ(ierr); // update H
      ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr); // update h and mask
      if ((doSkip == PETSC_TRUE) && (skipCountDown > 0))
        skipCountDown--;
      ierr = verbPrintf(2,grid.com, "h"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }

#if (PISM_LOG_EVENTS)
PetscLogEventEnd(massbalEVENT,0,0,0,0);
#endif
    
    ierr = additionalAtEndTimestep(); CHKERRQ(ierr);

    // end the flag line and report a summary
    ierr = verbPrintf(2,grid.com, " %d%c  +%6.5f\n", skipCountDown, adaptReasonFlag,
                      dt / secpera); CHKERRQ(ierr);
    ierr = summary(tempAgeStep,reportHomolTemps); CHKERRQ(ierr);

    ierr = updateViewers(); CHKERRQ(ierr);

    if (endOfTimeStepHook() != 0) break;
  }
  
  ierr = forcingCleanup(); CHKERRQ(ierr);  // puts back bed and Ts (removes offsets)

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


// note no range checking in these two:
int IceModel::intMask(PetscScalar maskvalue) {
  return static_cast<int>(floor(maskvalue + 0.5));
}


int IceModel::modMask(PetscScalar maskvalue) {
  int intmask = static_cast<int>(floor(maskvalue + 0.5));
  if (intmask > MASK_FLOATING) {
    return intmask - 4;
  } else {
    return intmask;
  }
}


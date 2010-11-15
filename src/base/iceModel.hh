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

#ifndef __iceModel_hh
#define __iceModel_hh

#include <signal.h>
#include <gsl/gsl_rng.h>
#include <petscsnes.h>
#include <petsctime.h>		// PetscGetTime()

#include "flowlaw_factory.hh"   // IceFlowLawFactory and friends
#include "materials.hh"         // IceBasalResistancePlasticLaw, SSAStrengthExtension
#include "pism_const.hh"
#include "enthalpyConverter.hh"
#include "grid.hh"
#include "iceModelVec.hh"
#include "NCVariable.hh"
#include "PISMBedSmoother.hh"
#include "PISMVars.hh"
#include "Timeseries.hh"
#include "PISMProf.hh"

#include "../earth/PISMBedDef.hh"
#include "../coupler/PISMOcean.hh"
#include "../coupler/PISMSurface.hh"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

//! The base class for PISM.  Contains all essential variables, parameters, and flags for modelling an ice sheet.
class IceModel {
public:
  // see iceModel.cc for implementation of constructor and destructor:
  IceModel(IceGrid &g, NCConfigVariable &config, NCConfigVariable &overrides);
  virtual ~IceModel(); // must be virtual merely because some members are virtual

  // see iMinit.cc
  virtual PetscErrorCode grid_setup();
  virtual PetscErrorCode init_physics();
  virtual PetscErrorCode init_couplers();
  virtual PetscErrorCode set_grid_from_options();
  virtual PetscErrorCode set_grid_defaults();
  virtual PetscErrorCode model_state_setup();
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode report_grid_parameters(); // should be moved into a method of IceGrid
  virtual PetscErrorCode allocate_internal_objects();
  virtual PetscErrorCode misc_setup();


  // see iceModel.cc
  PetscErrorCode init();
  virtual PetscErrorCode run();
  virtual PetscErrorCode step(bool do_mass_continuity, 
                              bool do_energy,
			      bool do_age,
			      bool do_skip,
			      bool use_ssa_when_grounded);
  virtual PetscErrorCode setExecName(const char *my_executable_short_name);
  virtual IceFlowLawFactory &getIceFlowLawFactory() { return iceFactory; }
  virtual IceFlowLaw *getIceFlowLaw() {return ice;}

  // see iMbootstrap.cc 
  virtual PetscErrorCode bootstrapFromFile(const char *fname);
  virtual PetscErrorCode readShelfStreamBCFromFile(const char *fname); // should be moved into the SSA 
  // stress balance module

  // see iMoptions.cc
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode set_output_size(string option, string description,
					 string default_value, set<string> &result);

  // see iMutil.cc
  virtual void attach_surface_model(PISMSurfaceModel *surf);
  virtual void attach_ocean_model(PISMOceanModel *ocean);
  virtual PetscErrorCode additionalAtStartTimestep();
  virtual PetscErrorCode additionalAtEndTimestep();
  virtual PetscErrorCode compute_cell_areas(); // is an initialization step; should go there

  // see iMIO.cc
  virtual PetscErrorCode initFromFile(const char *);
  virtual PetscErrorCode writeFiles(const char* default_filename);
  virtual PetscErrorCode write_model_state(const char *filename);
  virtual PetscErrorCode write_variables(const char* filename, set<string> vars,
					 nc_type nctype);
  virtual PetscErrorCode write_extra_fields(const char *filename); // should be deprecated

protected:

  IceGrid               &grid;

  NCConfigVariable      mapping, //!< grid projection (mapping) parameters
    &config,			 //!< configuration flags and parameters
    &overrides;			 //!< flags and parameters overriding config, see -config_override
  NCGlobalAttributes    global_attributes;

  IceFlowLawFactory     iceFactory;
  IceFlowLaw            *ice;
  IceBasalResistancePlasticLaw *basal;
  SSAStrengthExtension  ssaStrengthExtend;

  EnthalpyConverter *EC;
 
  PISMBedSmoother   *sia_bed_smoother; // should probably go inside the SIA stress balance module

  PISMSurfaceModel *surface;
  PISMOceanModel   *ocean;
  PISMBedDef       *beddef;

  //! \brief A dictionary with pointers to IceModelVecs below, for passing them
  //! from the IceModel core to other components (such as surface and ocean models)
  PISMVars variables;

  // state variables
  IceModelVec2S
        vh,		//!< ice surface elevation; ghosted
        vH,		//!< ice thickness; ghosted
        vdHdt,		//!< \f$ \frac{dH}{dt} \f$; ghosted to simplify the code computing it
        vtauc,		//!< yield stress for basal till (plastic or pseudo-plastic model); no ghosts
        vHmelt,		//!< thickness of the basal meltwater; ghosted (because of the diffusion)
        vbmr,	    //!< rate of production of basal meltwater (ice-equivalent); no ghosts
        vLongitude,	//!< Longitude; ghosted to compute cell areas
        vLatitude,	//!< Latitude; ghosted to compute cell areas
        vbed,		//!< bed topography; ghosted
        vuplift,	//!< bed uplift rate; ghosted to simplify the code computing it
        vGhf,		//!< geothermal flux; no ghosts
        vRb,		//!< basal frictional heating on regular grid; no ghosts
        vtillphi,	//!< friction angle for till under grounded ice sheet; no ghosts
    acab,		//!< accumulation/ablation rate; no ghosts
    artm,		//!< ice temperature at the ice surface but below firn; no ghosts
    shelfbtemp,		//!< ice temperature at the shelf base; no ghosts
    shelfbmassflux,	//!< ice mass flux into the ocean at the shelf base; no ghosts
    cell_area;		//!< cell areas (computed using the WGS84 datum)

  IceModelVec2Stag uvbar, //!< ubar and vbar on staggered grid; ubar at i+1/2, vbar at j+1/2
                   vWork2dStag; //!< used for diffusivity and averaged ice hardness on the staggered grid
 
  IceModelVec2V vel_basal,	//!< basal sliding velocities on standard grid; ghosted
    vel_bar; //!< vertically-averaged horizontal velocity on standard grid; ghosted

  IceModelVec2Mask vMask; //!< mask for flow type with values SHEET, DRAGGING, FLOATING

  IceModelVec3
        u3, v3, w3,	//!< velocity of ice; m s-1 (ghosted)
        Sigma3, 	//!< strain-heating term in conservation of energy model; J s-1 m-3 (no ghosts)
        T3,		//!< absolute temperature of ice; K (ghosted)
        Enth3,          //!< enthalpy; J / kg (ghosted)
        tau3;		//!< age of ice; s (ghosted because it is evaluated on the staggered-grid)

  IceModelVec3Bedrock
        Tb3;		//!< temperature of lithosphere (bedrock) under ice or ocean; K (no ghosts)

  // parameters
  PetscReal   dt,     //!< mass continuity time step, s
              dtTempAge,  //!< enthalpy/temperature and age time-steps in seconds
              maxdt_temporary, dt_force,
              CFLviolcount,    //!< really is just a count, but PetscGlobalSum requires this type
              dt_from_diffus, dt_from_cfl, CFLmaxdt, CFLmaxdt2D,
              gDmax,		// global max of the diffusivity
              gmaxu, gmaxv, gmaxw,  // global maximums on 3D grid of abs value of vel components
              gdHdtav,  //!< average value in map-plane (2D) of dH/dt, where there is ice; m s-1
    total_sub_shelf_ice_flux,
    total_basal_ice_flux,
    total_surface_ice_flux,
    nonneg_rule_flux,
    ocean_kill_flux,
    float_kill_flux,
    dvoldt;  //!< d(total ice volume)/dt; m3 s-1
  PetscInt    skipCountDown;

  // physical parameters used frequently enough to make looking up via
  // config.get() a hassle:
  PetscScalar standard_gravity;
  // Initialized in the IceModel constructor from the configuration file;
  // SHOULD NOT be hard-wired.

  // flags
  bool leaveNuHAloneSSA;        // should be removed
  PetscTruth  updateHmelt,
              holdTillYieldStress,
              shelvesDragToo,
              doAdaptTimeStep, 
              realAgeForGrainSize, // is always equal to PETSC_FALSE; should we let it go?
              ssaSystemToASCIIMatlab,
              allowAboveMelting,
              computeSIAVelocities;
  char        adaptReasonFlag;

  char        ssaMatlabFilePrefix[PETSC_MAX_PATH_LEN];
  string      stdout_flags, stdout_ssa;

  string executable_short_name;
  
protected:
  // see iceModel.cc
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode deallocate_internal_objects();
  virtual void setConstantNuHForSSA(PetscScalar); // move into the SSA module

  // see iMadaptive.cc
  virtual PetscErrorCode computeMax3DVelocities();
  virtual PetscErrorCode computeMax2DSlidingSpeed();
  virtual PetscErrorCode adaptTimeStepDiffusivity();
  virtual PetscErrorCode determineTimeStep(const bool doTemperatureCFL);
  virtual PetscErrorCode countCFLViolations(PetscScalar* CFLviol);

  // see iMbasal.cc: all relate to grounded SSA
  virtual PetscErrorCode initBasalTillModel();
  virtual PetscErrorCode computePhiFromBedElevation();
  virtual PetscScalar    getBasalWaterPressure(
                           PetscScalar thk, PetscScalar bwat, PetscScalar bmr,
                           PetscScalar frac, PetscScalar hmelt_max) const;
  virtual PetscErrorCode updateYieldStressUsingBasalWater();
  virtual PetscScalar    basalDragx(PetscScalar **tauc, PISMVector2 **uv,
				    PetscInt i, PetscInt j) const;
  virtual PetscScalar    basalDragy(PetscScalar **tauc, PISMVector2 **uv,
				    PetscInt i, PetscInt j) const;
  virtual PetscErrorCode diffuseHmelt();

  // see iMbeddef.cc
  PetscScalar last_bed_def_update;
  virtual PetscErrorCode bed_def_setup();
  virtual PetscErrorCode bed_def_step(bool &bed_changed);

  // see iMbootstrap.cc 
  virtual PetscErrorCode putTempAtDepth();
  virtual PetscErrorCode bootstrapSetBedrockColumnTemp(PetscInt i, PetscInt j,
						       PetscScalar Ttopbedrock,
						       PetscScalar geothermflux,
						       PetscScalar bed_thermal_k);
  virtual PetscErrorCode setMaskSurfaceElevation_bootstrap();

  // see iMdefaults.cc
  PetscErrorCode setDefaults();

  // see iMenthalpy.cc
  virtual PetscErrorCode setEnth3FromT3_ColdIce();
  virtual PetscErrorCode setEnth3FromT3AndLiqfrac3(IceModelVec3 &Liqfrac3);
  virtual PetscErrorCode setLiquidFracFromEnthalpy(IceModelVec3 &useForLiquidFrac);
  virtual PetscErrorCode setCTSFromEnthalpy(IceModelVec3 &useForCTS);
  virtual PetscErrorCode getEnthalpyCTSColumn(PetscScalar p_air, //!< atmospheric pressure
					      PetscScalar thk,	 //!< ice thickness
					      PetscInt ks,	 //!< index of the level just below the surface
					      const PetscScalar *Enth, //!< Enthalpy
					      const PetscScalar *w,    //!< vert. velocity
					      PetscScalar *lambda,     //!< lambda
					      PetscScalar **Enth_s//!< enthalpy for the pr.-melting temp.
                                              );   
  virtual PetscErrorCode enthalpyAndDrainageStep(
                PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                PetscScalar* bulgeCount);

  // see iMgeometry.cc
  virtual PetscErrorCode computeDrivingStress(IceModelVec2S &vtaudx, IceModelVec2S &vtaudy);
  virtual PetscErrorCode updateSurfaceElevationAndMask();
  virtual PetscErrorCode update_mask();
  virtual PetscErrorCode update_surface_elevation();
  virtual PetscErrorCode massContExplicitStep();

  // see iMgrainsize.cc
  virtual PetscScalar    grainSizeVostok(PetscScalar age) const;

  // see iMIO.cc
  virtual PetscErrorCode set_time_from_options();
  virtual PetscErrorCode dumpToFile(const char *filename);
  virtual PetscErrorCode regrid();

  // see iMmatlab.cc
  virtual PetscErrorCode writeSSAsystemMatlab(IceModelVec2S vNuH[2]);

  // see iMreport.cc
  virtual PetscErrorCode computeFlowUbarStats(
                       PetscScalar *gUbarmax, PetscScalar *gUbarSIAav,
                       PetscScalar *gUbarstreamav, PetscScalar *gUbarshelfav,
                       PetscScalar *gicegridfrac, PetscScalar *gSIAgridfrac,
                       PetscScalar *gstreamgridfrac, PetscScalar *gshelfgridfrac);
  virtual PetscErrorCode volumeArea(
                       PetscScalar& gvolume,PetscScalar& garea,
                       PetscScalar& gvolSIA, PetscScalar& gvolstream, 
                       PetscScalar& gvolshelf);
  virtual PetscErrorCode energyStats(
                       PetscScalar iarea,
                       PetscScalar &gmeltfrac, PetscScalar &gtemp0);
  virtual PetscErrorCode ageStats(PetscScalar ivol, PetscScalar &gorigfrac);
  virtual PetscErrorCode summary(bool tempAndAge);
  virtual PetscErrorCode summaryPrintLine(
              PetscTruth printPrototype, bool tempAndAge,
              PetscScalar year, PetscScalar delta_t, 
              PetscScalar volume, PetscScalar area,
              PetscScalar meltfrac, PetscScalar H0, PetscScalar T0);

  // see iMreport.cc;  methods for computing diagnostic quantities:
  // spatially-varying:
  virtual PetscErrorCode compute_by_name(string name, IceModelVec* &result);
  virtual PetscErrorCode compute_bwp(IceModelVec2S &result);
  virtual PetscErrorCode compute_cbar(IceModelVec2S &result);
  virtual PetscErrorCode compute_cbase(IceModelVec2S &result, IceModelVec2S &tmp);
  virtual PetscErrorCode compute_cflx(IceModelVec2S &result, IceModelVec2S &cbar);
  virtual PetscErrorCode compute_csurf(IceModelVec2S &result, IceModelVec2S &tmp);
  virtual PetscErrorCode compute_cts(IceModelVec3 &useForCTS);
  virtual PetscErrorCode compute_dhdt(IceModelVec2S &result);
  virtual PetscErrorCode compute_diffusivity(IceModelVec2S &result);
  virtual PetscErrorCode compute_enthalpybase(IceModelVec2S &result);
  virtual PetscErrorCode compute_enthalpysurf(IceModelVec2S &result);
  virtual PetscErrorCode compute_hardav(IceModelVec2S &result);
  virtual PetscErrorCode compute_liqfrac(IceModelVec3 &useForLiqfrac);
  virtual PetscErrorCode compute_schoofs_theta(IceModelVec2S &result);
  virtual PetscErrorCode compute_taud(IceModelVec2S &result, IceModelVec2S &tmp);
  virtual PetscErrorCode compute_temp(IceModelVec3 &result);
  virtual PetscErrorCode compute_temp_pa(IceModelVec3 &result);
  virtual PetscErrorCode compute_tempbase(IceModelVec2S &result);
  virtual PetscErrorCode compute_tempicethk(IceModelVec2S &result);
  virtual PetscErrorCode compute_tempicethk_basal(IceModelVec2S &result);
  virtual PetscErrorCode compute_temppabase(IceModelVec3 &hasPATemp,
                                            IceModelVec2S &result);
  virtual PetscErrorCode compute_tempsurf(IceModelVec2S &result);
  virtual PetscErrorCode compute_thksmooth(IceModelVec2S &result);
  virtual PetscErrorCode compute_topgsmooth(IceModelVec2S &result);
  virtual PetscErrorCode compute_uvelbase(IceModelVec2S &result);
  virtual PetscErrorCode compute_uvelsurf(IceModelVec2S &result);
  virtual PetscErrorCode compute_vvelbase(IceModelVec2S &result);
  virtual PetscErrorCode compute_vvelsurf(IceModelVec2S &result);
  virtual PetscErrorCode compute_wvel(IceModelVec3 &result);
  virtual PetscErrorCode compute_wvelbase(IceModelVec3 &wvel, IceModelVec2S &result);
  virtual PetscErrorCode compute_wvelsurf(IceModelVec3 &wvel, IceModelVec2S &result);

  // profiling, etc:
  virtual PetscErrorCode compute_rank(IceModelVec2S &result);
  virtual PetscErrorCode compute_proc_ice_area(IceModelVec2S &result);
  virtual PetscErrorCode compute_bueler_brown_f(IceModelVec2S &result);
  // scalar:
  virtual PetscErrorCode ice_mass_bookkeeping();
  virtual PetscErrorCode compute_ice_volume(PetscScalar &result);
  virtual PetscErrorCode compute_ice_volume_temperate(PetscScalar &result);
  virtual PetscErrorCode compute_ice_volume_cold(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area_temperate(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area_cold(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area_grounded(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area_floating(PetscScalar &result);
  virtual PetscErrorCode compute_ice_enthalpy(PetscScalar &result);
  virtual PetscErrorCode compute_by_name(string name, PetscScalar &result);

  // see iMsia.cc
  virtual PetscErrorCode surfaceGradientSIA();

  virtual PetscErrorCode surface_gradient_eta();
  virtual PetscErrorCode surface_gradient_haseloff();
  virtual PetscErrorCode surface_gradient_mahaffy();

  virtual PetscErrorCode velocitySIAStaggered();
  virtual PetscErrorCode basalSlidingHeatingSIA();
  virtual PetscErrorCode velocities2DSIAToRegular();
  virtual PetscErrorCode SigmaSIAToRegular();
  virtual PetscErrorCode horizontalVelocitySIARegular();
  virtual PetscScalar    basalVelocitySIA( // not recommended, generally
                             PetscScalar x, PetscScalar y, PetscScalar H, PetscScalar T,
                             PetscScalar alpha, PetscScalar mu, PetscScalar min_T) const;

  // see iMssa.cc
  virtual PetscScalar    bueler_brown_f(PetscScalar v_squared);
  virtual PetscErrorCode allocateSSAobjects();
  virtual PetscErrorCode destroySSAobjects();
  virtual PetscErrorCode initSSA();
  virtual PetscErrorCode velocitySSA(PetscInt *numiter);
  virtual PetscErrorCode velocitySSA(IceModelVec2S vNuH[2], PetscInt *numiter);
  virtual PetscErrorCode computeEffectiveViscosity(IceModelVec2S vNuH[2], PetscReal epsilon);
  virtual PetscErrorCode testConvergenceOfNu(IceModelVec2S vNuH[2], IceModelVec2S vNuHOld[2],
                                             PetscReal *norm, PetscReal *normChange);
  virtual PetscErrorCode assembleSSAMatrix(bool includeBasalShear, IceModelVec2S vNuH[2], Mat A);
  virtual PetscErrorCode assembleSSARhs(Vec rhs);
  virtual PetscErrorCode trivialMoveSSAXtoIMV2V();
  virtual PetscErrorCode broadcastSSAVelocity(bool updateVelocityAtDepth);
  virtual PetscErrorCode correctSigma();
  virtual PetscErrorCode correctBasalFrictionalHeating();
  virtual PetscErrorCode compute_hardav_staggered(IceModelVec2Stag &result);
  
  // see iMtemp.cc
  virtual PetscErrorCode energyStep();
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);
  virtual PetscErrorCode ageStep();
  virtual bool checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                      PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE);
  virtual PetscErrorCode excessToFromBasalMeltLayer(
                      const PetscScalar rho_c, const PetscScalar z, const PetscScalar dz,
                      PetscScalar *Texcess, PetscScalar *Hmelt);

  // see iMutil.cc
  virtual int            endOfTimeStepHook();
  virtual PetscErrorCode stampHistoryCommand();
  virtual PetscErrorCode stampHistoryEnd();
  virtual PetscErrorCode stampHistory(string);
  virtual PetscErrorCode check_maximum_thickness();
  virtual PetscErrorCode check_maximum_thickness_hook(const int old_Mz);
  virtual bool           issounding(const PetscInt i, const PetscInt j);

  // see iMvelocity.cc
  virtual PetscErrorCode velocity(bool updateSIAVelocityAtDepth);    
  virtual PetscErrorCode vertVelocityFromIncompressibility();

protected:
  // working space (a convenience)
  static const PetscInt nWork2d=6;
  IceModelVec2S vWork2d[nWork2d];
  // 3D working space (with specific purposes)
  IceModelVec3 vWork3d;
  IceModelVec3 Sigmastag3[2], Istag3[2], DeltaStag3[2];

  // for saving SSA velocities for restart
  bool have_ssa_velocities;	//!< use ubar_ssa and vbar_ssa from a previous
				//! run if true, otherwise set them to zero in
				//! IceModel::initSSA()
  IceModelVec2V vel_ssa, vel_ssa_old;

  // SSA solve vars; note pieces of the SSA Velocity routine are defined in iMssa.cc
  KSP SSAKSP;
  Mat SSAStiffnessMatrix;
  Vec SSAX, SSARHS;  // Global vectors for solution of the linear system
  DA  SSADA;

  // Set of variables to put in the output file:
  set<string> output_vars;

  // This is related to the snapshot saving feature
  string snapshots_filename;
  bool save_snapshots, snapshots_file_is_ready, split_snapshots;
  vector<double> snapshot_times;
  set<string> snapshot_vars;
  unsigned int current_snapshot;
  PetscErrorCode init_snapshots();
  PetscErrorCode write_snapshot();

  // scalar time-series
  bool save_ts;			//! true if the user requested time-series output
  string ts_filename;		//! file to write time-series to
  vector<double> ts_times;	//! times requested
  unsigned int current_ts;	//! index of the current time
  set<string> ts_vars;		//! variables requested
  vector<DiagnosticTimeseries*> timeseries;
  PetscErrorCode init_timeseries();
  PetscErrorCode create_timeseries();
  PetscErrorCode flush_timeseries();
  PetscErrorCode write_timeseries();
  PetscErrorCode ts_max_timestep(double t_years, double& dt_years);

  // spatially-varying time-series
  bool save_extra, extra_file_is_ready, split_extra;
  string extra_filename;
  vector<double> extra_times;
  unsigned int current_extra;
  set<string> extra_vars;
  PetscErrorCode init_extras();
  PetscErrorCode write_extras();
  PetscErrorCode extras_max_timestep(double t_years, double& dt_years);

  // automatic backups
  double backup_interval;
  string backup_filename;
  PetscReal last_backup_time;
  PetscErrorCode init_backups();
  PetscErrorCode write_backup();

  // diagnostic viewers; see iMviewers.cc
  virtual PetscErrorCode init_viewers();
  virtual PetscErrorCode update_viewers();
  virtual PetscErrorCode update_nu_viewers(IceModelVec2S vNu[2]);
  set<string> map_viewers, slice_viewers, sounding_viewers;
  PetscInt     id, jd;	     // sounding indices
  bool view_diffusivity, view_log_nuH, view_nuH;

private:
  // for event logging (profiling); see run() and velocity()
  int siaEVENT, ssaEVENT, velmiscEVENT, beddefEVENT, massbalEVENT, tempEVENT;
  PetscLogDouble start_time;

  PISMProf *prof;
  int event_step,		//!< total time spent doing time-stepping
    event_velocity,		//!< total velocity computation
    event_vel_inc,		//!< vert. velocity from incompressibility
    event_sia,			//!< SIA computation
    event_sia_2d,               //!< 2D SIA computation
    event_sia_3d,               //!< filling u3 and v3 with SIA velocities
    event_ssa,			//!< SSA computation
    event_ssa_3d,               //!< combining SIA and SSA and filling u3 and v3
    event_energy,		//!< energy balance computation
    event_vel_com,		//!< velocity communication
    event_thk_com,		//!< thickness communication
    event_mass,			//!< mass balance computation
    event_age,			//!< age computation
    event_beddef,		//!< bed deformation step
    event_output;		//!< time spent writing the output file
};

#endif /* __iceModel_hh */


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

#include "materials.hh"
#include "pism_const.hh"
#include "grid.hh"
#include "iceModelVec.hh"
#include "NCVariable.hh"
#include "PISMVars.hh"
#include "Timeseries.hh"

#include "../earth/deformation.hh"

#include "../coupler/PISMOcean.hh"
#include "../coupler/PISMSurface.hh"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

struct SSASNESNode {
  PetscScalar u, v;
};
struct SSASNESCtx {
  DA               ssada;
  IceGrid          *grid;
  IceBasalResistancePlasticLaw *basal;
  IceModelVec2     ctxH,
                   ctxtauc,
                   ctxtaudx,
                   ctxtaudy,
                   ctxNu[2];
  IceModelVec2Mask ctxMask;
  Vec              ctxBV;
  PetscScalar      schoofReg,
                   constantHardness;
  bool             useConstantHardness,
                   useConstantNu,
                   useStoredNu,
                   useIceBasalResistancePlasticLaw;
};


// two structs needed in iMinverse.cc, iMinverseMat.cc
struct InverseModelCtx {
  // all the fields involved in the inverse model which determines till friction angle phi
  IceModelVec2 *usIn,     // x component of observed surf vel
               *vsIn,     // y component of observed surf vel
               *invMask,  // either 0, 1, or 2; 0 if no valid velocity, 1 if valid, 2 if valid
                          // and stencil width
               *usSIA,
               *vsSIA, 
               *fofv,
               *taubxComputed,
               *taubyComputed,
               *effPressureN,
               *oldtillphi;
};
struct RegPoissonCtx {
  // describes Poisson-like problem solved when inverting surface velocities
  //   using a regularization
  DA           da;              // must be first in struct
  IceGrid      *grid;
  PetscReal    epsilon;
  IceModelVec2 *f,              // = f(x,y) in PDE
               *g;              // = g(x,y) in PDE
};


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
  virtual PetscErrorCode report_grid_parameters();
  virtual PetscErrorCode allocate_internal_objects();
  virtual PetscErrorCode misc_setup();


  // see iceModel.cc
  PetscErrorCode init();
  virtual PetscErrorCode run();
  virtual PetscErrorCode step(bool do_mass_conserve, bool do_temp,
			      bool do_age, bool do_skip,
			      bool do_bed_deformation,
			      bool do_plastic_till);
  virtual PetscErrorCode diagnosticRun();
  virtual PetscErrorCode setExecName(const char *my_executable_short_name);
  virtual IceFlowLawFactory &getIceFlowLawFactory() { return iceFactory; }
  virtual IceFlowLaw *getIceFlowLaw() {return ice;}

  // see iMbootstrap.cc 
  virtual PetscErrorCode bootstrapFromFile(const char *fname);
  virtual PetscErrorCode readShelfStreamBCFromFile(const char *fname);

  // see iMoptions.cc
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode set_output_size(string option, string description,
					 string default_value, set<string> &result);

  // see iMutil.cc
  virtual void attach_surface_model(PISMSurfaceModel *surf);
  virtual void attach_ocean_model(PISMOceanModel *ocean);
  virtual PetscErrorCode additionalAtStartTimestep();
  virtual PetscErrorCode additionalAtEndTimestep();

  // see iMIO.cc
  virtual PetscErrorCode initFromFile(const char *);
  virtual PetscErrorCode writeFiles(const char* default_filename);
  virtual PetscErrorCode write_model_state(const char *filename);
  virtual PetscErrorCode write_variables(const char* filename, set<string> vars,
					 nc_type nctype);
  virtual PetscErrorCode write_extra_fields(const char *filename);

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

  PISMSurfaceModel *surface;
  PISMOceanModel   *ocean;

  //! \brief A dictionary with pointers to IceModelVecs below, for passing them
  //! from the IceModel core to other components (such as surface and ocean models)
  PISMVars variables;

  InverseModelCtx       inv;
  
  // state variables
  IceModelVec2
        vh,		//!< ice surface elevation
        vH,		//!< ice thickness
        vdHdt,		//!< \f$ \frac{dH}{dt} \f$
        vtauc,		//!< yield stress for basal till (plastic or pseudo-plastic model)
        vHmelt,		//!< thickness of the basal meltwater
        vbasalMeltRate,	//!< rate of production of basal meltwater (ice-equivalent)
        vLongitude,	//!< Longitude
        vLatitude,	//!< Latitude 
        vbed,		//!< bed topography
        vuplift,	//!< bed uplift rate
        vGhf,		//!< geothermal flux
        vRb,		//!< basal frictional heating on regular grid
        vtillphi,	//!< friction angle for till under grounded ice sheet
        vuvbar[2],	//!< ubar and vbar on staggered grid; ubar at i+1/2, vbar at j+1/2
        vub, vvb,	//!< basal velocities on standard grid
        vubar, vvbar,	//!< vertically-averaged horizontal velocity on standard grid
    acab,
    artm,
    shelfbtemp,
    shelfbmassflux;

  IceModelVec2Mask vMask; //!< mask for flow type with values SHEET, DRAGGING, FLOATING

  IceModelVec3
        u3, v3, w3,	//!< velocity of ice; m s-1
        Sigma3, 	//!< strain-heating term in conservation of energy model; J s-1 m-3
        T3,		//!< absolute temperature of ice; K
        Enth3,          //!< FIXME: not functional in IceModel; check initialization when IceEnthalpyModel moved
        tau3;		//!< age of ice; s

  IceModelVec3Bedrock
        Tb3;		//!< temperature of lithosphere (bedrock) under ice or ocean; K

  // parameters
  PetscReal   dt, dtTempAge,  // current mass cont. and temp/age; time steps in seconds
              maxdt_temporary, dt_force,
              CFLviolcount,    //!< really is just a count, but PetscGlobalSum requires this type
              dt_from_diffus, dt_from_cfl, CFLmaxdt, CFLmaxdt2D, gDmax,
              gmaxu, gmaxv, gmaxw,  // global maximums on 3D grid of abs value of vel components
              gdHdtav,  //!< average value in map-plane (2D) of dH/dt, where there is ice; m s-1
    total_sub_shelf_ice_flux,
    total_basal_ice_flux,
    total_surface_ice_flux,
    dvoldt;  //!< d(total ice volume)/dt; m3 s-1
  PetscInt    skipCountDown;

  // physical parameters used frequently enough to make looking up via
  // config.get() a hassle:
  PetscScalar standard_gravity;
  // Initialized in the IceModel constructor from the configuration file;
  // SHOULD NOT be hard-wired.

  // flags
  bool leaveNuHAloneSSA;
  PetscTruth  updateHmelt,
              holdTillYieldStress,
              useConstantTillPhi,
              shelvesDragToo,
              doAdaptTimeStep, 
              realAgeForGrainSize,
              ssaSystemToASCIIMatlab,
              reportPATemps,
              allowAboveMelting,
              computeSIAVelocities,
              transformForSurfaceGradient,
              doColdIceMethods; // if true, use cold ice internals but read and write
                                //        additional enthalpy fields to and from file
                                // FIXME: check initialization when IceEnthalpyModel moved
  char        adaptReasonFlag;

  char        ssaMatlabFilePrefix[PETSC_MAX_PATH_LEN];
  string      stdout_flags, stdout_ssa;

  string executable_short_name;
  
protected:
  // see iceModel.cc
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode deallocate_internal_objects();
  virtual void setConstantNuHForSSA(PetscScalar);

  // see iMadaptive.cc
  virtual PetscErrorCode computeMaxDiffusivity(bool update_diffusivity_viewer);
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
  virtual PetscScalar    basalDragx(
                           PetscScalar **tauc, PetscScalar **u, PetscScalar **v,
                           PetscInt i, PetscInt j) const;
  virtual PetscScalar    basalDragy(
                           PetscScalar **tauc, PetscScalar **u, PetscScalar **v,
                           PetscInt i, PetscInt j) const;
  virtual PetscErrorCode diffuseHmelt();

  // see iMbeddef.cc
  Vec            g2natural;
  VecScatter     top0ctx; // possibly useful general tool for putting Vecs on processor zero
  PetscErrorCode createScatterToProcZero(Vec& samplep0);
  PetscErrorCode destroyScatterToProcZero();
  BedDeformLC    bdLC;
  PetscScalar    lastBedDefUpdateYear;
  IceModelVec2   vbedlast;
  IceModelVec2   vHlast;	//!< used for simple pointwise isostasy and to compute uplift
  Vec            Hp0, bedp0,                       // vecs on proc zero for
                 Hstartp0, bedstartp0, upliftp0;   // passing to bdLC
  virtual PetscErrorCode bedDefSetup();
  virtual PetscErrorCode bedDefCleanup();
  virtual PetscErrorCode bedDefStepIfNeeded();
  virtual PetscErrorCode bed_def_step_iso();

  // see iMbootstrap.cc 
  virtual PetscErrorCode putTempAtDepth();
  virtual PetscErrorCode bootstrapSetBedrockColumnTemp(PetscInt i, PetscInt j,
						       PetscScalar Ttopbedrock,
						       PetscScalar geothermflux,
						       PetscScalar bed_thermal_k);
  virtual PetscErrorCode setMaskSurfaceElevation_bootstrap();

  // see iMdefaults.cc
  PetscErrorCode setDefaults();

  // see iMgeometry.cc
  virtual PetscErrorCode computeDrivingStress(IceModelVec2 &vtaudx, IceModelVec2 &vtaudy);
  virtual PetscErrorCode updateSurfaceElevationAndMask();
  virtual PetscErrorCode massContExplicitStep();

  // see iMgrainsize.cc
  virtual PetscScalar    grainSizeVostok(PetscScalar age) const;

  // see iMinverse.cc
  virtual PetscErrorCode invertSurfaceVelocities(const char *filename);
  virtual PetscErrorCode createInvFields();
  virtual PetscErrorCode destroyInvFields();
  virtual PetscErrorCode writeInvFields(const char *filename);
  virtual PetscErrorCode readObservedSurfVels(const char *filename);
  virtual PetscErrorCode smoothObservedSurfVels(const PetscInt passes);
  virtual PetscErrorCode computeSIASurfaceVelocity();
  virtual PetscErrorCode getGforInverse(
                const PetscScalar x, const PetscScalar UsuSIAdiffsqr, 
                const PetscScalar UsuSIAdiffdotuSIA, const PetscScalar uSIAsqr,
                PetscScalar &G, PetscScalar &Gprime);
  virtual PetscErrorCode computeFofVforInverse();
  virtual PetscErrorCode removeSIApart();
  virtual PetscErrorCode getEffectivePressureForInverse();
  virtual PetscErrorCode getVfromUforInverse(
                const PetscScalar U_x, const PetscScalar U_y,
                PetscScalar &V_x, PetscScalar &V_y, PetscScalar &magVsqr);
  virtual PetscErrorCode computeTFAFromBasalShearNoReg(
                const PetscScalar phi_low, const PetscScalar phi_high);

  // see iMinverseMat.cc
  virtual PetscErrorCode computeBasalShearFromSSA();
  virtual PetscErrorCode fillRegPoissonData(RegPoissonCtx &user);
  virtual PetscErrorCode computeTFAFromBasalShear(
                const PetscScalar phi_low, const PetscScalar phi_high,
                const PetscScalar invRegEps, const char *invfieldsfilename);

  // see iMIO.cc
  virtual PetscErrorCode set_time_from_options();
  virtual PetscErrorCode dumpToFile(const char *filename);
  virtual PetscErrorCode regrid();

  // see iMmatlab.cc
  virtual PetscErrorCode writeSSAsystemMatlab(IceModelVec2 vNuH[2]);

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
                       PetscScalar iarea, bool useHomoTemp, 
                       PetscScalar &gmeltfrac, PetscScalar &gtemp0);
  virtual PetscErrorCode ageStats(PetscScalar ivol, PetscScalar &gorigfrac);
  virtual PetscErrorCode summary(bool tempAndAge, bool useHomoTemp);
  virtual PetscErrorCode summaryPrintLine(
              PetscTruth printPrototype, bool tempAndAge,
              PetscScalar year, PetscScalar delta_t, 
              PetscScalar volume, PetscScalar area,
              PetscScalar meltfrac, PetscScalar H0, PetscScalar T0);

  // see iMreport.cc;  methods for computing diagnostic quantities:
  // spatially-varying:
  virtual PetscErrorCode compute_by_name(string name, IceModelVec* &result);
  virtual PetscErrorCode compute_bwp(IceModelVec2 &result);
  virtual PetscErrorCode compute_cbar(IceModelVec2 &result);
  virtual PetscErrorCode compute_cbase(IceModelVec2 &result, IceModelVec2 &tmp);
  virtual PetscErrorCode compute_cflx(IceModelVec2 &result, IceModelVec2 &cbar);
  virtual PetscErrorCode compute_csurf(IceModelVec2 &result, IceModelVec2 &tmp);
  virtual PetscErrorCode compute_dhdt(IceModelVec2 &result);
  virtual PetscErrorCode compute_hardav(IceModelVec2 &result);
  virtual PetscErrorCode compute_rank(IceModelVec2 &result);
  virtual PetscErrorCode compute_taud(IceModelVec2 &result, IceModelVec2 &tmp);
  virtual PetscErrorCode compute_temp_pa(IceModelVec3 &useForPATemp);
  virtual PetscErrorCode compute_tempbase(IceModelVec2 &result);
  virtual PetscErrorCode compute_temppabase(IceModelVec3 &hasPATemp,
                                            IceModelVec2 &result);
  virtual PetscErrorCode compute_tempsurf(IceModelVec2 &result);
  virtual PetscErrorCode compute_uvelbase(IceModelVec2 &result);
  virtual PetscErrorCode compute_uvelsurf(IceModelVec2 &result);
  virtual PetscErrorCode compute_vvelbase(IceModelVec2 &result);
  virtual PetscErrorCode compute_vvelsurf(IceModelVec2 &result);
  virtual PetscErrorCode compute_wvelbase(IceModelVec2 &result);
  virtual PetscErrorCode compute_wvelsurf(IceModelVec2 &result);
  // scalar:
  virtual PetscErrorCode ice_mass_bookkeeping();
  virtual PetscErrorCode compute_ice_volume(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area_grounded(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area_floating(PetscScalar &result);
  virtual PetscErrorCode compute_by_name(string name, PetscScalar &result);

  // see iMsia.cc
  virtual PetscErrorCode surfaceGradientSIA();
  virtual PetscErrorCode velocitySIAStaggered();
  virtual PetscErrorCode basalSlidingHeatingSIA();
  virtual PetscErrorCode velocities2DSIAToRegular();
  virtual PetscErrorCode SigmaSIAToRegular();
  virtual PetscErrorCode horizontalVelocitySIARegular();
  virtual PetscScalar    basalVelocitySIA( // not recommended, generally
                             PetscScalar x, PetscScalar y, PetscScalar H, PetscScalar T,
                             PetscScalar alpha, PetscScalar mu, PetscScalar min_T) const;

  // see iMssa.cc
  virtual PetscErrorCode initSSA();
  virtual PetscErrorCode velocitySSA(PetscInt *numiter);
  virtual PetscErrorCode velocitySSA(IceModelVec2 vNuH[2], PetscInt *numiter);
  virtual PetscErrorCode computeEffectiveViscosity(IceModelVec2 vNuH[2], PetscReal epsilon);
  virtual PetscErrorCode testConvergenceOfNu(IceModelVec2 vNuH[2], IceModelVec2 vNuHOld[2],
                                             PetscReal *norm, PetscReal *normChange);
  virtual PetscErrorCode assembleSSAMatrix(bool includeBasalShear, IceModelVec2 vNuH[2], Mat A);
  virtual PetscErrorCode assembleSSARhs(Vec rhs);
  virtual PetscErrorCode moveVelocityToDAVectors(Vec x);
  virtual PetscErrorCode broadcastSSAVelocity(bool updateVelocityAtDepth);
  virtual PetscErrorCode correctSigma();
  virtual PetscErrorCode correctBasalFrictionalHeating();

  // see iMssaSNES.cc; UNDER DEVELOPMENT
  virtual PetscErrorCode mapUVbarSSAToSSASNESVec(DA ssasnesda, Vec &ssasnesX);
  virtual PetscErrorCode mapSSASNESVecToUVbarSSA(DA ssasnesda, Vec ssasnesX);
  virtual PetscErrorCode setbdryvalSSA(DA ssasnesda, Vec &ssasnesBV);
  virtual PetscErrorCode solvefeedback(SNES snes, Vec residual);
  virtual PetscErrorCode getNuFromNuH(IceModelVec2 vNuH[2], SSASNESCtx *user);
  virtual PetscErrorCode velocitySSA_SNES(IceModelVec2 vNuH[2], PetscInt *its);

  // see iMtemp.cc; uses columnSystem.{hh|cc}
  virtual PetscErrorCode temperatureStep();
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
  Vec g2;			//!< Global work vector; used by the bed
				//!< deformation code only. (FIXME!)
  IceModelVec2 vWork2d[nWork2d];
  // 3D working space (with specific purposes)
  IceModelVec3 Tnew3, taunew3;
  IceModelVec3 Sigmastag3[2], Istag3[2];

  // for saving SSA velocities for restart
  int have_ssa_velocities;	//!< use vubarSSA and vvbarSSA from a previous
				//! run if 1, otherwise set them to zero in
				//! IceModel::initSSA()
  IceModelVec2 vubarSSA, vvbarSSA;

  // SSA solve vars; note pieces of the SSA Velocity routine are defined in iMssa.cc
  KSP SSAKSP;
  Mat SSAStiffnessMatrix;
  Vec SSAX, SSARHS;  // Global vectors for solution of the linear system
  Vec SSAXLocal; // We need a local copy of the solution to map back to a DA based vector
  VecScatter SSAScatterGlobalToLocal;

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

  // diagnostic viewers; see iMviewers.cc
  virtual PetscErrorCode init_viewers();
  virtual PetscErrorCode update_viewers();
  virtual PetscErrorCode update_nu_viewers(IceModelVec2 vNu[2], IceModelVec2[2], bool);
  set<string> map_viewers, slice_viewers, sounding_viewers;
  PetscInt     id, jd;	     // sounding indices
  bool view_diffusivity, view_log_nuH, view_nuH;

private:
  // for event logging (profiling); see run() and velocity()
  int siaEVENT, ssaEVENT, velmiscEVENT, beddefEVENT, massbalEVENT, tempEVENT;
};

#endif /* __iceModel_hh */


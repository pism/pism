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

#ifndef __iceModel_hh
#define __iceModel_hh

#include <signal.h>
#include <netcdf.h>
#include <gsl/gsl_rng.h>
#include <petscsnes.h>

#include "materials.hh"
#include "ssaJed/iceShelfExtension.hh"  // only used for Jed's ssa_external in base/ssa/
#include "pism_const.hh"
#include "grid.hh"
#include "beddefLC.hh"
#include "iceModelVec.hh"

#include "../coupler/forcing.hh"
#include "../coupler/pccoupler.hh"

// With this forward declaration, we don't have to recompile all of
// Pism if the SSA header "pismssa.hh" changes.
typedef struct _p_SSA *SSA;

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond
 

// see iMnames.cc
struct titleNname {
  char title[100]; // these short titles appear on PETSc graphical viewers and 
                   //   in Matlab output file
  char name[30];   // these names are for Matlab output vars
};

struct SSASNESNode {
  PetscScalar u, v;
};
struct SSASNESCtx {
  DA               ssada;
  IceGrid          *grid;
  PlasticBasalType *basal;
  IceModelVec2     ctxH,
                   ctxMask,
                   ctxtauc,
                   ctxtaudx,
                   ctxtaudy,
                   ctxNu[2];
  Vec              ctxBV;
  PetscScalar      schoofReg,
                   constantHardness;
  PetscTruth       useConstantHardness,
                   useConstantNu,
                   useStoredNu,
                   usePlasticBasalType;
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
  IceModel(IceGrid &g);
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
  virtual PetscErrorCode diagnosticRun();
  virtual PetscErrorCode setExecName(const char *my_executable_short_name);
  virtual IceFactory &getIceFactory() { return iceFactory; }
  virtual IceType *getIce() {return ice;}

  // see iMbootstrap.cc 
  virtual PetscErrorCode bootstrapFromFile(const char *fname);
  virtual PetscErrorCode readShelfStreamBCFromFile(const char *fname);

  // see iMoptions.cc
  virtual PetscErrorCode setFromOptions();
  
  // see iMtests.cc; FIXME: these should go in a derived class
  virtual PetscErrorCode testIceModelVec3();
  virtual PetscErrorCode testIceModelVec3Bedrock();

  // see iMutil.cc
  virtual PetscErrorCode attachAtmospherePCC(PISMAtmosphereCoupler &aPCC);
  virtual PetscErrorCode attachOceanPCC(PISMOceanCoupler &oPCC);
  virtual PetscErrorCode additionalAtStartTimestep();
  virtual PetscErrorCode additionalAtEndTimestep();

  // see iMIO.cc
  bool hasSuffix(const char* fname, const char* suffix) const;  
  virtual PetscErrorCode initFromFile(const char *);
  virtual PetscErrorCode writeFiles(const char* default_filename, 
                                    const PetscTruth forceFullDiagnostics = PETSC_FALSE);
  virtual PetscErrorCode write_model_state(const char filename[]);
  virtual PetscErrorCode write_extra_fields(const char filename[]);

protected:

  IceGrid               &grid;
  int initial_Mz; //!< the number of vertical grid levels at the start of the
		  //! run; used by the grid extension code

  PolarStereoParams     psParams;
  NCConfigVariable      config;
  
  LocalInterpCtx        *bootstrapLIC;

  IceFactory            iceFactory;
  IceType               *ice;
  PlasticBasalType      *basal;
  BasalTypeSIA          *basalSIA;
  BedrockThermalType    bed_thermal;
  DeformableEarthType   bed_deformable;
  SeaWaterType          ocean;
  FreshWaterType        porewater;
  SSAStrengthExtension  ssaStrengthExtend;

  PISMAtmosphereCoupler *atmosPCC;
  IceInfoNeededByAtmosphereCoupler info_atmoscoupler;

  PISMOceanCoupler      *oceanPCC;
  IceInfoNeededByOceanCoupler info_oceancoupler;

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
        vMask,		//!< mask for flow type with values SHEET, DRAGGING, FLOATING
        vGhf,		//!< geothermal flux
        vRb,		//!< basal frictional heating on regular grid
        vtillphi,	//!< friction angle for till under grounded ice sheet
        vuvbar[2],	//!< ubar and vbar on staggered grid; ubar at i+1/2, vbar at j+1/2
        vub, vvb,	//!< basal velocities on standard grid
        vubar, vvbar;	//!< vertically-averaged horizontal velocity on standard grid

  IceModelVec3
        u3, v3, w3,	//!< velocity of ice; m s-1
        Sigma3, 	//!< strain-heating term in conservation of energy model; J s-1 m-3
        T3,		//!< absolute temperature of ice; K
        tau3;		//!< age of ice; s

  IceModelVec3Bedrock
        Tb3;		//!< temperature of lithosphere (bedrock) under ice or ocean; K

  // parameters
  PetscReal   maxdt, muSliding, enhancementFactor, initial_age_years_default,
              dt, dtTempAge,  // current mass cont. and temp/age; time steps in seconds
              dt_force,
              ssaRelativeTolerance, ssaEpsilon, betaShelvesDragToo,
              plastic_till_c_0, plastic_till_mu, plastic_till_pw_fraction, plasticRegularization,
              tauc_default_value, pseudo_plastic_q, pseudo_plastic_uthreshold,
              startYear, endYear, run_year_default, maxdt_temporary,
              ssaIntervalYears, bedDefIntervalYears, adaptTimeStepRatio,
              CFLviolcount,    //!< really is just a count, but PetscGlobalSum requires this type
              dt_from_diffus, dt_from_cfl, CFLmaxdt, CFLmaxdt2D, gDmax,
              gmaxu, gmaxv, gmaxw,  // global maximums on 3D grid of abs value of vel components
              gdHdtav,  //!< average value in map-plane (2D) of dH/dt, where there is ice; m s-1
              dvoldt,  //!< d(total ice volume)/dt; m3 s-1
              min_temperature_for_SIA_sliding, Hmelt_max, globalMinAllowedTemp,
              constantGrainSize;
  PetscInt    skipCountDown,
              skipMax,
              noSpokesLevel,
              maxLowTempCount,
              ssaMaxIterations;

  // flags
  PetscTruth  doMassConserve, doTemp, doBedDef, doBedIso,
              thermalBedrock, includeBMRinContinuity, updateHmelt,
              isDrySimulation, holdTillYieldStress, useConstantTillPhi,
              useSSAVelocity, doPlasticTill, doPseudoPlasticTill,
              doSuperpose, useConstantNuHForSSA,
              computeSurfGradInwardSSA,
              leaveNuHAloneSSA, shelvesDragToo,
              yearsStartRunEndDetermined, doAdaptTimeStep, doOceanKill, floatingIceKilled,
              realAgeForGrainSize,
              showViewers, ssaSystemToASCIIMatlab, doSkip, reportHomolTemps,
              allowAboveMelting,
              createVecs_done, createViewers_done, createBasal_done,
              computeSIAVelocities, transformForSurfaceGradient;
  char        adaptReasonFlag;

  // file names
  char         ssaMatlabFilePrefix[PETSC_MAX_PATH_LEN];

  // viewer and sounding
  static const PetscInt   tnN = 75;
  static const            titleNname tn[tnN];  // see iMnames.cc
  PetscViewer             runtimeViewers[tnN]; // see iMviewers.cc
  char         diagnostic[PETSC_MAX_PATH_LEN], // see iMviewers.cc
               diagnosticBIG[PETSC_MAX_PATH_LEN], // see iMviewers.cc
               matlabOutVars[PETSC_MAX_PATH_LEN]; // see iMmatlab.cc
  PetscDrawLG  kspLG;
  PetscInt     id, jd, kd;
  Vec          Td, wd, ud, vd, Sigmad, taud; // hold soundings

  string history; // history of commands used to generate this instance of IceModel
  string executable_short_name;
  
protected:
  // see iceModel.cc
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode destroyVecs();
  virtual void setMaxTimeStepYears(PetscScalar); //!< use this value for adaptive stepping
  virtual void setAdaptTimeStepRatio(PetscScalar);
  virtual PetscErrorCode setStartYear(PetscScalar);
  virtual PetscErrorCode setEndYear(PetscScalar);
  virtual void setInitialAgeYears(PetscScalar d);
  virtual void setAllGMaxVelocities(PetscScalar);
  virtual void setConstantNuHForSSA(PetscScalar);

  // see iMadaptive.cc
  virtual PetscErrorCode computeMaxDiffusivity(bool updateDiffusViewer);
  virtual PetscErrorCode computeMax3DVelocities();
  virtual PetscErrorCode computeMax2DSlidingSpeed();
  virtual PetscErrorCode adaptTimeStepDiffusivity();
  virtual PetscErrorCode determineTimeStep(const bool doTemperatureCFL);

  // see iMbasal.cc
  virtual PetscErrorCode initBasalTillModel();
  virtual PetscErrorCode computePhiFromBedElevation();
  virtual PetscScalar    getEffectivePressureOnTill(
            const PetscScalar thk, const PetscScalar melt_thk);
  virtual PetscErrorCode updateYieldStressFromHmelt();
  virtual PetscErrorCode diffuseHmelt();
  virtual PetscScalar    basalVelocity(PetscScalar x, PetscScalar y,
                                       PetscScalar H, PetscScalar T, PetscScalar alpha, PetscScalar mu) const;
  virtual PetscScalar basalDragx(PetscScalar **tauc, PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;
  virtual PetscScalar basalDragy(PetscScalar **tauc, PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;

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
  virtual PetscErrorCode bootstrapSetBedrockColumnTemp(const PetscInt i, const PetscInt j,
                            const PetscScalar Ttopbedrock, const PetscScalar geothermflux);
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
  virtual PetscErrorCode write3DPlusToFile(const char filename[]);
  virtual PetscErrorCode regrid();

  // see iMmatlab.cc
  virtual bool           matlabOutWanted(const char name);
  virtual PetscErrorCode VecView_g2ToMatlab(PetscViewer v, 
                                    const char *varname, const char *shorttitle);
  virtual PetscErrorCode write2DToMatlab(PetscViewer v, const char singleCharName, 
                                 IceModelVec2 &l2, const PetscScalar scale);
  virtual PetscErrorCode write2DToMatlab(PetscViewer v, const char singleCharName, 
                                 Vec l2, const PetscScalar scale);
  virtual PetscErrorCode writeSliceToMatlab(PetscViewer v, const char singleCharName, 
                                    IceModelVec3 &imv3, const PetscScalar scale);
  virtual PetscErrorCode writeSurfaceValuesToMatlab(PetscViewer v, const char singleCharName, 
                                            IceModelVec3 &imv3, const PetscScalar scale);
  virtual PetscErrorCode writeSpeed2DToMatlab(PetscViewer v, const char scName, 
                          IceModelVec2 &lu, IceModelVec2 &lv, const PetscScalar scale, 
                          const PetscTruth doLog, const PetscScalar log_missing);
  virtual PetscErrorCode writeSpeedSurfaceValuesToMatlab(PetscViewer v, const char scName, 
                          IceModelVec3 &imv3_u, IceModelVec3 &imv3_v, const PetscScalar scale, 
                          const PetscTruth doLog, const PetscScalar log_missing);
  virtual PetscErrorCode writeLog2DToMatlab(PetscViewer v, const char scName, 
                          IceModelVec2 &l, const PetscScalar scale, const PetscScalar thresh,
                          const PetscScalar log_missing);
  virtual PetscErrorCode writeSoundingToMatlab(PetscViewer v, const char scName, 
                          IceModelVec3 &imv3,
                          const PetscScalar scale, const PetscTruth doTandTb);
  virtual PetscErrorCode writeMatlabVars(const char *fname);
  virtual PetscErrorCode writeSSAsystemMatlab(IceModelVec2 vNuH[2]);

  // see iMnames.cc; note tn is statically-initialized in iMnames.cc
  int cIndex(const char singleCharName);

  // see iMreport.cc
  virtual PetscErrorCode computeFlowUbarStats
                      (PetscScalar *gUbarmax, PetscScalar *gUbarSIAav,
                       PetscScalar *gUbarstreamav, PetscScalar *gUbarshelfav,
                       PetscScalar *gicegridfrac, PetscScalar *gSIAgridfrac,
                       PetscScalar *gstreamgridfrac, PetscScalar *gshelfgridfrac);
  virtual PetscErrorCode volumeArea(PetscScalar& gvolume,PetscScalar& garea,
                            PetscScalar& gvolSIA, PetscScalar& gvolstream, 
                            PetscScalar& gvolshelf);
  virtual PetscErrorCode summary(bool tempAndAge, bool useHomoTemp);
  virtual PetscErrorCode summaryPrintLine(
              const PetscTruth printPrototype, const PetscTruth tempAndAge,
              const PetscScalar year, const PetscScalar dt, 
              const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
              const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0);

  // see iMsia.cc
  virtual PetscErrorCode surfaceGradientSIA();
  virtual PetscErrorCode velocitySIAStaggered();
  virtual PetscErrorCode basalSlidingHeatingSIA();
  virtual PetscErrorCode velocities2DSIAToRegular();
  virtual PetscErrorCode SigmaSIAToRegular();
  virtual PetscErrorCode horizontalVelocitySIARegular();

  // see iMssa.cc
  virtual PetscErrorCode initSSA();
  virtual PetscErrorCode velocitySSA(PetscInt *numiter);
  virtual PetscErrorCode velocitySSA(IceModelVec2 vNuH[2], PetscInt *numiter);
  virtual PetscErrorCode computeEffectiveViscosity(IceModelVec2 vNuH[2], PetscReal epsilon);
  virtual PetscErrorCode testConvergenceOfNu(IceModelVec2 vNuH[2], IceModelVec2 vNuHOld[2],
                                             PetscReal *norm, PetscReal *normChange);
  virtual PetscErrorCode assembleSSAMatrix(const bool includeBasalShear, IceModelVec2 vNuH[2], Mat A);
  virtual PetscErrorCode assembleSSARhs(bool surfGradInward, Vec rhs);
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
  virtual PetscErrorCode temperatureAgeStep();
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);
  virtual PetscErrorCode ageStep(PetscScalar* CFLviol);
  virtual bool checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                      PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE);
  virtual PetscErrorCode excessToFromBasalMeltLayer(
                      const PetscScalar rho_c, const PetscScalar z, const PetscScalar dz,
                      PetscScalar *Texcess, PetscScalar *Hmelt);

  // see iMutil.cc
  virtual PetscErrorCode getMagnitudeOf2dVectorField(IceModelVec2 &vfx, IceModelVec2 &vfy,
						     IceModelVec2 &vmag);
  virtual int endOfTimeStepHook();
  virtual PetscErrorCode stampHistoryCommand();
  virtual PetscErrorCode stampHistoryEnd();
  virtual PetscErrorCode stampHistory(string);
  virtual PetscErrorCode stampHistoryAdd(string);
  virtual PetscErrorCode check_maximum_thickness();

  // see iMvelocity.cc
  virtual PetscErrorCode velocity(bool updateSIAVelocityAtDepth);    
  virtual PetscErrorCode vertVelocityFromIncompressibility();
    
  // see iMviewers.cc
  int isViewer(char name);
  virtual PetscErrorCode updateSoundings();
  virtual PetscErrorCode updateOneSounding(const char scName, Vec l, const PetscScalar scale);
  virtual PetscErrorCode getViewerDims(const PetscInt target_size, 
                               const PetscScalar Lx, const PetscScalar Ly,
                               PetscInt *xdim, PetscInt *ydim);
  virtual PetscErrorCode createOneViewerIfDesired(const char singleCharName);
  virtual PetscErrorCode createOneViewerIfDesired(const char singleCharName, const char* title);
  virtual PetscErrorCode createOneViewerIfDesired(PetscViewer* v, 
                                          const char singleCharName, const char* title);
  virtual PetscErrorCode createViewers();
  virtual PetscErrorCode update2DViewer(const char scName, Vec l2, const PetscScalar scale);
  virtual PetscErrorCode update2DViewer(const char scName, IceModelVec2 &l2, const PetscScalar scale);
  virtual PetscErrorCode updateSliceViewer(const char scName, IceModelVec3 &imv3, const PetscScalar scale);
  virtual PetscErrorCode updateSurfaceValuesViewer(const char scName, 
                   IceModelVec3 &imv3, const PetscScalar scale);
  virtual PetscErrorCode updateSpeed2DViewer(const char scName, IceModelVec2 &lu, IceModelVec2 &lv, 
                   const PetscScalar scale, const PetscTruth doLog, 
                   const PetscScalar log_missing);
  virtual PetscErrorCode updateSpeedSurfaceValuesViewer(const char scName, 
                   IceModelVec3 &imv3_u, IceModelVec3 &imv3_v, 
                   const PetscScalar scale, const PetscTruth doLog, 
                   const PetscScalar log_missing);
  virtual PetscErrorCode updateLog2DViewer(const char scName, Vec l,
                   const PetscScalar scale, const PetscScalar thresh, 
                   const PetscScalar log_missing);
  virtual PetscErrorCode updateViewers();  // it calls updateSoundings()
  virtual PetscErrorCode updateNuViewers(IceModelVec2 vNu[2], IceModelVec2 vNuOld[2], bool updateNu_tView);
  virtual PetscErrorCode destroyViewers();

protected:
  // working space (a convenience)
  static const PetscInt nWork2d=6;
  Vec g2;			//!< Global work vector
  IceModelVec2 vWork2d[nWork2d];
  // 3D working space (with specific purposes)
  IceModelVec3 Tnew3, taunew3;
  IceModelVec3 Sigmastag3[2], Istag3[2];

protected:
  int have_ssa_velocities;	//!< use vubarSSA and vvbarSSA from a previous
				//! run if 1, otherwise set them to zero in
				//! IceModel::initSSA()
  IceModelVec2 vubarSSA, vvbarSSA;

private:
  // for event logging (profiling); see run() and velocity()
  int siaEVENT, ssaEVENT, velmiscEVENT, beddefEVENT, massbalEVENT, tempEVENT;

protected:
  // Pieces of the SSA Velocity routine defined in iMssa.cc.
  KSP SSAKSP;
  Mat SSAStiffnessMatrix;
  Vec SSAX, SSARHS;  // Global vectors for solution of the linear system
  Vec SSAXLocal; // We need a local copy of the solution to map back to a DA based vector
  VecScatter SSAScatterGlobalToLocal;

  // Jed's External SSA solver:  "If non-NULL, we are using the velocity solver in
  // src/base/ssa and all the SSA* objects above are not used.  Eventually we should
  // move SSA legacy stuff out of IceModel (either by putting into an implementation
  // of SSA or by just deleting it)."
  // Bueler's comments:  The type "SSA" is not documented and is defined as opaquely
  // as possible by going through src/ssaJed/pismssa.h and src/ssaJed/ssaimpl.h.  This seems
  // to be Jed's decision to start using PIMPL instead of C++ like the rest of PISM.
  // The actual implementation of this SSA is a finite element method is in 
  // src/ssaJed/ssafe.c.  It is very promising.  It is essentially impossible for me to
  // understand without an analysis of the design principles of dohp.  Go to dohp.org
  // and read Karniadakis and Sherwin.
  SSA ssa;

  // Jed's shelfExtension object: a strength extension that factors the nu*H coefficient
  // of the SSA equations so that it can use your IceType.  Essentially unnecessary because
  // the extension is the replacement for an unimplemented stress boundary condition.
  IceShelfExtension  shelfExtensionJed;  
                                      

  // This is related to the snapshot saving feature
  char snapshots_filename[PETSC_MAX_PATH_LEN];
  bool save_snapshots, file_is_ready, split_snapshots;

  // equally spaced snapshots
  PetscTruth save_at_equal_intervals;
  double first_snapshot, snapshot_dt, next_snapshot, last_snapshot;

  // unequally spaced snapshots
  static const int max_n_snapshots = 200;
  double save_at[max_n_snapshots];
  int n_snapshots, current_snapshot;

  PetscErrorCode init_snapshots_from_options();
  PetscErrorCode write_snapshot();
};

#endif /* __iceModel_hh */


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

#ifndef __iceModel_hh
#define __iceModel_hh

#include <signal.h>
#include <netcdf.h>
#include <gsl/gsl_rng.h>
#include <petscda.h>
#include <petscksp.h>
#include "materials.hh"
#include "pism_const.hh"
#include "grid.hh"
#include "forcing.hh"
#include "beddefLC.hh"
#include "iceModelVec.hh"

// remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond
 

// see iMnames.cc
struct titleNname {
  char title[100]; // these short titles appear on PETSc graphical viewers and 
                   //   in Matlab output file
  char name[30];   // these names are for Matlab output vars
};


// see iMIO.cc
struct PolarStereoParams {
  // these are "double" and not "float" ultimately because of how ncgen works
  double svlfp, // straight_vertical_longitude_from_pole; defaults to 0
         lopo,  // latitude_of_projection_origin; defaults to 90
         sp;    // standard_parallel; defaults to -71
};


// see iMinverseMat.cc
struct RegPoissonTaucCtx {
  // describes Poisson-like problem solved when inverting surface velocities
  //   using a regularization
  DA          da;             // must be first in struct
  IceGrid     *grid;
  PetscReal   epsilon;
  Vec         f;              // = f(x,y) in PDE
  Vec         g;              // = g(x,y) in PDE
};
extern PetscErrorCode RegPoissonTaucFunctionLocal(
                DALocalInfo *info, PetscScalar **x, PetscScalar **F,
                RegPoissonTaucCtx *user);
extern PetscErrorCode RegPoissonTaucJacobianLocal(
                DALocalInfo *info, PetscScalar **x, Mat jac,
                RegPoissonTaucCtx *user);

//! The base class for PISM.  Contains all essential variables, parameters, and flags for modelling an ice sheet.
class IceModel {
public:
  // see iceModel.cc for implementation of constructor and destructor:
  IceModel(IceGrid &g, IceType *i);
  virtual ~IceModel(); // must be virtual merely because some members are virtual

  // see iceModel.cc
  virtual PetscErrorCode run();
  virtual PetscErrorCode diagnosticRun();
  virtual PetscErrorCode additionalAtStartTimestep();
  virtual PetscErrorCode additionalAtEndTimestep();
  virtual PetscErrorCode setExecName(const char *my_executable_short_name);

  // see iMbootstrap.cc 
  virtual PetscErrorCode bootstrapFromFile(const char *fname);
  virtual PetscErrorCode readShelfStreamBCFromFile(const char *fname);

  // see iMinverse.cc
  virtual PetscErrorCode invertSurfaceVelocities();

  // see iMoptions.cc
  virtual PetscErrorCode setFromOptions();

  // see iMtests.cc; FIXME: these should go in a derived class
  virtual PetscErrorCode testIceModelVec3();
  virtual PetscErrorCode testIceModelVec3Bedrock();

  // see iMutil.cc
  virtual PetscErrorCode initFromOptions(PetscTruth doHook = PETSC_TRUE);

  // see iMIO.cc
  bool hasSuffix(const char* fname, const char* suffix) const;  
  virtual PetscErrorCode initFromFile(const char *);
  virtual PetscErrorCode writeFiles(const char* default_filename, 
                                    const PetscTruth forceFullDiagnostics = PETSC_FALSE);
  virtual PetscErrorCode write_model_state(const char filename[]);

protected:
   static const int MASK_SHEET;
   static const int MASK_DRAGGING;
   static const int MASK_FLOATING;
   static const int MASK_FLOATING_OCEAN0;

  IceGrid            &grid;

  PolarStereoParams  psParams;
  
  LocalInterpCtx     *bootstrapLIC;

  IceType               *ice;
  PlasticBasalType      *basal;
  BasalTypeSIA          *basalSIA;
  BedrockThermalType    bed_thermal;
  DeformableEarthType   bed_deformable;
  SeaWaterType          ocean;
  FreshWaterType        porewater;

  // state variables
  IceModelVec2 vh,		//!< ice surface elevation
    vH,				//!< ice thickness
    vdHdt,			//!< \frac{dH}{dt}
    vtauc,			//!< yield stress for basal till (plastic or pseudo-plastic model)
    vHmelt, 			//!< thickness of the basal meltwater
    vbasalMeltRate,		//!< basal meltwater production rate
  /*!< rate of production of basal meltwater (ice-equivalent) */
    vLongitude,		   //!< Longitude
    vLatitude,		   //!< Latitude 
    vbed,		   //!< bed topography
    vuplift,		   //!< bed uplift rate
    vub, vvb,		   //!< basal velocities on standard grid
    vuvbar[2],		   //!< u bar and v bar on staggered grid
  /*!< 2D; vuvbar[0] and vuvbar[1] are u bar and v bar on staggered grid, */
    vubar, vvbar,	     	//!< vertically-averaged vels on standard grid
  /*!< vertically-averaged vels (u bar and v bar) on standard grid */
    vMask,		   //!< mask for flow type
    vTs,		   //!< surface temperature
    vAccum,		   //!< accumulation
    vAccumSnow,		   //!< see iMpdd.cc; vAccum is net (including ablation)
    vGhf,		   //!< geothermal flux
    vRb,		   //!< basal frictional heating on regular grid
    vtillphi;		   //!< friction angle for till under grounded ice sheet

  IceModelVec3        u3, v3, w3, Sigma3, T3, tau3;  // note gs3 was one of these
  IceModelVec3Bedrock Tb3;

  // parameters
  PetscScalar maxdt, muSliding, enhancementFactor, initial_age_years_default;
  PetscScalar dt, dtTempAge, dt_force; // current mass cont. and temp/age 
                                       //   time steps in seconds
  PetscScalar constantNuHForSSA, constantHardnessForSSA, min_thickness_SSA,
              regularizingVelocitySchoof, regularizingLengthSchoof,
              ssaRelativeTolerance, ssaEpsilon, betaShelvesDragToo;
  PetscScalar plastic_till_c_0, plastic_till_mu, plastic_till_pw_fraction, plasticRegularization,
              tauc_default_value, pseudo_plastic_q, pseudo_plastic_uthreshold;
  PetscScalar startYear, endYear, run_year_default, maxdt_temporary;
  PetscScalar ssaIntervalYears, bedDefIntervalYears, adaptTimeStepRatio;
  PetscScalar CFLviolcount;    // really is just a count, but PetscGlobalSum requires this type
  PetscScalar dt_from_diffus, dt_from_cfl, CFLmaxdt, CFLmaxdt2D, gDmax;
  PetscScalar gmaxu, gmaxv, gmaxw;  // global maximums on 3D grid for abs value 
                                    // of 3D components of velocities
  PetscScalar gdHdtav, dvoldt; // average value in map-plane (2D) of dH/dt (where there is ice) 
                               //   [units m/s] and d(volume)/dt [units m^3/s]
  PetscScalar min_temperature_for_SIA_sliding, Hmelt_max, globalMinAllowedTemp, 
              oceanHeatFlux, seaLevel;
  PetscScalar constantGrainSize;
  PetscInt    skipCountDown, skipMax, noSpokesLevel, maxLowTempCount,
              flowLawNumber,ssaMaxIterations;

  // flags
  PetscTruth  doMassConserve, doTemp, doBedDef, doBedIso;
  PetscTruth  initialized_p, thermalBedrock, includeBMRinContinuity, updateHmelt,
              isDrySimulation, holdTillYieldStress, useConstantTillPhi;
  PetscTruth  useSSAVelocity, doPlasticTill, doPseudoPlasticTill,
              doSuperpose, useConstantNuHForSSA, 
              useConstantHardnessForSSA, computeSurfGradInwardSSA,
              leaveNuHAloneSSA, shelvesDragToo;
  PetscTruth  yearsStartRunEndDetermined, doAdaptTimeStep, doOceanKill, floatingIceKilled;
  PetscTruth  flowLawUsesGrainSize, realAgeForGrainSize;
  PetscTruth  showViewers, ssaSystemToASCIIMatlab, doSkip, reportHomolTemps,
              allowAboveMelting;
  PetscTruth  createVecs_done, createViewers_done, createBasal_done;
  PetscTruth  computeSIAVelocities, transformForSurfaceGradient;
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
  Vec          Td, wd, ud, vd, Sigmad, taud;

  char *history; // history of commands used to generate this 
		 // instance of IceModel
  size_t history_size;	// the number of bytes allocated for the
			// history string
  char executable_short_name[PETSC_MAX_PATH_LEN];
  
protected:
  // see iceModel.cc
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode destroyVecs();
  virtual void setMaxTimeStepYears(PetscScalar); // use this value for adaptive stepping
  virtual void setAdaptTimeStepRatio(PetscScalar);
  virtual PetscErrorCode setStartYear(PetscScalar);
  virtual PetscErrorCode setEndYear(PetscScalar);
  virtual void setInitialAgeYears(PetscScalar d);
  virtual void setAllGMaxVelocities(PetscScalar);
  virtual void setConstantNuHForSSA(PetscScalar);
  virtual PetscTruth isInitialized() const;
  virtual int intMask(PetscScalar);
  virtual int modMask(PetscScalar);

  // see iMadaptive.cc
  virtual PetscErrorCode computeMaxDiffusivity(bool updateDiffusViewer);
  virtual PetscErrorCode computeMax3DVelocities();
  virtual PetscErrorCode computeMax2DSlidingSpeed();
  virtual PetscErrorCode adaptTimeStepDiffusivity();
  virtual PetscErrorCode determineTimeStep(const bool doTemperatureCFL);

  // see iMbasal.cc
  virtual PetscErrorCode initBasalTillModel();
  virtual PetscScalar    getEffectivePressureOnTill(
            const PetscScalar thk, const PetscScalar melt_thk);
  virtual PetscErrorCode updateYieldStressFromHmelt();
  virtual PetscScalar    basalVelocity(const PetscScalar x, const PetscScalar y, 
            const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
            const PetscScalar mu);
  virtual PetscScalar basalDragx(PetscScalar **tauc, PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;
  virtual PetscScalar basalDragy(PetscScalar **tauc, PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;

  // see iMbeddef.cc: possibly useful general tool for putting Vecs on processor zero
  // (e.g. also used by derived class IceHEINOModel)
  Vec            g2natural;
  VecScatter     top0ctx;
  PetscTruth     top0ctx_created;
  PetscErrorCode createScatterToProcZero(Vec& samplep0);
  PetscErrorCode destroyScatterToProcZero();
  PetscErrorCode putLocalOnProcZero(Vec& vlocal, Vec& onp0);
  PetscErrorCode getLocalFromProcZero(Vec& onp0, Vec& vlocal);

  // see iMbeddef.cc
  BedDeformLC    bdLC;
  PetscScalar    lastBedDefUpdateYear;
  IceModelVec2   vbedlast;
  IceModelVec2   vHlast;	// used for simple pointwise isostasy and to
				// compute uplift
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

  // see iMforcing.cc
  IceSheetForcing  *dTforcing, *dSLforcing;
  PetscScalar    TsOffset;
  virtual PetscErrorCode initForcingFromOptions();
  virtual PetscErrorCode updateForcing();
  virtual PetscErrorCode forcingCleanup();

  // see iMgeometry.cc
  virtual PetscErrorCode computeDrivingStress(IceModelVec2 vtaudx, IceModelVec2 vtaudy);
  virtual PetscErrorCode updateSurfaceElevationAndMask();
  virtual PetscErrorCode massContExplicitStep();

  // see iMgrainsize.cc
//  PetscErrorCode computeGrainSize_PseudoAge(const PetscScalar H, const PetscInt Mz, 
//                          PetscScalar *w, PetscScalar *age_wspace, PetscScalar **gs);
  virtual PetscScalar    grainSizeVostok(PetscScalar age) const;

  // see iMinverse.cc
  virtual PetscErrorCode computeSIASurfaceVelocity(
	        IceModelVec2 &usSIA_out, IceModelVec2 &vsSIA_out);
  virtual PetscErrorCode getGforInverse(
                const PetscScalar x, const PetscScalar UsuSIAdiffsqr, 
                const PetscScalar UsuSIAdiffdotuSIA, const PetscScalar uSIAsqr,
                PetscScalar &G, PetscScalar &Gprime);
  virtual PetscErrorCode computeFofVforInverse(
                IceModelVec2 us_in, IceModelVec2 vs_in, 
                IceModelVec2 usSIA_in, IceModelVec2 vsSIA_in, 
                IceModelVec2 &fofv_out);
  virtual PetscErrorCode removeSIApart(
                IceModelVec2 us_in, IceModelVec2 vs_in, 
                IceModelVec2 usSIA_in, IceModelVec2 vsSIA_in, IceModelVec2 fofv_in,
		IceModelVec2 &ub_out, IceModelVec2 &vb_out);
  virtual PetscErrorCode computeBasalShearFromSSA(
                IceModelVec2 ub_in, IceModelVec2 vb_in, 
                IceModelVec2 &taubx_out, IceModelVec2 &tauby_out);
  virtual PetscErrorCode computeYieldStressFromBasalShearUsingPseudoPlastic(
                const PetscScalar invRegEps, const PetscTruth invShowFG,
                IceModelVec2 ub_in, IceModelVec2 vb_in,
	        IceModelVec2 taubx_in, IceModelVec2 tauby_in, 
                IceModelVec2 &tauc_out);
  virtual PetscErrorCode fillRegPoissonTaucData(
                IceModelVec2 ub_in, IceModelVec2 vb_in,
	        IceModelVec2 taubx_in, IceModelVec2 tauby_in,
	        RegPoissonTaucCtx &user);
  virtual PetscErrorCode computeTFAFromYieldStress(
                const PetscScalar phi_low, const PetscScalar phi_high,
                IceModelVec2 tauc_in, IceModelVec2 &tfa_out);
  virtual PetscErrorCode writeInvFields(const char *filename,
                IceModelVec2 us, IceModelVec2 vs,
                IceModelVec2 usSIA, IceModelVec2 vsSIA,
                IceModelVec2 taubxComputed, IceModelVec2 taubyComputed,
                IceModelVec2 fofv, IceModelVec2 taucComputed);

  // see iMIO.cc
  virtual PetscErrorCode warnUserOptionsIgnored(const char *fname);
  virtual PetscErrorCode setStartRunEndYearsFromOptions(const PetscTruth grid_p_year_VALID);
  virtual PetscErrorCode dumpToFile(const char *filename);
  virtual PetscErrorCode write3DPlusToFile(const char filename[]);
  virtual PetscErrorCode regrid(const char *fname);

  // see iMmatlab.cc
  virtual bool           matlabOutWanted(const char name);
  virtual PetscErrorCode VecView_g2ToMatlab(PetscViewer v, 
                                    const char *varname, const char *shorttitle);
  virtual PetscErrorCode write2DToMatlab(PetscViewer v, const char singleCharName, 
                                 IceModelVec2 l2, const PetscScalar scale);
  virtual PetscErrorCode write2DToMatlab(PetscViewer v, const char singleCharName, 
                                 Vec l2, const PetscScalar scale);
  virtual PetscErrorCode writeSliceToMatlab(PetscViewer v, const char singleCharName, 
                                    IceModelVec3 imv3, const PetscScalar scale);
  virtual PetscErrorCode writeSurfaceValuesToMatlab(PetscViewer v, const char singleCharName, 
                                            IceModelVec3 imv3, const PetscScalar scale);
  virtual PetscErrorCode writeSpeed2DToMatlab(PetscViewer v, const char scName, 
                          IceModelVec2 lu, IceModelVec2 lv, const PetscScalar scale, 
                          const PetscTruth doLog, const PetscScalar log_missing);
  virtual PetscErrorCode writeSpeedSurfaceValuesToMatlab(PetscViewer v, const char scName, 
                          IceModelVec3 imv3_u, IceModelVec3 imv3_v, const PetscScalar scale, 
                          const PetscTruth doLog, const PetscScalar log_missing);
  virtual PetscErrorCode writeLog2DToMatlab(PetscViewer v, const char scName, 
                          IceModelVec2 l, const PetscScalar scale, const PetscScalar thresh,
                          const PetscScalar log_missing);
  virtual PetscErrorCode writeSoundingToMatlab(PetscViewer v, const char scName, 
                          IceModelVec3 imv3,
                          const PetscScalar scale, const PetscTruth doTandTb);
  virtual PetscErrorCode writeMatlabVars(const char *fname);
  virtual PetscErrorCode writeSSAsystemMatlab(IceModelVec2 vNuH[2]);

  // see iMnames.cc; note tn is statically-initialized in iMnames.cc
  int cIndex(const char singleCharName);

  // see iMoptions.cc
  virtual PetscErrorCode determineSpacingTypeFromOptions(
                      const PetscTruth forceEqualIfNoOption);

  // see iMpdd.cc (positive degree day model)
  gsl_rng     *pddRandGen;
  char        monthlyTempsFile[PETSC_MAX_PATH_LEN];
  Vec*        vmonthlyTs;                  // will allocate 12 2D Vecs 
  PetscTruth  doPDD, doPDDTrueRand, pddStuffCreated, pddRandStuffCreated,
              pddHaveMonthlyTempData;
  PetscScalar pddStdDev, pddFactorSnow, pddFactorIce, pddRefreezeFrac, 
              pddSummerWarming, pddSummerPeakDay;
  virtual PetscErrorCode initPDDFromOptions();
  virtual PetscErrorCode readMonthlyTempDataPDD();
  virtual PetscErrorCode updateSurfaceBalanceFromPDD();
  virtual PetscErrorCode putBackSnowAccumPDD();
  virtual PetscErrorCode PDDCleanup();
  virtual double         CalovGreveIntegrand(const double Tac);
  virtual PetscScalar getTemperatureFromYearlyCycle(
                  const PetscScalar summer_warming, const PetscScalar Tma, 
                  const PetscScalar day, const PetscInt i, const PetscInt j) const;
  virtual PetscScalar getSummerWarming(
                  const PetscScalar elevation, const PetscScalar latitude, 
                  const PetscScalar Tma) const;
  virtual double getSurfaceBalanceFromSnowAndPDD(
                     const double snowrate, const double mydt, const double pdds);

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
  virtual PetscErrorCode testConvergenceOfNu(IceModelVec2 vNuH[2], IceModelVec2 vNuHOld[2], PetscReal *, PetscReal *);
  virtual PetscErrorCode assembleSSAMatrix(const bool includeBasalShear, IceModelVec2 vNuH[2], Mat A);
  virtual PetscErrorCode assembleSSARhs(bool surfGradInward, Vec rhs);
  virtual PetscErrorCode moveVelocityToDAVectors(Vec x);
  virtual PetscErrorCode broadcastSSAVelocity(bool updateVelocityAtDepth);
  virtual PetscErrorCode correctSigma();
  virtual PetscErrorCode correctBasalFrictionalHeating();

  // see iMtemp.cc
  virtual PetscErrorCode temperatureAgeStep();
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount);
  virtual PetscErrorCode ageStep(PetscScalar* CFLviol);
  virtual PetscErrorCode diffuseHmelt();
  virtual PetscErrorCode solveTridiagonalSystem(
           const PetscScalar* L, const PetscScalar* D, const PetscScalar* U,
           PetscScalar* x, const PetscScalar* rhs, PetscScalar* work, const int n) const;
  virtual bool checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                      PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE);
  virtual PetscErrorCode excessToFromBasalMeltLayer(
                      const PetscScalar rho_c, const PetscScalar z, const PetscScalar dz,
                      PetscScalar *Texcess, PetscScalar *Hmelt);
  virtual PetscErrorCode getMzMbzForTempAge(PetscInt &ta_Mz, PetscInt &ta_Mbz);
  virtual PetscErrorCode getVertLevsForTempAge(const PetscInt ta_Mz, const PetscInt ta_Mbz,
                                       PetscScalar &ta_dzEQ, PetscScalar &ta_dzbEQ, 
                                       PetscScalar *ta_zlevEQ, PetscScalar *ta_zblevEQ);

  // see iMutil.cc
  virtual PetscErrorCode getMagnitudeOf2dVectorField(IceModelVec2 vfx, IceModelVec2 vfy,
						     IceModelVec2 vmag);
  virtual int endOfTimeStepHook();
  virtual PetscErrorCode afterInitHook();
  virtual PetscErrorCode stampHistoryCommand();
  virtual PetscErrorCode stampHistoryEnd();
  virtual PetscErrorCode stampHistory(const char*);
  virtual PetscErrorCode stampHistoryAdd(const char*);
  virtual PetscErrorCode thicknessTooLargeCheck();

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
  virtual PetscErrorCode update2DViewer(const char scName, IceModelVec2 l2, const PetscScalar scale);
  virtual PetscErrorCode updateSliceViewer(const char scName, IceModelVec3 imv3, const PetscScalar scale);
  virtual PetscErrorCode updateSurfaceValuesViewer(const char scName, IceModelVec3 imv3, const PetscScalar scale);
  virtual PetscErrorCode updateSpeed2DViewer(const char scName, IceModelVec2 lu, IceModelVec2 lv, 
                   const PetscScalar scale, const PetscTruth doLog, 
                   const PetscScalar log_missing);
  virtual PetscErrorCode updateSpeedSurfaceValuesViewer(const char scName, IceModelVec3 imv3_u, IceModelVec3 imv3_v, 
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
  Vec g2;    // Global work vector
  IceModelVec2 vWork2d[nWork2d];
  // 3D working space (with specific purposes)
  IceModelVec3 Tnew3, taunew3;
  IceModelVec3 Sigmastag3[2], Istag3[2];

protected:
  int have_ssa_velocities;	// use vubarSSA and vvbarSSA from a previous
				// run if 1, otherwise set them to zero in
				// IceModel::initSSA
  IceModelVec2 vubarSSA, vvbarSSA;

private:
#if (PISM_LOG_EVENTS)
  // for event logging; see run() and velocity()
  int siaEVENT, ssaEVENT, velmiscEVENT, 
      beddefEVENT, pddEVENT, massbalEVENT, tempEVENT;
#endif

  // Pieces of the SSA Velocity routine defined in iMssa.cc.
  // Note these do not initialize correctly for derived classes if made
  // "private" however, derived classes should not need access to the details
  // of the linear system which uses these
  KSP SSAKSP;
  Mat SSAStiffnessMatrix;
  Vec SSAX, SSARHS;  // Global vectors for solution of the linear system
  Vec SSAXLocal; // We need a local copy of the solution to map back to a DA based vector
  VecScatter SSAScatterGlobalToLocal;

protected:
  // This is related to the snapshot saving feature
  char snapshots_filename[PETSC_MAX_PATH_LEN];
  bool save_snapshots, file_is_ready;

  // equally spaced snapshots
  PetscTruth save_at_equal_intervals;
  double first_snapshot, snapshot_dt, next_snapshot, last_snapshot;

  // unequally spaced snapshots
  static const int max_n_snapshots = 20;
  double save_at[max_n_snapshots];
  int n_snapshots, current_snapshot;

  PetscErrorCode init_snapshots_from_options();
  PetscErrorCode write_snapshot();
};

#endif /* __iceModel_hh */


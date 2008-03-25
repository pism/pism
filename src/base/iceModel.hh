// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
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
#include "grid.hh"
#include "materials.hh"
#include "forcing.hh"
#include "beddefLC.hh"
#include "iceModelVec.hh"

// following simply unclutters iceModel.hh:
#include "iceModelpreamble.hh" 

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
  PetscErrorCode setExecName(const char *my_executable_short_name);

  // see iMbootstrap.cc 
  PetscErrorCode bootstrapFromFile_netCDF(const char *fname);
  PetscErrorCode readShelfStreamBCFromFile_netCDF(const char *fname);

  // see iMoptions.cc
  virtual PetscErrorCode setFromOptions();

  // see iMtests.cc; FIXME: these should go in a derived class
  PetscErrorCode testIceModelVec3();
  PetscErrorCode testIceModelVec3Bedrock();

  // see iMutil.cc
  virtual PetscErrorCode initFromOptions();
  virtual PetscErrorCode initFromOptions(PetscTruth doHook);

  // see iMIO.cc
  PetscErrorCode initFromFile_netCDF(const char *);
  PetscErrorCode writeFiles(const char* defaultbasename);
  PetscErrorCode writeFiles(const char* defaultbasename, 
                            const PetscTruth forceFullDiagnostics);

protected:
   static const int MASK_SHEET;
   static const int MASK_DRAGGING;
   static const int MASK_FLOATING;
   static const int MASK_FLOATING_OCEAN0;

  IceGrid        &grid;

  PolarStereoParams  psParams;
  
  NCTool         nct;

  IceType               *ice;
//  IceType               &ice;
  BasalType             *basal;
  BedrockThermalType    bed_thermal;
  DeformableEarthType   bed_deformable;
  SeaWaterType          ocean;
  FreshWaterType        porewater;

  // state variables
  Vec vh, vH, vbed,             // 2D vectors; Mx x My
    vAccum, vTs,                // accumulation, surface temp
    vAccumSnow,                 // see iMpdd.cc; vAccum is net (including ablation)
    vMask,                      // mask for flow type
    vGhf,                       // geothermal flux
    vHmelt, vbasalMeltRate,     // thickness and rate of production of basal
                                // meltwater (ice-equivalent)
    vdHdt, vuplift,                 
    vub, vvb,                   // basal vels on standard grid
    vRb,                        // basal frictional heating on regular grid
    vubar, vvbar,               // vertically-averaged vels (u bar and v bar) on
                                // standard grid
    vtauc, vtillphi, vbeta;     // basal fields; plastic/viscous coefficients
  Vec*         vuvbar;                  // 2D; vuvbar[0] and vuvbar[1] are 
                                        //   u bar and v bar on staggered grid,
  Vec          vLongitude, vLatitude;

  IceModelVec3        u3, v3, w3, Sigma3, T3, tau3;  // note gs3 was one of these
  IceModelVec3Bedrock Tb3;

  // parameters
  PetscInt    flowLawNumber;
  PetscScalar maxdt, muSliding, enhancementFactor, initial_age_years_default;
  PetscScalar dt, dtTempAge, dt_force; // current mass cont. and temp/age 
                                       //   time steps in seconds
  PetscScalar constantNuForSSA, constantHardnessForSSA, min_thickness_SSA,
              regularizingVelocitySchoof, regularizingLengthSchoof,
              ssaRelativeTolerance, ssaEpsilon, beta_default_drag_SSA;
  PetscInt    ssaMaxIterations;
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
  PetscScalar isothermalFlux_n_exponent, min_temperature_for_SIA_sliding, Hmelt_max;
  PetscInt    tempskipCountDown, tempskipMax, noSpokesLevel, maxLowTempCount;
  PetscScalar globalMinAllowedTemp, oceanHeatFlux;

  // flags
  PetscTruth  doMassConserve, doTemp, doBedDef, doBedIso, flowLawUsesGrainSize;
  PetscTruth  initialized_p, thermalBedrock, includeBMRinContinuity, updateHmelt,
              isDrySimulation, holdTillYieldStress, useConstantTillPhi;
  PetscTruth  useSSAVelocity, doPlasticTill, doPseudoPlasticTill,
              doSuperpose, pureSuperpose, useConstantNuForSSA, 
              useConstantHardnessForSSA, computeSurfGradInwardSSA, leaveNuAloneSSA;
  PetscTruth  yearsStartRunEndDetermined, doAdaptTimeStep, doOceanKill, allowAboveMelting;
  PetscTruth  realAgeForGrainSize;
  PetscTruth  showViewers, ssaSystemToASCIIMatlab, doTempSkip, reportHomolTemps;
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

  char history[HISTORY_STRING_LENGTH]; // history of commands used to generate this 
                                       // instance of IceModel
  char executable_short_name[PETSC_MAX_PATH_LEN];
  
protected:
  // see iceModel.cc
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode destroyVecs();
  void setMaxTimeStepYears(PetscScalar); // use this value for adaptive stepping
  void setAdaptTimeStepRatio(PetscScalar);
  PetscErrorCode setStartYear(PetscScalar);
  PetscErrorCode setEndYear(PetscScalar);
  void setInitialAgeYears(PetscScalar d);
  void setAllGMaxVelocities(PetscScalar);
  void setNoViewers();
  void setConstantNuForSSA(PetscScalar);
  void setIsothermalFlux(PetscTruth use, PetscScalar n, PetscScalar A);
  PetscTruth isInitialized() const;
  PetscErrorCode updateSurfaceElevationAndMask();
  PetscErrorCode massBalExplicitStep();
  int intMask(PetscScalar);
  int modMask(PetscScalar);

  // see iMbasal.cc (basal mechanical model: either linear viscosity or
  //                 perfectly plastic according to doPlasticTill) 
  PetscErrorCode initBasalTillModel();
  virtual PetscErrorCode updateYieldStressFromHmelt();
  virtual PetscScalar basalVelocity(const PetscScalar x, const PetscScalar y, 
       const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
       const PetscScalar mu);
  virtual PetscScalar basalDragx(PetscScalar **beta, PetscScalar **tauc,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;
  virtual PetscScalar basalDragy(PetscScalar **beta, PetscScalar **tauc,
                                 PetscScalar **u, PetscScalar **v,
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
  Vec            vHlast, vbedlast;  // used for simple pointwise isostasy and to
                                   // compute uplift
  Vec            Hp0, bedp0,                       // vecs on proc zero for
                 Hstartp0, bedstartp0, upliftp0;   // passing to bdLC
  PetscErrorCode bedDefSetup();
  PetscErrorCode bedDefCleanup();
  PetscErrorCode bedDefStepIfNeeded();
  PetscErrorCode bed_def_step_iso();

  // see iMbootstrap.cc 
  PetscErrorCode putTempAtDepth();
  PetscErrorCode bootstrapSetBedrockColumnTemp(const PetscInt i, const PetscInt j,
                            const PetscScalar Ttopbedrock, const PetscScalar geothermflux);
  PetscErrorCode setMaskSurfaceElevation_bootstrap();

  // see iMdefaults.cc
  virtual PetscErrorCode setDefaults();

  // see iMforcing.cc
  IceSheetForcing  *dTforcing, *dSLforcing;
  PetscScalar    TsOffset, bedSLOffset;
  PetscErrorCode initForcingFromOptions();
  PetscErrorCode updateForcing();
  PetscErrorCode forcingCleanup();

  // see iMgrainsize.cc
  PetscErrorCode computeGrainSize_PseudoAge(const PetscScalar H, const PetscInt Mz, 
                          PetscScalar *w, PetscScalar *age_wspace, PetscScalar **gs);
  PetscScalar    grainSizeVostok(PetscScalar age) const;

  // see iMIO.cc
  bool hasSuffix(const char* fname, const char* suffix) const;
  PetscErrorCode warnUserOptionsIgnored(const char *fname);
  PetscErrorCode setStartRunEndYearsFromOptions(const PetscTruth grid_p_year_VALID);
  PetscErrorCode dumpToFile_netCDF(const char *fname);
  PetscErrorCode dumpToFile_diagnostic_netCDF(const char *diag_fname);
  PetscErrorCode regrid_netCDF(const char *fname);

  // see iMmatlab.cc
  bool           matlabOutWanted(const char name);
  PetscErrorCode VecView_g2ToMatlab(PetscViewer v, 
                                    const char *varname, const char *shorttitle);
  PetscErrorCode write2DToMatlab(PetscViewer v, const char singleCharName, 
                                 Vec l2, const PetscScalar scale);
  PetscErrorCode writeSliceToMatlab(PetscViewer v, const char singleCharName, 
                                    IceModelVec3 imv3, const PetscScalar scale);
  PetscErrorCode writeSurfaceValuesToMatlab(PetscViewer v, const char singleCharName, 
                                            IceModelVec3 imv3, const PetscScalar scale);
  PetscErrorCode writeSpeed2DToMatlab(PetscViewer v, const char scName, 
                          Vec lu, Vec lv, const PetscScalar scale, 
                          const PetscTruth doLog, const PetscScalar log_missing);
  PetscErrorCode writeSpeedSurfaceValuesToMatlab(PetscViewer v, const char scName, 
                          IceModelVec3 imv3_u, IceModelVec3 imv3_v, const PetscScalar scale, 
                          const PetscTruth doLog, const PetscScalar log_missing);
  PetscErrorCode writeLog2DToMatlab(PetscViewer v, const char scName, 
                          Vec l, const PetscScalar scale, const PetscScalar thresh,
                          const PetscScalar log_missing);
  PetscErrorCode writeSoundingToMatlab(PetscViewer v, const char scName, IceModelVec3 imv3,
                          const PetscScalar scale, const PetscTruth doTandTb);
  PetscErrorCode writeMatlabVars(const char *fname);
  PetscErrorCode writeSSAsystemMatlab(Vec vNu[2]);

  // see iMnames.cc; note tn is statically-initialized in iMnames.cc
  int cIndex(const char singleCharName);

  // see iMoptions.cc
  PetscErrorCode determineSpacingTypeFromOptions();

  // see iMpdd.cc (positive degree day model)
  gsl_rng     *pddRandGen;
  char        monthlyTempsFile[PETSC_MAX_PATH_LEN];
  Vec*        vmonthlyTs;                  // will allocate 12 2D Vecs 
  PetscTruth  doPDD, doPDDTrueRand, pddStuffCreated, pddRandStuffCreated,
              pddHaveMonthlyTempData;
  PetscScalar pddStdDev, pddFactorSnow, pddFactorIce, pddRefreezeFrac, 
              pddSummerWarming, pddSummerPeakDay;
  PetscErrorCode initPDDFromOptions();
  PetscErrorCode readMonthlyTempDataPDD();
  PetscErrorCode updateSurfaceBalanceFromPDD();
  PetscErrorCode putBackSnowAccumPDD();
  PetscErrorCode PDDCleanup();
  double         CalovGreveIntegrand(const double Tac);
  virtual PetscScalar getTemperatureFromYearlyCycle(
                  const PetscScalar summer_warming, const PetscScalar Tma, 
                  const PetscScalar day, const PetscInt i, const PetscInt j) const;
  virtual PetscScalar getSummerWarming(
                  const PetscScalar elevation, const PetscScalar latitude, 
                  const PetscScalar Tma) const;
  virtual double getSurfaceBalanceFromSnowAndPDD(
                     const double snowrate, const double mydt, const double pdds);

  // see iMreport.cc
  // note setVerbosityLevel(), verbosityLevelFromOptions(), and verbPrintf()
  // are all in iMreport.cc, but they are not IceModel methods
  PetscErrorCode computeFlowUbarStats
                      (PetscScalar *gUbarmax, PetscScalar *gUbarSIAav,
                       PetscScalar *gUbarstreamav, PetscScalar *gUbarshelfav,
                       PetscScalar *gicegridfrac, PetscScalar *gSIAgridfrac,
                       PetscScalar *gstreamgridfrac, PetscScalar *gshelfgridfrac);
  PetscErrorCode volumeArea(PetscScalar& gvolume,PetscScalar& garea,
                            PetscScalar& gvolSIA, PetscScalar& gvolstream, 
                            PetscScalar& gvolshelf);
  virtual PetscErrorCode summary(bool tempAndAge, bool useHomoTemp);
  virtual PetscErrorCode summaryPrintLine(
              const PetscTruth printPrototype, const PetscTruth tempAndAge,
              const PetscScalar year, const PetscScalar dt, 
              const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
              const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0);

  // see iMsia.cc
  PetscErrorCode surfaceGradientSIA();
  PetscErrorCode velocitySIAStaggered();
  PetscErrorCode basalSIA();
  PetscErrorCode velocities2DSIAToRegular();
  PetscErrorCode SigmaSIAToRegular();
  PetscErrorCode horizontalVelocitySIARegular();

  // see iMssa.cc
  PetscErrorCode initSSA();
  PetscErrorCode velocitySSA(PetscInt *numiter);
  PetscErrorCode velocitySSA(Vec vNu[2], PetscInt *numiter);
  PetscErrorCode setupGeometryForSSA(const PetscScalar minH);
  PetscErrorCode cleanupGeometryAfterSSA(const PetscScalar minH);
  virtual PetscErrorCode computeEffectiveViscosity(Vec vNu[2], PetscReal epsilon);
  PetscErrorCode testConvergenceOfNu(Vec vNu[2], Vec vNuOld[2], PetscReal *, PetscReal *);
  PetscErrorCode assembleSSAMatrix(Vec vNu[2], Mat A);
  PetscErrorCode assembleSSARhs(bool surfGradInward, Vec rhs);
  PetscErrorCode moveVelocityToDAVectors(Vec x);
  PetscErrorCode broadcastSSAVelocity(bool updateVelocityAtDepth);
  PetscErrorCode correctSigma();
  PetscErrorCode correctBasalFrictionalHeating();

  // see iMtemp.cc
  PetscErrorCode temperatureAgeStep();
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount);
  PetscErrorCode ageStep(PetscScalar* CFLviol);
  PetscErrorCode solveTridiagonalSystem(
           const PetscScalar* L, const PetscScalar* D, const PetscScalar* U,
           PetscScalar* x, const PetscScalar* rhs, PetscScalar* work, const int n) const;
  bool checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                      PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE);
  PetscErrorCode excessToFromBasalMeltLayer(
                      const PetscScalar rho_c, const PetscScalar z, const PetscScalar dz,
                      PetscScalar *Texcess, PetscScalar *Hmelt);
  PetscErrorCode getMzMbzForTempAge(PetscInt &ta_Mz, PetscInt &ta_Mbz);
  PetscErrorCode getVertLevsForTempAge(const PetscInt ta_Mz, const PetscInt ta_Mbz,
                                       PetscScalar &ta_dzEQ, PetscScalar &ta_dzbEQ, 
                                       PetscScalar *ta_zlevEQ, PetscScalar *ta_zblevEQ);

  // see iMutil.cc
  virtual int endOfTimeStepHook();
  virtual PetscErrorCode afterInitHook();
  PetscErrorCode stampHistoryCommand();
  PetscErrorCode stampHistoryEnd();
  PetscErrorCode stampHistory(const char*);
  PetscErrorCode stampHistoryAdd(const char*);
  PetscErrorCode computeMaxDiffusivity(bool updateDiffusViewer);
  PetscErrorCode computeBasalDrivingStress(Vec myVec);
  PetscErrorCode adaptTimeStepDiffusivity();
  virtual PetscErrorCode determineTimeStep(const bool doTemperatureCFL);

  // see iMvelocity.cc
  virtual PetscErrorCode velocity(bool updateSIAVelocityAtDepth);    
  PetscErrorCode vertVelocityFromIncompressibility();
  PetscScalar    capBasalMeltRate(const PetscScalar bMR);
  PetscErrorCode computeMax3DVelocities();
  PetscErrorCode computeMax2DSlidingSpeed();
    
  // see iMviewers.cc
  int isViewer(char name);
  PetscErrorCode updateSoundings();
  PetscErrorCode updateOneSounding(const char scName, Vec l, const PetscScalar scale);
  PetscErrorCode getViewerDims(const PetscInt target_size, const PetscScalar Lx, const PetscScalar Ly,
                               PetscInt *xdim, PetscInt *ydim);
  PetscErrorCode createOneViewerIfDesired(const char singleCharName);
  PetscErrorCode createOneViewerIfDesired(const char singleCharName, const char* title);
  PetscErrorCode createOneViewerIfDesired(PetscViewer* v, 
                                          const char singleCharName, const char* title);
  PetscErrorCode createViewers();
  PetscErrorCode update2DViewer(const char scName, Vec l2, const PetscScalar scale);
  PetscErrorCode updateSliceViewer(const char scName, IceModelVec3 imv3, const PetscScalar scale);
  PetscErrorCode updateSurfaceValuesViewer(const char scName, IceModelVec3 imv3, const PetscScalar scale);
  PetscErrorCode updateSpeed2DViewer(const char scName, Vec lu, Vec lv, 
                   const PetscScalar scale, const PetscTruth doLog, 
                   const PetscScalar log_missing);
  PetscErrorCode updateSpeedSurfaceValuesViewer(const char scName, IceModelVec3 imv3_u, IceModelVec3 imv3_v, 
                   const PetscScalar scale, const PetscTruth doLog, 
                   const PetscScalar log_missing);
  PetscErrorCode updateLog2DViewer(const char scName, Vec l,
                   const PetscScalar scale, const PetscScalar thresh, 
                   const PetscScalar log_missing);
  PetscErrorCode updateViewers();  // it calls updateSoundings()
  PetscErrorCode updateNuViewers(Vec vNu[2], Vec vNuOld[2], bool updateNu_tView);
  PetscErrorCode destroyViewers();

protected:
  // working space (a convenience)
  static const PetscInt nWork2d=6;
  Vec g2;    // Global work vector
  Vec* vWork2d;

private:
  // 3D working space (with specific purposes)
  IceModelVec3 Tnew3, taunew3;
  IceModelVec3 Sigmastag3[2], Istag3[2];
  
#if (LOG_PISM_EVENTS)
  // for event logging; see run() and velocity()
  int siaEVENT, ssaEVENT, velmiscEVENT, 
      beddefEVENT, pddEVENT, massbalEVENT, tempEVENT;
#endif

  // Pieces of the SSA Velocity routine defined in iMssa.cc.
  // Note these do not initialize correctly for derived classes if made
  // "private" however, derived classes should not need access to the details
  // of the linear system which uses these
  Vec vubarSSA, vvbarSSA;
  KSP SSAKSP;
  Mat SSAStiffnessMatrix;
  Vec SSAX, SSARHS;  // Global vectors for solution of the linear system
  Vec SSAXLocal; // We need a local copy of the solution to map back to a DA based vector
  VecScatter SSAScatterGlobalToLocal;
};

#endif /* __iceModel_hh */


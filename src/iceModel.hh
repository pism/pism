// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
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
#include "beddefLC.hh"

using namespace std;

// this structure is needed in regrid only (see iMutil.cc)
struct InterpCtx {
  Mat A;
  DA dac;
  DA daf;
  Vec gc;
  Vec gf;
  Vec coarse;
  Vec fine;
};

struct MaskInterp {
  int number_allowed;
  int allowed_levels[50];// must be strictly increasing
};

// this structure is to save polar stereographic projection
// defaults are for Antarctic ice sheet
struct PolarStereoParams {
  // these are "double" and not "float" ultimately because of how ncgen works
  double svlfp; // straight_vertical_longitude_from_pole; defaults to 0
  double lopo;  // latitude_of_projection_origin; defaults to 90
  double sp;    // standard_parallel; defaults to -71
};

// this is a utility to get a IceType based on options from user, for initialization
// of an instance of IceModel below; see iceModel.cc
PetscErrorCode getFlowLawFromUser(MPI_Comm com, IceType* &ice, PetscInt &flowLawNum);


// this utility prints only when verbosityLevel >= thresh; see iMutil.cc
extern PetscInt verbosityLevel;
PetscErrorCode setVerbosityLevel(PetscInt level);
PetscErrorCode verbosityLevelFromOptions();
PetscErrorCode verbPrintf(const int thresh, MPI_Comm comm,const char format[],...);


class IceModel {
public:
  // see iceModel.cc
  IceModel(IceGrid &g, IceType &i);
  virtual ~IceModel();
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode destroyVecs();
  void setDoTimeStep(PetscTruth);
  void setTimeStepYears(PetscScalar); // use a constant time step
  void setMaxTimeStepYears(PetscScalar); // use this value for adaptive stepping
  void setAdaptTimeStepRatio(PetscScalar);
  PetscErrorCode setStartYear(PetscScalar);
  PetscErrorCode setEndYear(PetscScalar);
  void setInitialAgeYears(PetscScalar d);
  void setShowViewers(PetscTruth);
  void setDoMassConserve(PetscTruth);
  void setDoTemp(PetscTruth);
  void setIncludeBMRinContinuity(PetscTruth);
  void setDoGrainSize(PetscTruth);
  void setDoBedDef(PetscTruth);
  void setDoBedIso(PetscTruth);
  void setAllGMaxVelocities(PetscScalar);
  void setThermalBedrock(PetscTruth);
  void setOceanKill(PetscTruth);
  void setNoViewers();
  void setUseMacayealVelocity(PetscTruth);
  void setDoSuperpose(PetscTruth);
  void setConstantNuForMacAyeal(PetscScalar);
  void setRegularizingVelocitySchoof(PetscScalar);
  void setRegularizingLengthSchoof(PetscScalar);
  void setMacayealEpsilon(PetscScalar);
  void setMacayealRelativeTolerance(PetscScalar);
  void setIsDrySimulation(PetscTruth);
  void setEnhancementFactor(PetscScalar);
  void setMuSliding(PetscScalar);
  void setGSIntervalYears(PetscScalar);
  void setBedDefIntervalYears(PetscScalar);
  void setNoSpokes(PetscInt);
  void setIsothermalFlux(PetscTruth use, PetscScalar n, PetscScalar A);
  void setIsothermalFlux(PetscTruth);
  PetscTruth isInitialized() const;

  // see iceModel.cc
  virtual PetscErrorCode run();
  virtual PetscErrorCode additionalAtStartTimestep();
  virtual PetscErrorCode additionalAtEndTimestep();
  virtual PetscErrorCode diagnosticRun();

  // see iMdefaults.cc
  virtual PetscErrorCode setDefaults();

  // see iMoptions.cc
  virtual PetscErrorCode setFromOptions();

  // see iMutil.cc
  virtual PetscErrorCode initFromOptions();
  virtual PetscErrorCode initFromOptions(PetscTruth doHook);

  // see iMviewers.cc
  PetscErrorCode setSoundingFromOptions();

  // see iMvelocity.cc
  PetscErrorCode velocity(bool updateSIAVelocityAtDepth);
    
  // see iMIO.cc
  PetscErrorCode initFromFile(const char *);
  virtual PetscErrorCode dumpToFile_Matlab(const char *);
  PetscErrorCode writeFiles(const char* basename);
  PetscErrorCode writeFiles(const char* basename, const PetscTruth fullDiagnostics);

  // see iMIOnetcdf.cc
  PetscErrorCode bootstrapFromFile_netCDF(const char *fname);
  PetscErrorCode initFromFile_netCDF(const char *fname);
  PetscErrorCode readShelfStreamBCFromFile_netCDF(const char *fname);
  PetscErrorCode dumpToFile_netCDF(const char *fname);
  PetscErrorCode dumpToFile_diagnostic_netCDF(const char *diag_fname);

protected:
   static const int MASK_SHEET;
   static const int MASK_DRAGGING;
   static const int MASK_FLOATING;
   static const int MASK_FLOATING_OCEAN0;

   static const PetscScalar DEFAULT_START_YEAR;
   static const PetscScalar DEFAULT_RUN_YEARS;

  //used in iMutil.cc
   static const PetscScalar DEFAULT_ADDED_TO_SLOPE_FOR_DIFF_IN_ADAPTIVE;
   static const PetscScalar DEFAULT_ADDED_TO_GDMAX_ADAPT;
   static const PetscScalar DEFAULT_ADAPT_TIMESTEP_RATIO;

  //used in iMIO.cc and iMIOnetcdf.cc
   static const PetscScalar DEFAULT_h_VALUE_MISSING;
   static const PetscScalar DEFAULT_H_VALUE_MISSING;
   static const PetscScalar DEFAULT_BED_VALUE_MISSING;
   static const PetscScalar DEFAULT_ACCUM_VALUE_MISSING;
   static const PetscScalar DEFAULT_SURF_TEMP_VALUE_MISSING;
   static const PetscScalar DEFAULT_GEOTHERMAL_FLUX_VALUE_MISSING;
   static const PetscScalar DEFAULT_ACCUMULATION_IN_OCEAN0;

  //used in iMvelocity.cc
   static const PetscScalar DEFAULT_MINH_MACAYEAL;  // m; minimum thickness for MacAyeal velocity computation
   static const PetscScalar DEFAULT_MIN_SHEET_TO_DRAGGING;   // m/a; critical SIA speed for switch SIA --> MacAyeal
   static const PetscScalar DEFAULT_MAX_SPEED_DRAGGING_TO_SHEET;  // m/a; crit Mac speed for switch MacAyeal --> SIA
   static const PetscScalar DEFAULT_MAX_SPEEDSIA_DRAGGING_TO_SHEET;    // m/a; crit SIA speed for switch MacAyeal --> SIA
   static const PetscScalar DEFAULT_MAX_VEL_FOR_CFL;  // max velocity used in CFL calculation if velocity is not actually calculated  (which would be weird)
   static const PetscScalar DEFAULT_MAXSLOPE_MACAYEAL; // no units/pure number; cap to avoid bad behavior
   static const PetscScalar DEFAULT_EPSILON_MACAYEAL;  // units?;  initial amount of (denominator) regularization in
  // computation of effective viscosity
   static const PetscScalar DEFAULT_EPSILON_MULTIPLIER_MACAYEAL;  // no units/pure number; epsilon goes up by this ratio when
  // previous value failed
   static const PetscScalar DEFAULT_VERT_VEL_MACAYEAL;  // temp evolution uses this value; incompressibility not satisfied
   static const PetscScalar DEFAULT_BASAL_DRAG_COEFF_MACAYEAL;
  static const PetscScalar DEFAULT_TAUC;
  //used in iMvelocity.C and iMutil.C
   static const PetscScalar DEFAULT_MIN_TEMP_FOR_SLIDING;  // note less than ice.meltingTemp;
  // if above this value then decide to slide
   static const PetscScalar DEFAULT_INITIAL_AGE_YEARS;  // age for ice to start age computation
   static const PetscScalar DEFAULT_GRAIN_SIZE;  // size of ice grains when assumed constant

  // used in iMtemp.cc
  static const PetscScalar DEFAULT_OCEAN_HEAT_FLUX;
  static const PetscScalar DEFAULT_MAX_HMELT;

  IceGrid        &grid;

  PolarStereoParams  psParams;

  IceType        &ice;
  BasalType      *basal;
  BedrockType    bedrock;
  SeaWaterType   ocean;
  FreshWaterType porewater;

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
    vtauc, vbeta;               // basal fields; plastic/viscous coefficients
  Vec*         vuvbar;                  // 2D; vuvbar[0] and vuvbar[1] are 
                                        //   u bar and v bar on staggered grid,
  Vec*         vDf;                     // vDf[0],vDf[1] are diffusivity on staggered grid
  Vec          vu, vv, vw,              // 3D: standard grid, Mx x My x Mz
               vSigma, vT, vgs, vtau;   //   strain-heating, temp, grain size, age
  Vec          vTb;                     // 3D bed: Mx x My x Mbz
  Vec          vLongitude, vLatitude;

  // parameters
  PetscScalar maxdt, muSliding, enhancementFactor;
  PetscScalar dt, dtTempAge;    // current mass cont. and temp/age time steps in seconds
  PetscScalar dt_force, maxdt_temporary; // might be set by additionalAt??Timestep()
  PetscScalar constantNuForMacAyeal, constantHardnessForMacAyeal,
              regularizingVelocitySchoof, regularizingLengthSchoof,
              macayealRelativeTolerance, macayealEpsilon;
  PetscInt    macayealMaxIterations;
  PetscScalar startYear, endYear;
  PetscScalar gsIntervalYears, bedDefIntervalYears, adaptTimeStepRatio;
  PetscScalar CFLviolcount;    // really is just a count, but PetscGlobalSum requires this type
  PetscScalar dt_from_diffus, dt_from_cfl, CFLmaxdt, CFLmaxdt2D, gDmax;
  PetscScalar gmaxu, gmaxv, gmaxw;  // global maximums on 3D grid for abs value 
                                    // of 3D components of velocities
  PetscScalar gdHdtav, dvoldt; // average value in map-plane (2D) of dH/dt (where there is ice) 
                               //   [units m/s] and d(volume)/dt [units m^3/s]
  PetscScalar isothermalFlux_n_exponent, isothermalFlux_A_softness;
  PetscInt    tempskipCountDown, tempskipMax, noSpokesLevel;

  // flags
  PetscTruth  doMassConserve, doTemp, doGrainSize, doBedDef, doBedIso;
  PetscTruth  initialized_p, thermalBedrock, includeBMRinContinuity, isDrySimulation;
  PetscTruth  useMacayealVelocity, doPlasticTill, doSuperpose, useConstantNuForMacAyeal, 
              useConstantHardnessForMacAyeal, computeSurfGradInwardMacAyeal;
  PetscTruth  yearsStartRunEndDetermined, doAdaptTimeStep, doOceanKill, allowAboveMelting;
  PetscTruth  showViewers, doTempSkip;
  PetscTruth  createVecs_done, createViewers_done, createBasal_done;
  PetscTruth  computeSIAVelocities, useIsothermalFlux;
  char        adaptReasonFlag;

  // viewer and sounding
  char         diagnostic[PETSC_MAX_PATH_LEN], diagnosticBIG[PETSC_MAX_PATH_LEN];
  PetscViewer  uvbarView[2], nuView[2], lognuView, NuView[2];
  PetscViewer  HView, hView, accumView, bedView, HmeltView, basalmeltView, maskView;
  PetscViewer  speedView, ubarView, vbarView, ghfView, upliftView, TsView;
  PetscViewer  T2View, TView, uView, vView, wView, SigmaView, SigmaMapView;
  PetscViewer  slidespeedView, RbView, gsView, gsMapView, betaView, taucView;
  PetscViewer  dhView, diffusView, tauView, tauMapView, umapView, vmapView, wmapView;
  PetscViewer  surfHorSpeedView, surfuView, surfvView, surfwView;
  PetscDrawLG  kspLG;
  PetscInt     id, jd, kd;
  Vec          Td, wd, ud, vd, Sigmad, gsd, taud;

  // working space (a convenience)
  Vec g2, g3, g3b;    // Global work vectors
  Vec* vWork3d;
  Vec* vWork2d;
  static const PetscInt nWork3d=6, nWork2d=10;

protected:
  // see iceModel.cc
  PetscErrorCode updateSurfaceElevationAndMask();
  PetscErrorCode massBalExplicitStep();
  int intMask(PetscScalar);
  int modMask(PetscScalar);

  // see iMutil.cc
  virtual int endOfTimeStepHook();
  virtual PetscErrorCode afterInitHook();
  PetscErrorCode stampHistoryCommand();
  PetscErrorCode stampHistoryEnd();
  PetscErrorCode stampHistory(const char*);
  PetscErrorCode stampHistoryString(const char*);
  PetscErrorCode computeFlowUbarStats
                      (PetscScalar *gUbarmax, PetscScalar *gUbarSIAav,
                       PetscScalar *gUbarstreamav, PetscScalar *gUbarshelfav,
                       PetscScalar *gicegridfrac, PetscScalar *gSIAgridfrac,
                       PetscScalar *gstreamgridfrac, PetscScalar *gshelfgridfrac);
  PetscErrorCode computeMaxDiffusivity(bool updateDiffusViewer);
  PetscErrorCode adaptTimeStepDiffusivity();
  virtual PetscErrorCode determineTimeStep(const bool doTemperatureCFL);
  PetscErrorCode volumeArea(PetscScalar& gvolume,PetscScalar& garea,
                            PetscScalar& gvolSIA, PetscScalar& gvolstream, 
                            PetscScalar& gvolshelf);
  virtual PetscErrorCode summary(bool,bool);
  virtual PetscErrorCode summaryPrintLine(
              const PetscTruth printPrototype, const PetscTruth tempAndAge,
              const PetscScalar year, const PetscScalar dt, 
              const PetscInt tempskipCount, const char adaptReason,
              const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
              const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0);
  PetscErrorCode getHorSliceOf3D(Vec v3D, Vec &gslice, PetscInt k);
  PetscErrorCode getSurfaceValuesOf3D(Vec v3D, Vec &g2D);

  // see iMregrid.cc (these use nc_util.cc heavily)
  PetscErrorCode regrid(const char *regridFile);
  PetscErrorCode regrid_netCDF(const char *fname);

  // see iMgrainsize.cc
  PetscErrorCode updateGrainSizeIfNeeded();
  PetscErrorCode updateGrainSizeNow();
  void setConstantGrainSize(PetscScalar d);
  PetscScalar grainSizeVostok(PetscScalar age) const;

  // see iMviewers.cc
  int isViewer(char name);
  PetscErrorCode initSounding();
  PetscErrorCode updateSoundings();
  PetscErrorCode updateViewers();  // it calls updateSoundings()
  PetscErrorCode createOneViewerIfDesired(PetscViewer *viewer, 
                                          char name, const char* title);
  PetscErrorCode createViewers();
  PetscErrorCode destroyViewers();

  // see iMbasal.cc (basal mechanical model: either linear viscosity or
  //                 perfectly plastic according to doPlasticTill) 
  PetscErrorCode initBasalTillModel();
  PetscErrorCode updateYieldStressFromHmelt();
  virtual PetscScalar basalVelocity(const PetscScalar x, const PetscScalar y, 
       const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
       const PetscScalar mu);
  virtual PetscScalar basalDragx(PetscScalar **beta, PetscScalar **tauc,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;
  virtual PetscScalar basalDragy(PetscScalar **beta, PetscScalar **tauc,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;

  // see iMpdd.cc (positive degree day model for ablation)
  gsl_rng     *pddRandGen;
  PetscTruth  doPDD, pddStuffCreated;
  PetscScalar pddStdDev, pddFactorSnow, pddFactorIce, pddRefreezeFrac, pddSummerWarming;
  static const PetscScalar DEFAULT_PDD_STD_DEV;
  static const PetscScalar DEFAULT_PDD_FACTOR_SNOW;
  static const PetscScalar DEFAULT_PDD_FACTOR_ICE;
  static const PetscScalar DEFAULT_PDD_REFREEZE_FRAC;
  static const PetscScalar DEFAULT_PDD_SUMMER_WARMING;
  PetscErrorCode initPDDFromOptions();
  bool IsIntegralYearPDD();
  virtual PetscScalar getTemperatureFromYearlyCycle(
                  const PetscScalar summer_warming, const PetscScalar Ta, const int day) const;
  virtual PetscScalar getSummerWarming(
                  const PetscScalar elevation, const PetscScalar latitude, const PetscScalar Ta) const;
  PetscErrorCode updateNetAccumFromPDD();
  PetscErrorCode putBackSnowAccumPDD();
  PetscErrorCode PDDCleanup();
  
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

  // see iMmacayeal.cc
  PetscErrorCode velocityMacayeal();

  PetscErrorCode setupForMacayeal(const PetscScalar minH);
  PetscErrorCode cleanupAfterMacayeal(const PetscScalar minH);
  virtual PetscErrorCode computeEffectiveViscosity(Vec vNu[2], PetscReal epsilon);
  PetscErrorCode testConvergenceOfNu(Vec vNu[2], Vec vNuOld[2], PetscReal *, PetscReal *);
  PetscErrorCode assembleMacayealMatrix(Vec vNu[2], Mat A);
  PetscErrorCode assembleMacayealRhs(bool surfGradInward, Vec rhs);
  PetscErrorCode moveVelocityToDAVectors(Vec x);
  PetscErrorCode broadcastMacayealVelocity();
  PetscErrorCode correctSigma();
  PetscErrorCode correctBasalFrictionalHeating();
  PetscErrorCode updateNuViewers(Vec vNu[2], Vec vNuOld[2], bool updateNu_tView);

  // see iMtemp.cc
  PetscErrorCode temperatureAgeStep();
  PetscErrorCode temperatureStep();
  PetscErrorCode ageStep(PetscScalar* CFLviol);
  PetscErrorCode solveTridiagonalSystem(
           const PetscScalar* L, const PetscScalar* D, const PetscScalar* U,
           PetscScalar* x, const PetscScalar* rhs, PetscScalar* work, const int n) const;
  bool checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                      PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE);
  PetscErrorCode excessToFromBasalMeltLayer(PetscScalar rho_c, PetscScalar z,
                                            PetscScalar *Texcess, PetscScalar *Hmelt);

  // see iMvelocity.cc
  PetscErrorCode velocitySIAStaggered(bool faststep);
  PetscErrorCode basalSIAConditionsToRegular();
  PetscErrorCode SigmaSIAToRegular();
  PetscErrorCode horizontalVelocitySIARegular();
  PetscErrorCode verticalVelocitySIARegular();
  PetscErrorCode smoothSigma();
  PetscErrorCode vertAveragedVelocityToRegular();
  PetscErrorCode computeMax3DVelocities();
  PetscScalar    capBasalMeltRate(const PetscScalar bMR);
  
  // see iMIO.cc
  bool hasSuffix(const char* fname, const char* suffix) const;
  PetscErrorCode LVecView(DA da, Vec l, Vec g, PetscViewer v);
  PetscErrorCode LVecLoad(DA da, Vec l, Vec g, PetscViewer v);
  PetscErrorCode setStartRunEndYearsFromOptions(const PetscTruth grid_p_year_VALID);

  // see iMIOnetcdf.cc
  PetscErrorCode putTempAtDepth();
  PetscErrorCode getIndZero(DA da, Vec vind, Vec vindzero, VecScatter ctx);
  PetscErrorCode nc_check(int stat);
  PetscErrorCode ncVarToDAVec(int ncid, int vid, DA da, Vec vecl,
                              Vec vecg, Vec vindzero);
  PetscErrorCode ncVarToDAVec(int ncid, int vid, DA da, Vec vecl,
                              Vec vecg, Vec vindzero, MaskInterp masktool);
  PetscErrorCode getFirstLast(int ncid, int vid, PetscScalar *gfirst, PetscScalar *glast);
  PetscErrorCode setMaskSurfaceElevation_bootstrap();
  PetscErrorCode setAccumInOcean();

  // see iMIOlegacy.cc
  // these are used in legacy pant and pismr runs and apply to "init.nc" which has non-standard
  // bootstrap conventions
  Vec    vbalvel;
  PetscErrorCode cleanJustNan_legacy();
  PetscErrorCode createMask_legacy(PetscTruth balVelRule);
  PetscErrorCode cleanInputData_legacy();
  PetscErrorCode reconfigure_legacy_Mbz();
  PetscErrorCode bootstrapFromFile_netCDF_legacyAnt(const char *fname);

private:
  // Pieces of the Macayeal Velocity routine defined in iMmacayeal.cc.
  // Note these do not initialize correctly for derived classes if made
  // "private" however, derived classes should not need access to the details
  // of the linear system which uses these
  KSP MacayealKSP;
  Mat MacayealStiffnessMatrix;
  Vec MacayealX, MacayealRHS;  // Global vectors for solution of the linear system
  Vec MacayealXLocal; // We need a local copy of the solution to map back to a DA based vector
  VecScatter MacayealScatterGlobalToLocal;

};

#endif /* __iceModel_hh */

// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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
#include "beddefLC.hh"

// following simply unclutters iceModel.hh:
#include "iceModelpreamble.hh" 

//! The base class for PISM.  Contains all essential variables, parameters, and flags for modelling an ice sheet.
/*!
@cond VERTCHANGE
  FIXME: add text about change of vertical variable
@endcond
 */
class IceModel {
public:
  // see iceModel.cc for implementation of constructor and destructor:
  IceModel(IceGrid &g, IceType &i);
  virtual ~IceModel(); // must be virtual merely because some members are virtual

  void setDoTimeStep(PetscTruth);
  void setTimeStepYears(PetscScalar); // use a constant time step
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

  // see iceModel.cc
  virtual PetscErrorCode run();
  virtual PetscErrorCode diagnosticRun();
  PetscScalar maxdt_temporary; // might be set by additionalAt??Timestep()
  virtual PetscErrorCode additionalAtStartTimestep();
  virtual PetscErrorCode additionalAtEndTimestep();

  // see iMdefaults.cc
  virtual PetscErrorCode setDefaults();

  // see iMoptions.cc
  virtual PetscErrorCode setFromOptions();

  // see iMutil.cc
  virtual PetscErrorCode initFromOptions();
  virtual PetscErrorCode initFromOptions(PetscTruth doHook);

  // see iMvelocity.cc
  virtual PetscErrorCode velocity(bool updateSIAVelocityAtDepth);
    
  // see iMIO.cc
  PetscErrorCode initFromFile(const char *);
  PetscErrorCode writeFiles(const char* defaultbasename);
  PetscErrorCode writeFiles(const char* defaultbasename, const PetscTruth forceFullDiagnostics);

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

  //used in iMvelocity.cc, iMsia.cc, iMssa.cc
   static const PetscScalar DEFAULT_MINH_SSA;  // m; minimum thickness for SSA velocity computation
   static const PetscScalar DEFAULT_MIN_SHEET_TO_DRAGGING;   // m/a; critical SIA speed for switch SIA --> SSA
   static const PetscScalar DEFAULT_MAX_SPEED_DRAGGING_TO_SHEET;  // m/a; crit SSA speed for switch SSA --> SIA
   static const PetscScalar DEFAULT_MAX_SPEEDSIA_DRAGGING_TO_SHEET;    // m/a; crit SIA speed for switch SSA --> SIA
   static const PetscScalar DEFAULT_MAX_VEL_FOR_CFL;  // max velocity used in CFL calculation if velocity is not actually calculated  (which would be weird)
   static const PetscScalar DEFAULT_MAXSLOPE_SSA; // no units/pure number; cap to avoid bad behavior
   static const PetscScalar DEFAULT_EPSILON_SSA;  // units?;  initial amount of (denominator) regularization in
  // computation of effective viscosity
   static const PetscScalar DEFAULT_EPSILON_MULTIPLIER_SSA;  // no units/pure number; epsilon goes up by this ratio when
  // previous value failed
   static const PetscScalar DEFAULT_VERT_VEL_SSA;  // temp evolution uses this value; incompressibility not satisfied
   static const PetscScalar DEFAULT_BASAL_DRAG_COEFF_SSA;
   static const PetscScalar DEFAULT_TAUC;
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
    vtauc, vtillphi, vbeta;     // basal fields; plastic/viscous coefficients
  Vec*         vuvbar;                  // 2D; vuvbar[0] and vuvbar[1] are 
                                        //   u bar and v bar on staggered grid,
  Vec          vu, vv, vw,              // 3D: standard grid, Mx x My x Mz
               vSigma, vT, vgs, vtau;   //   strain-heating, temp, grain size, age
  Vec          vTb;                     // 3D bed: Mx x My x Mbz
  Vec          vLongitude, vLatitude;

  // parameters
  PetscScalar maxdt, muSliding, enhancementFactor;
  PetscScalar dt, dtTempAge;    // current mass cont. and temp/age time steps in seconds
  PetscScalar dt_force;
  PetscScalar constantNuForSSA, constantHardnessForSSA,
              regularizingVelocitySchoof, regularizingLengthSchoof,
              ssaRelativeTolerance, ssaEpsilon;
  PetscInt    ssaMaxIterations;
  PetscScalar plastic_till_c_0, plastic_till_mu, plastic_till_pw_fraction, plasticRegularization;
  PetscScalar startYear, endYear;
  PetscScalar ssaIntervalYears, gsIntervalYears, bedDefIntervalYears, adaptTimeStepRatio;
  PetscScalar CFLviolcount;    // really is just a count, but PetscGlobalSum requires this type
  PetscScalar dt_from_diffus, dt_from_cfl, CFLmaxdt, CFLmaxdt2D, gDmax;
  PetscScalar gmaxu, gmaxv, gmaxw;  // global maximums on 3D grid for abs value 
                                    // of 3D components of velocities
  PetscScalar gdHdtav, dvoldt; // average value in map-plane (2D) of dH/dt (where there is ice) 
                               //   [units m/s] and d(volume)/dt [units m^3/s]
  PetscScalar isothermalFlux_n_exponent;
  PetscInt    tempskipCountDown, tempskipMax, noSpokesLevel, maxLowTempCount;
  PetscScalar globalMinAllowedTemp;

  // flags
  PetscTruth  doMassConserve, doTemp, doGrainSize, doBedDef, doBedIso;
  PetscTruth  initialized_p, thermalBedrock, includeBMRinContinuity, updateHmelt,
              isDrySimulation, holdTillYieldStress, useConstantTillPhi;
  PetscTruth  useSSAVelocity, doPlasticTill, doSuperpose, pureSuperpose, useConstantNuForSSA, 
              useConstantHardnessForSSA, computeSurfGradInwardSSA, leaveNuAloneSSA;
  PetscTruth  yearsStartRunEndDetermined, doAdaptTimeStep, doOceanKill, allowAboveMelting;
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
  Vec          Td, wd, ud, vd, Sigmad, gsd, taud;

protected:
  // see iceModel.cc
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode destroyVecs();
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

  // see iMgrainsize.cc
  PetscErrorCode updateGrainSizeIfNeeded();
  PetscErrorCode updateGrainSizeNow();
  void setConstantGrainSize(PetscScalar d);
  PetscScalar grainSizeVostok(PetscScalar age) const;

  // see iMIO.cc
  bool hasSuffix(const char* fname, const char* suffix) const;
  PetscErrorCode setStartRunEndYearsFromOptions(const PetscTruth grid_p_year_VALID);

  // see iMIOnetcdf.cc
  PetscErrorCode putTempAtDepth();
  PetscErrorCode getIndZero(DA da, Vec vind, Vec vindzero, VecScatter ctx);
  PetscErrorCode ncVarToDAVec(int ncid, int vid, DA da, Vec vecl,
                              Vec vecg, Vec vindzero);
  PetscErrorCode ncVarToDAVec(int ncid, int vid, DA da, Vec vecl,
                              Vec vecg, Vec vindzero, MaskInterp masktool);
  PetscErrorCode getFirstLast(int ncid, int vid, PetscScalar *gfirst, PetscScalar *glast);
  PetscErrorCode setMaskSurfaceElevation_bootstrap();
  PetscErrorCode regrid_netCDF(const char *fname);

  // see iMIOlegacy.cc
  // these are used in legacy pant and pismr runs and apply to "init.nc" which has non-standard
  // bootstrap conventions
  Vec    vbalvel;
  PetscErrorCode cleanJustNan_legacy();
  PetscErrorCode createMask_legacy(PetscTruth balVelRule);
  PetscErrorCode cleanInputData_legacy();
  PetscErrorCode reconfigure_legacy_Mbz();
  PetscErrorCode bootstrapFromFile_netCDF_legacyAnt(const char *fname);

  // see iMmatlab.cc
  bool           matlabOutWanted(const char name);
  PetscErrorCode VecView_g2ToMatlab(PetscViewer v, 
                                    const char *varname, const char *shorttitle);
  PetscErrorCode write2DToMatlab(PetscViewer v, const char singleCharName, 
                                 Vec l2, const PetscScalar scale);
  PetscErrorCode writeSliceToMatlab(PetscViewer v, const char singleCharName, 
                                    Vec l3, const PetscScalar scale);
  PetscErrorCode writeSurfaceValuesToMatlab(PetscViewer v, const char singleCharName, 
                                            Vec l3, const PetscScalar scale);
  PetscErrorCode writeSpeed2DToMatlab(PetscViewer v, const char scName, 
                          Vec lu, Vec lv, const PetscScalar scale, 
                          const PetscTruth doLog, const PetscScalar log_missing);
  PetscErrorCode writeSpeedSurfaceValuesToMatlab(PetscViewer v, const char scName, 
                          Vec lu, Vec lv, const PetscScalar scale, 
                          const PetscTruth doLog, const PetscScalar log_missing);
  PetscErrorCode writeLog2DToMatlab(PetscViewer v, const char scName, 
                          Vec l, const PetscScalar scale, const PetscScalar thresh,
                          const PetscScalar log_missing);
  PetscErrorCode writeSoundingToMatlab(PetscViewer v, const char scName, Vec l,
                          const PetscScalar scale, const PetscTruth doTandTb);
  PetscErrorCode writeMatlabVars(const char *fname);
  PetscErrorCode writeSSAsystemMatlab(Vec vNu[2]);

  // see iMnames.cc; note tn is statically-initialized in iMnames.cc
  int cIndex(const char singleCharName);

  // see iMpdd.cc (positive degree day model for ablation)
  gsl_rng     *pddRandGen;
  PetscTruth  doPDD, doPDDTrueRand, pddStuffCreated, pddRandStuffCreated;
  PetscScalar pddStdDev, pddFactorSnow, pddFactorIce, pddRefreezeFrac, 
              pddSummerWarming, pddSummerPeakDay;
  static const PetscScalar DEFAULT_PDD_STD_DEV;
  static const PetscScalar DEFAULT_PDD_FACTOR_SNOW;
  static const PetscScalar DEFAULT_PDD_FACTOR_ICE;
  static const PetscScalar DEFAULT_PDD_REFREEZE_FRAC;
  static const PetscScalar DEFAULT_PDD_SUMMER_WARMING;
  static const PetscScalar DEFAULT_PDD_SUMMER_PEAK_DAY;
  PetscErrorCode initPDDFromOptions();
  PetscErrorCode updateSurfaceBalanceFromPDD();
  PetscErrorCode putBackSnowAccumPDD();
  PetscErrorCode PDDCleanup();
  double         CalovGreveIntegrand(const double Tac);
  virtual PetscScalar getTemperatureFromYearlyCycle(
                  const PetscScalar summer_warming, const PetscScalar Tma, const PetscScalar day) const;
  virtual PetscScalar getSummerWarming(
                  const PetscScalar elevation, const PetscScalar latitude, const PetscScalar Tma) const;
  virtual double getSurfaceBalanceFromSnowAndPDD(
                     const double snowrate, const double mydt, const double pdds);

  // see iMreport.cc
  // note setVerbosityLevel(), verbosityLevelFromOptions(), and verbPrintf()
  // are all in iMreport.cc
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
  virtual PetscErrorCode temperatureStep();
  PetscErrorCode ageStep(PetscScalar* CFLviol);
  PetscErrorCode solveTridiagonalSystem(
           const PetscScalar* L, const PetscScalar* D, const PetscScalar* U,
           PetscScalar* x, const PetscScalar* rhs, PetscScalar* work, const int n) const;
  bool checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                      PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE);
  PetscErrorCode excessToFromBasalMeltLayer(PetscScalar rho_c, PetscScalar z,
                                            PetscScalar *Texcess, PetscScalar *Hmelt);

  // see iMutil.cc
  virtual int endOfTimeStepHook();
  virtual PetscErrorCode afterInitHook();
  PetscErrorCode stampHistoryCommand();
  PetscErrorCode stampHistoryEnd();
  PetscErrorCode stampHistory(const char*);
  PetscErrorCode stampHistoryString(const char*);
  PetscErrorCode computeMaxDiffusivity(bool updateDiffusViewer);
  PetscErrorCode computeBasalDrivingStress(Vec myVec);
  PetscErrorCode adaptTimeStepDiffusivity();
  virtual PetscErrorCode determineTimeStep(const bool doTemperatureCFL);
  PetscErrorCode getHorSliceOf3D(Vec v3D, Vec &gslice, PetscInt k);
  PetscErrorCode getSurfaceValuesOf3D(Vec v3D, Vec &g2D);

  // see iMvelocity.cc
  PetscErrorCode vertVelocityFromIncompressibility();
  PetscScalar    capBasalMeltRate(const PetscScalar bMR);
  PetscErrorCode smoothSigma();
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
  PetscErrorCode updateSliceViewer(const char scName, Vec l3, const PetscScalar scale);
  PetscErrorCode updateSurfaceValuesViewer(const char scName, Vec l3, const PetscScalar scale);
  PetscErrorCode updateSpeed2DViewer(const char scName, Vec lu, Vec lv, 
                   const PetscScalar scale, const PetscTruth doLog, 
                   const PetscScalar log_missing);
  PetscErrorCode updateSpeedSurfaceValuesViewer(const char scName, Vec lu, Vec lv, 
                   const PetscScalar scale, const PetscTruth doLog, 
                   const PetscScalar log_missing);
  PetscErrorCode updateLog2DViewer(const char scName, Vec l,
                   const PetscScalar scale, const PetscScalar thresh, 
                   const PetscScalar log_missing);
  PetscErrorCode updateViewers();  // it calls updateSoundings()
  PetscErrorCode updateNuViewers(Vec vNu[2], Vec vNuOld[2], bool updateNu_tView);
  PetscErrorCode destroyViewers();

private:
  // working space (a convenience)
  Vec g2, g3, g3b;    // Global work vectors
  Vec* vWork3d;
  Vec* vWork2d;
  static const PetscInt nWork3d=4, nWork2d=6;

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

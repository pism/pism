// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
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

#include <petscda.h>
#include <petscksp.h>

#if (WITH_NETCDF)
#include <netcdf.h>
#endif

#if (WITH_FFTW)
#include <fftw3.h>
#endif

#include "grid.hh"
#include "materials.hh"


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


// this is a utility to get a IceType based on options from user, for initialization
// of an instance of IceModel below
PetscErrorCode getFlowLawFromUser(MPI_Comm com, IceType* &ice, PetscInt &flowLawNum);


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
  PetscErrorCode setRunYears(PetscScalar);
  void setInitialAgeYears(PetscScalar d);
  void setShowViewers(PetscTruth);
  void setDoMassBal(PetscTruth);
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
  void setConstantNuForMacAyeal(PetscScalar);
  void setMacayealEpsilon(PetscScalar);
  void setMacayealRelativeTolerance(PetscScalar);
  void setIsDrySimulation(PetscTruth);
  void setEnhancementFactor(PetscScalar);
  void setMuSliding(PetscScalar);
  void setTempskip(PetscInt);
  void setGSIntervalYears(PetscScalar);
  void setBedDefIntervalYears(PetscScalar);
  void setBeVerbose(PetscTruth);
  void setAllowRegridding(PetscTruth);
  void setNoSpokes(PetscInt);
  void setIsothermalFlux(PetscTruth use, PetscScalar n, PetscScalar A);
  void setIsothermalFlux(PetscTruth);
  PetscTruth isInitialized() const;

  virtual PetscErrorCode run();

  // see iMutil.cc for default empty versions
  virtual PetscErrorCode additionalAtStartTimestep();
  virtual PetscErrorCode additionalAtEndTimestep();

  // see iMdefaults.cc
  virtual PetscErrorCode setDefaults();

  // see iMoptions.cc
  virtual PetscErrorCode setFromOptions();

  // see iMviewers.cc
  PetscErrorCode setSoundingFromOptions();

  // see iMvelocity.cc
  PetscErrorCode velocity(bool updateSIAVelocityAtDepth);
    
  // see iMIO.cc
  PetscErrorCode initFromFile(const char *);
  PetscErrorCode dumpToFile(const char *);
  virtual PetscErrorCode dumpToFile_Matlab(const char *);
  PetscErrorCode writeFiles(const char* basename);
  PetscErrorCode writeFiles(const char* basename, const char* formats);
  PetscErrorCode initFromOptions();

  // see iMIOnetcdf.cc
#if (WITH_NETCDF)
  PetscErrorCode initFromFile_netCDF(const char *fname);
#endif

protected:
   static const int MASK_SHEET;
   static const int MASK_DRAGGING;
   static const int MASK_FLOATING;
   static const int MASK_FLOATING_OCEAN0;

  //used in iMutil.cc
   static const PetscScalar DEFAULT_ADDED_TO_SLOPE_FOR_DIFF_IN_ADAPTIVE;
   static const PetscScalar DEFAULT_ADDED_TO_GDMAX_ADAPT;
   static const PetscScalar DEFAULT_ADAPT_TIMESTEP_RATIO;

  //used in iMIO.cc
   static const PetscScalar DEFAULT_h_VALUE_MISSING;
   static const PetscScalar DEFAULT_H_VALUE_MISSING;
   static const PetscScalar DEFAULT_BED_VALUE_MISSING;
   static const PetscScalar DEFAULT_ACCUM_VALUE_MISSING;
   static const PetscScalar DEFAULT_SURF_TEMP_VALUE_MISSING;

  //used in iMvelocity.cc
   static const PetscScalar DEFAULT_MINH_MACAYEAL;  // m; minimum thickness for MacAyeal velocity computation
   static const PetscScalar DEFAULT_MIN_SHEET_TO_DRAGGING;   // m/a; critical SIA speed for switch SIA --> MacAyeal
   static const PetscScalar DEFAULT_MAX_SPEED_DRAGGING_TO_SHEET;  // m/a; crit Mac speed for switch MacAyeal --> SIA
   static const PetscScalar DEFAULT_MAX_SPEEDSIA_DRAGGING_TO_SHEET;    // m/a; crit SIA speed for switch MacAyeal --> SIA
   static const PetscScalar DEFAULT_MAX_VEL_FOR_CFL;  // max velocity used in CFL calculation if velocity is not actually calculated  (which would be weird)
   static const PetscScalar DEFAULT_MAXSLOPE_MACAYEAL; // no units/pure number; cap to avoid bad behavior
   static const PetscInt    DEFAULT_MAX_ITERATIONS_MACAYEAL;
   static const PetscScalar DEFAULT_EPSILON_MACAYEAL;  // units?;  initial amount of (denominator) regularization in
  // computation of effective viscosity
   static const PetscScalar DEFAULT_EPSILON_MULTIPLIER_MACAYEAL;  // no units/pure number; epsilon goes up by this ratio when
  // previous value failed
   static const PetscScalar DEFAULT_VERT_VEL_MACAYEAL;  // temp evolution uses this value; incompressibility not satisfied
   static const PetscScalar DEFAULT_BASAL_DRAG_COEFF_MACAYEAL;
  //used in iMvelocity.C and iMutil.C
   static const PetscScalar DEFAULT_MIN_TEMP_FOR_SLIDING;  // note less than ice.meltingTemp;
  // if above this value then decide to slide
   static const PetscScalar DEFAULT_INITIAL_AGE_YEARS;  // age for ice to start age computation
   static const PetscScalar DEFAULT_GRAIN_SIZE;  // size of ice grains when assumed constant

  // used in iMtemp.cc
  static const PetscScalar DEFAULT_OCEAN_HEAT_FLUX;
  static const PetscScalar DEFAULT_MAX_HMELT;

  IceGrid      &grid;
  IceType      &ice;
  BedrockType  bedrock;
  DumbOceanType ocean;

  Vec          vh, vH, vbed,            // 2D vectors; Mx x My
               vAccum, vTs,             // accumulation, surface temp
               vMask,                   // mask for flow type
               vGhf,                    // geothermal flux
               vHmelt, vbasalMeltRate,  // thickness and rate of production of 
                                        // basal melt water (ice-equivalent)
               vuplift,                 
               vub, vvb,                // basal vels on standard grid
               vRb,                     // basal frictional heating on regular grid
               vubar, vvbar;            // vertically-averaged vels 
                                        //   (u bar and v bar) on standard grid
  Vec*         vuvbar;                  // 2D; vuvbar[0] and vuvbar[1] are 
                                        //   u bar and v bar on staggered grid
  Vec          vu, vv, vw,              // 3D: standard grid, Mx x My x Mz
               vSigma, vT, vgs, vtau;   //   strain-heating, temp, grain size, age
  Vec          vTb;                     // 3D bed: Mx x My x Mbz

  // parameters and flags
  PetscScalar muSliding, enhancementFactor;
  PetscScalar dt, dtTempAge;    // current mass cont. and temp/age time steps in seconds
  PetscScalar maxdt;
  PetscScalar dt_force, maxdt_temporary; // might be set by additionalAt??Timestep()
  char        adaptReasonFlag;
  PetscTruth  thermalBedrock, includeBMRinContinuity, isDrySimulation;
  PetscTruth  useMacayealVelocity, initialized_p, useConstantNuForMacAyeal;
  PetscScalar constantNuForMacAyeal, macayealRelativeTolerance, macayealEpsilon;
  PetscScalar adaptTimeStepRatio;
  PetscScalar startYear, endYear;
  PetscTruth  relativeEndYear, doAdaptTimeStep, doOceanKill, allowAboveMelting;
  PetscTruth  doMassBal, doTemp, doGrainSize, doBedDef, doBedIso;
  PetscTruth  showViewers, allowRegridding, beVerbose;
  PetscTruth  createVecs_done, createViewers_done;
  PetscInt    tempskip, noSpokesLevel;
  PetscScalar gsIntervalYears, bedDefIntervalYears;
  PetscScalar CFLviolcount;    // really is just a count, but PetscGlobalSum requires this type
  PetscScalar CFLmaxdt, gDmax;
  PetscScalar gmaxu, gmaxv, gmaxw;  // global maximums on 3D grid for abs value 
                                    // of 3D components of velocities
  PetscScalar gdHdtav, dvoldt; // average value in map-plane (2D) of dH/dt (where there is ice) 
                               //   [units m/s] and d(volume)/dt [units m^3/s]
  PetscTruth  useIsothermalFlux;
  PetscScalar isothermalFlux_n_exponent;
  PetscScalar isothermalFlux_A_softness;

  // viewer and sounding
  char         diagnostic[PETSC_MAX_PATH_LEN], diagnosticBIG[PETSC_MAX_PATH_LEN];
  PetscViewer  uvbarView[2], nuView[2], NuView[2];
  PetscViewer  HView, hView, accumView, bedView, HmeltView, basalmeltView, maskView;
  PetscViewer  speedView, ubarView, vbarView, ghfView, upliftView, TsView;
  PetscViewer  T2View, TView, uView, vView, wView, SigmaView, SigmaMapView;
  PetscViewer  slidespeedView, RbView, gsView, gsMapView;
  PetscViewer  dhView, diffusView, tauView, tauMapView, umapView, vmapView, wmapView;
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
  PetscErrorCode vPetscPrintf(MPI_Comm comm,const char format[],...);  // beVerbose modulated printf
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
  PetscErrorCode adaptTimeStepCFL();
  virtual PetscErrorCode determineTimeStep(const bool doTemperatureCFL);
  PetscErrorCode volumeArea(PetscScalar& gvolume,PetscScalar& garea,
                            PetscScalar& gvolSIA, PetscScalar& gvolstream, 
                            PetscScalar& gvolshelf);
  PetscErrorCode summary(bool,bool);
  PetscErrorCode checkForSymmetry(Vec vec, PetscReal *normx, PetscReal *normy,
                                   PetscInt stagger);

  // see iMregrid.cc
  PetscErrorCode regrid(const char *regridFile);
  PetscErrorCode getInterpCtx(const DA dac, const DA daf,
                              const IceModel &cmodel, InterpCtx &interpCtx);
  PetscErrorCode destroyInterpCtx(InterpCtx &i);
  PetscErrorCode regridVar(const char *vars, char c, const InterpCtx &i,
                           const Vec src, Vec dest);

  // see iMgrainsize.cc
  PetscErrorCode updateGrainSizeIfNeeded();
  PetscErrorCode updateGrainSizeNow();
  void setConstantGrainSize(PetscScalar d);
  PetscScalar grainSizeVostok(PetscScalar age) const;

  // see iMviewers.cc
  int isViewer(char name);
  PetscErrorCode initSounding();
  PetscErrorCode updateViewers();
  PetscErrorCode createOneViewerIfDesired(PetscViewer *viewer, 
                                          char name, const char* title);
  PetscErrorCode createViewers();
  PetscErrorCode destroyViewers();

  // see iMbasal.cc
  virtual PetscScalar basal(const PetscScalar x, const PetscScalar y, 
       const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
       const PetscScalar mu);
  PetscScalar basalDrag(const PetscScalar u, const PetscScalar v) const;

  // This is an unfortunate kludge, but I won't rewrite the whole Macayeal code
  // just to implement a goofy verification case.
  virtual PetscScalar basalDragx(PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;
  virtual PetscScalar basalDragy(PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const;

  // see iMbeddef.cc: possibly useful general tool for putting Vecs on processor zero
  Vec         g2natural;
  VecScatter  top0ctx;
  PetscErrorCode createScatterToProcZero(Vec& samplep0);
  PetscErrorCode destroyScatterToProcZero();
  PetscErrorCode putOnProcZero(Vec& vlocal, Vec& onp0);
  PetscErrorCode getFromProcZero(Vec& onp0, Vec& vlocal);

  // see iMbeddef.cc
  PetscScalar   lastBedDefUpdateYear;
  Vec           vHlast, vbedlast;  // used for simple pointwise isostasy and to
                                   // compute uplift
  // used for Lingle&Clark model:
  Vec           Hstartp0, bedstartp0, platep0fat, vleft, vright, lrmEp0;  
#if (WITH_FFTW)
  fftw_complex  *bdin, *bdout;  // note these are 2D arrays but must be sequential
  fftw_plan     bdplanfor,bdplanback;
#endif
  PetscScalar    volumeChange(Vec Hp0, Vec Hstartp0);
  PetscErrorCode conv2_same(Vec vA, const PetscInt mA, const PetscInt nA, 
                            Vec vB, const PetscInt mB, const PetscInt nB,
                            Vec &vresult);
  PetscErrorCode bedDefSetup();
  PetscErrorCode bedDefCleanup();
  PetscErrorCode bedDefStepIfNeeded();
  PetscErrorCode bed_def_step_iso(PetscScalar dtBedDef);
  PetscErrorCode bed_uplift_init_lc();
  PetscErrorCode bed_def_step_lc(PetscScalar dtBedDef);

  // see iMmacayeal.cc
  PetscErrorCode velocityMacayeal();

  PetscErrorCode setupForMacayeal(const PetscScalar minH, const PetscTruth adjustMask);
  PetscErrorCode cleanupAfterMacayeal(const PetscScalar minH);
  virtual PetscErrorCode computeEffectiveViscosity(Vec vNu[2], PetscReal epsilon);
  PetscErrorCode testConvergenceOfNu(Vec vNu[2], Vec vNuOld[2], PetscReal *, PetscReal *);
  PetscErrorCode assembleMacayealMatrix(Vec vNu[2], Mat A, Vec rhs);
  PetscErrorCode moveVelocityToDAVectors(Vec x);
  PetscErrorCode broadcastMacayealVelocity();
  PetscErrorCode correctSigma();
  PetscErrorCode correctBasalFrictionalHeating();

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
  PetscErrorCode velocitySIAStaggered();
  PetscErrorCode basalSIAConditionsToRegular();
  PetscErrorCode SigmaSIAToRegular();
  PetscErrorCode horizontalVelocitySIARegular();
  PetscErrorCode verticalVelocitySIARegular();
  PetscErrorCode smoothSigma();
  PetscErrorCode vertAveragedVelocityToRegular();
  PetscErrorCode computeMaxVelocities();
  PetscScalar    capBasalMeltRate(const PetscScalar bMR);
  
  // see iMIO.cc
  bool hasSuffix(const char* fname, const char* suffix) const;
  PetscErrorCode LVecView(DA da, Vec l, Vec g, PetscViewer v);
  PetscErrorCode LVecLoad(DA da, Vec l, Vec g, PetscViewer v);

  // see iMIOnetcdf.cc
  Vec    vbalvel;
  PetscErrorCode createMask(PetscTruth balVelRule);
  PetscErrorCode putTempAtDepth();
  PetscErrorCode getIndZero(DA da, Vec vind, Vec vindzero, VecScatter ctx);
  PetscErrorCode cleanInputData();
#if (WITH_NETCDF)
//   PetscErrorCode ncVarToDAVec(const NcVar *v, DA da, Vec vecl,
//                               Vec vecg, Vec vindzero);
  PetscErrorCode ncVarToDAVec(int ncid, int vid, DA da, Vec vecl,
                              Vec vecg, Vec vindzero);
#endif

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

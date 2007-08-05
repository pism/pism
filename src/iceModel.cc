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

#include <cmath>
#include <cstring>
#include <petscda.h>

#include "iceModel.hh"
#include "pism_signal.h"

PetscInt verbosityLevel;

// following numerical values have some significance; see updateSurfaceElevationAndMask() below
const int IceModel::MASK_SHEET = 1;
const int IceModel::MASK_DRAGGING = 2;
const int IceModel::MASK_FLOATING = 3;
// (modMask(mask[i][j]) == MASK_FLOATING) is criteria for floating; ..._OCEAN0 only used if -ocean_kill 
const int IceModel::MASK_FLOATING_OCEAN0 = 7;

PetscErrorCode getFlowLawFromUser(MPI_Comm com, IceType* &ice, PetscInt &flowLawNum) {
    PetscErrorCode ierr;
    PetscTruth     flowlawSet = PETSC_FALSE, useGK = PETSC_FALSE;

    ierr = PetscOptionsGetInt(PETSC_NULL, "-law", &flowLawNum, &flowlawSet); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-gk", &useGK); CHKERRQ(ierr);
    if (useGK==PETSC_TRUE) {
      flowlawSet = PETSC_TRUE;
      flowLawNum = 4;
    }

    ierr = verbPrintf((flowlawSet == PETSC_TRUE) ? 3 : 5,com, 
        "  [using flow law %d (where 0=Paterson-Budd,1=cold P-B,2=warm P-B,3=Hooke,4=Goldsby-Kohlstedt)]\n",
        flowLawNum); CHKERRQ(ierr);
    
    switch (flowLawNum) {
      case 0: // Paterson-Budd
        ice = new ThermoGlenIce;  
        break;
      case 1: // cold part of P-B
        ice = new ThermoGlenArrIce;  
        break;
      case 2: // warm part of P-B
        ice = new ThermoGlenArrIceWarm;  
        break;
      case 3: // Hooke
        ice = new ThermoGlenIceHooke;
        break;
      case 4: // Goldsby Kohlstedt
        ice = new HybridIce;  
        break;
      case 5: // Goldsby Kohlstedt stripped down
        ice = new HybridIceStripped;  
        break;
      default:
        SETERRQ(1,"\nflow law number for to initialize IceModel must be 0,1,2,3,4,5\n");
    }
    return 0;
}


IceModel::IceModel(IceGrid &g, IceType &i): grid(g), ice(i) {
  PetscErrorCode ierr;

  pism_signal = 0;
  signal(SIGTERM, pism_signal_handler);
  signal(SIGUSR1, pism_signal_handler);

  createBasal_done = PETSC_FALSE;
  top0ctx_created = PETSC_FALSE;
  createVecs_done = PETSC_FALSE;
  createViewers_done = PETSC_FALSE;
  ierr = initIceParam(grid.com, &grid.p);
  if (ierr != 0) {
    verbPrintf(1,grid.com, "Error setting IceParams (in IceGrid).\n");
    PetscEnd();
  }        
  ierr = setDefaults();
  if (ierr != 0) {
    verbPrintf(1,grid.com, "Error setting defaults.\n");
    PetscEnd();
  }
  
  psParams.svlfp = 0.0;  // default polar stereographic projection settings
  psParams.lopo = 90.0;
  psParams.sp = -71.0;
}


IceModel::~IceModel() {
  if (createVecs_done == PETSC_TRUE) {
    destroyVecs();
  }
  if (createViewers_done == PETSC_TRUE) {
    destroyViewers();
  }
  if (createBasal_done == PETSC_TRUE) {
    delete basal;
  }
}


PetscErrorCode IceModel::createVecs() {
  PetscErrorCode ierr;

  if (createVecs_done == PETSC_TRUE) {
    ierr = destroyVecs(); CHKERRQ(ierr);
  }
  
  ierr = DACreateLocalVector(grid.da3, &vu); CHKERRQ(ierr);
  ierr = VecDuplicate(vu, &vv); CHKERRQ(ierr);
  ierr = VecDuplicate(vu, &vw); CHKERRQ(ierr);
  ierr = VecDuplicate(vu, &vSigma); CHKERRQ(ierr);
  ierr = VecDuplicate(vu, &vT); CHKERRQ(ierr);
  ierr = VecDuplicate(vu, &vtau); CHKERRQ(ierr);
  ierr = VecDuplicate(vu, &vgs); CHKERRQ(ierr);

  ierr = DACreateLocalVector(grid.da3b, &vTb); CHKERRQ(ierr);

  ierr = DACreateLocalVector(grid.da2, &vh); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vH); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vbed); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vAccum); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vTs); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vMask); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vGhf); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vubar); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vvbar); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vub); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vvb); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vRb); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vHmelt); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vbasalMeltRate); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vuplift); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vdHdt); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vbeta); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vtauc); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vLongitude); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &vLatitude); CHKERRQ(ierr);

  ierr = VecDuplicateVecs(vh, 2, &vDf); CHKERRQ(ierr);
  ierr = VecDuplicateVecs(vh, 2, &vuvbar); CHKERRQ(ierr);

  ierr = VecDuplicateVecs(vh, nWork2d, &vWork2d); CHKERRQ(ierr);
  ierr = VecDuplicateVecs(vu, nWork3d, &vWork3d); CHKERRQ(ierr);

  ierr = DACreateGlobalVector(grid.da2, &g2); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(grid.da3, &g3); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(grid.da3b, &g3b); CHKERRQ(ierr);

  const PetscInt M = 2 * grid.p->Mx * grid.p->My;
  ierr = MatCreateMPIAIJ(grid.com, PETSC_DECIDE, PETSC_DECIDE, M, M,
                         13, PETSC_NULL, 13, PETSC_NULL,
                         &MacayealStiffnessMatrix); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com, PETSC_DECIDE, M, &MacayealX); CHKERRQ(ierr);
  ierr = VecDuplicate(MacayealX, &MacayealRHS); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, M, &MacayealXLocal);
  ierr = VecScatterCreate(MacayealX, PETSC_NULL, MacayealXLocal, PETSC_NULL,
                          &MacayealScatterGlobalToLocal); CHKERRQ(ierr);
  ierr = KSPCreate(grid.com, &MacayealKSP); CHKERRQ(ierr);

  createVecs_done = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceModel::destroyVecs() {
  PetscErrorCode ierr;

  ierr = bedDefCleanup(); CHKERRQ(ierr);
  ierr = PDDCleanup(); CHKERRQ(ierr);

  ierr = VecDestroy(vu); CHKERRQ(ierr);
  ierr = VecDestroy(vv); CHKERRQ(ierr);
  ierr = VecDestroy(vw); CHKERRQ(ierr);
  ierr = VecDestroy(vSigma); CHKERRQ(ierr);
  ierr = VecDestroy(vT); CHKERRQ(ierr);
  ierr = VecDestroy(vtau); CHKERRQ(ierr);
  ierr = VecDestroy(vgs); CHKERRQ(ierr);

  ierr = VecDestroy(vh); CHKERRQ(ierr);
  ierr = VecDestroy(vH); CHKERRQ(ierr);
  ierr = VecDestroy(vbed); CHKERRQ(ierr);
  ierr = VecDestroy(vAccum); CHKERRQ(ierr);
  ierr = VecDestroy(vTs); CHKERRQ(ierr);
  ierr = VecDestroy(vMask); CHKERRQ(ierr);
  ierr = VecDestroy(vGhf); CHKERRQ(ierr);
  ierr = VecDestroy(vubar); CHKERRQ(ierr);
  ierr = VecDestroy(vvbar); CHKERRQ(ierr);
  ierr = VecDestroy(vub); CHKERRQ(ierr);
  ierr = VecDestroy(vvb); CHKERRQ(ierr);
  ierr = VecDestroy(vRb); CHKERRQ(ierr);
  ierr = VecDestroy(vHmelt); CHKERRQ(ierr);
  ierr = VecDestroy(vbasalMeltRate); CHKERRQ(ierr);
  ierr = VecDestroy(vuplift); CHKERRQ(ierr);
  ierr = VecDestroy(vdHdt); CHKERRQ(ierr);
  ierr = VecDestroy(vbeta); CHKERRQ(ierr);
  ierr = VecDestroy(vtauc); CHKERRQ(ierr);
  ierr = VecDestroy(vLongitude); CHKERRQ(ierr);
  ierr = VecDestroy(vLatitude); CHKERRQ(ierr);

  ierr = VecDestroyVecs(vuvbar, 2); CHKERRQ(ierr);
  ierr = VecDestroyVecs(vDf, 2); CHKERRQ(ierr);
  ierr = VecDestroyVecs(vWork3d, nWork3d); CHKERRQ(ierr);
  ierr = VecDestroyVecs(vWork2d, nWork2d); CHKERRQ(ierr);

  ierr = VecDestroy(g2); CHKERRQ(ierr);
  ierr = VecDestroy(g3); CHKERRQ(ierr);
  ierr = VecDestroy(g3b); CHKERRQ(ierr);

  ierr = KSPDestroy(MacayealKSP); CHKERRQ(ierr);
  ierr = MatDestroy(MacayealStiffnessMatrix); CHKERRQ(ierr);
  ierr = VecDestroy(MacayealX); CHKERRQ(ierr);
  ierr = VecDestroy(MacayealRHS); CHKERRQ(ierr);
  ierr = VecDestroy(MacayealXLocal); CHKERRQ(ierr);
  ierr = VecScatterDestroy(MacayealScatterGlobalToLocal); CHKERRQ(ierr);

  return 0;
}


void IceModel::setTimeStepYears(PetscScalar y) {
  dt = y * secpera;
  doAdaptTimeStep = PETSC_FALSE;
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
  VecSet(vtau, d*secpera);
}

void IceModel::setShowViewers(PetscTruth show_viewers) {
  showViewers = show_viewers;
}

void IceModel::setDoMassConserve(PetscTruth do_mb) {
  doMassConserve = do_mb;
}

void IceModel::setDoTemp(PetscTruth do_temp) {
  doTemp = do_temp;
}

void IceModel::setIncludeBMRinContinuity(PetscTruth includeit) {
  includeBMRinContinuity = includeit;
}

void IceModel::setDoGrainSize(PetscTruth do_gs) {
  doGrainSize = do_gs;
}

void IceModel::setDoBedDef(PetscTruth do_bd) {
  doBedDef = do_bd;
}

void IceModel::setDoBedIso(PetscTruth do_iso) {
  doBedIso = do_iso;
}

void IceModel::setIsDrySimulation(PetscTruth is_dry) {
  isDrySimulation = is_dry;
}

void IceModel::setAllGMaxVelocities(PetscScalar uvw_for_cfl) {
  gmaxu=uvw_for_cfl;
  gmaxv=uvw_for_cfl;
  gmaxw=uvw_for_cfl;
}

void IceModel::setThermalBedrock(PetscTruth tb) {
  thermalBedrock = tb;
}

void IceModel::setOceanKill(PetscTruth ok) {
  doOceanKill = ok;
}

void IceModel::setUseMacayealVelocity(PetscTruth umv) {
  useMacayealVelocity = umv;
}

void IceModel::setDoSuperpose(PetscTruth ds) {
  doSuperpose = ds;
}

void IceModel::setConstantNuForMacAyeal(PetscScalar nu) {
  useConstantNuForMacAyeal = PETSC_TRUE;
  constantNuForMacAyeal = nu;
}

void IceModel::setRegularizingVelocitySchoof(PetscScalar rvS) {
  regularizingVelocitySchoof = rvS;
}

void IceModel::setRegularizingLengthSchoof(PetscScalar rLS) {
  regularizingLengthSchoof = rLS;
}

void IceModel::setMacayealEpsilon(PetscScalar meps) {
  macayealEpsilon = meps;
}

void IceModel::setMacayealRelativeTolerance(PetscScalar mrc) {
  macayealRelativeTolerance = mrc;
}

void IceModel::setEnhancementFactor(PetscScalar e) {
  enhancementFactor = e;
}

void IceModel::setMuSliding(PetscScalar mu) {
  muSliding = mu;
}

void IceModel::setGSIntervalYears(PetscScalar years) {
  gsIntervalYears = years;
}

void IceModel::setBedDefIntervalYears(PetscScalar years) {
  bedDefIntervalYears = years;
}

void IceModel::setIsothermalFlux(PetscTruth use) {
  useIsothermalFlux = use;
}

void IceModel::setNoSpokes(PetscInt level) {
  noSpokesLevel = level;
}

void IceModel::setIsothermalFlux(PetscTruth use, PetscScalar n, PetscScalar A) {
  setIsothermalFlux(use);
  isothermalFlux_n_exponent = n;
  isothermalFlux_A_softness = A;
}


PetscTruth IceModel::isInitialized() const {
  return initialized_p;
}


PetscErrorCode IceModel::updateSurfaceElevationAndMask() {
  // should be called when either ice thickness or bed elevation change, to 
  // maintain consistency of geometry
  PetscErrorCode ierr;
  PetscScalar **h, **bed, **H, **mask, ***T;
  const int MASK_GROUNDED_TO_DETERMINE = 999;

  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // take this opportunity to check that H[i][j] >= 0
      if (H[i][j] < 0) {
        SETERRQ2(1,"Thickness negative at point i=%d, j=%d",i,j);
      }

      const PetscScalar hgrounded = bed[i][j] + H[i][j];

      if (isDrySimulation == PETSC_TRUE) {
        h[i][j] = hgrounded;
        // don't update mask; potentially one would want to do MacAyeal
        //   dragging ice shelf in dry case and/or ignor mean sea level elevation
      } else {

        const PetscScalar hfloating = (1-ice.rho/ocean.rho) * H[i][j];
        if (modMask(mask[i][j]) == MASK_FLOATING) {
          // check whether you are actually floating or grounded
          if (hgrounded > hfloating+1.0) {
            mask[i][j] = MASK_GROUNDED_TO_DETERMINE;
            h[i][j] = hgrounded; // actually grounded so update h
          } else {
            h[i][j] = hfloating; // actually floating so update h
          }
        } else { // deal with grounded ice according to mask
          if (hgrounded > hfloating-1.0) {
            h[i][j] = hgrounded; // actually grounded so update h
          } else {
            mask[i][j] = MASK_FLOATING;
            h[i][j] = hfloating; // actually floating so update h
          }
        }

        if (intMask(mask[i][j]) == MASK_GROUNDED_TO_DETERMINE) {
          if (useMacayealVelocity != PETSC_TRUE) {
            mask[i][j] = MASK_SHEET;
          } else {
            // if frozen to bed or essentially frozen to bed then make it SHEET
            if (T[i][j][0] + ice.beta_CC_grad * H[i][j] 
                         < DEFAULT_MIN_TEMP_FOR_SLIDING) { 
              mask[i][j] = MASK_SHEET;
            } else {
              // determine type of grounded ice by vote-by-neighbors
              //   (BOX stencil neighbors!):
              const PetscScalar neighmasksum = 
                modMask(mask[i-1][j+1]) + modMask(mask[i][j+1]) + modMask(mask[i+1][j+1]) +
                modMask(mask[i-1][j])   +                       + modMask(mask[i+1][j])  +
                modMask(mask[i-1][j-1]) + modMask(mask[i][j-1]) + modMask(mask[i+1][j-1]);
              // make SHEET if either all neighbors are SHEET or at most one is 
              //   DRAGGING; if any are floating then ends up DRAGGING:
              if (neighmasksum <= (7*MASK_SHEET + MASK_DRAGGING + 0.1)) { 
                mask[i][j] = MASK_SHEET;
              } else { // otherwise make DRAGGING
                mask[i][j] = MASK_DRAGGING;
              }
            }
          }
        }
        
      }

    }
  }

  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vh, INSERT_VALUES, vh); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::massBalExplicitStep() {
  const PetscScalar   dx = grid.p->dx, dy = grid.p->dy;
  PetscErrorCode ierr;
  PetscScalar **H, **Hnew, **uvbar[2];
  PetscScalar **u, **v, **accum, **basalMeltRate, **mask;
  Vec vHnew = vWork2d[0];

#if (MARGIN_TRICK)
  ierr = verbPrintf(4,grid.com,"  MARGIN_TRICK massBalExplicitStep() ..."); CHKERRQ(ierr);
#endif

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = VecCopy(vH, vHnew); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vHnew, &Hnew); CHKERRQ(ierr);

  PetscScalar icecount = 0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0)  icecount++;
      PetscScalar divQ;
      if (intMask(mask[i][j]) == MASK_SHEET) { 
        // staggered grid Div(Q) for Q = D grad h
#if (MARGIN_TRICK_TWO)
        /* switching between MARGIN_TRICK_TWO=1 and MARGIN_TRICK_TWO=0 versions seems to make no 
           significant difference
           in tests B, C, and G but there does seem to be a difference in EISMINT II exp F */
        const PetscScalar He = (H[i+1][j] > 0.0) ? 0.5*(H[i][j] + H[i+1][j])
                            : ( (H[i-1][j] > 0.0) ? 1.5*H[i-1][j] - 0.5*H[i][j] : H[i][j] );
        const PetscScalar Hw = (H[i-1][j] > 0.0) ? 0.5*(H[i-1][j] + H[i][j])
                            : ( (H[i+1][j] > 0.0) ? 1.5*H[i+1][j] - 0.5*H[i][j] : H[i][j] );
        const PetscScalar Hn = (H[i][j+1] > 0.0) ? 0.5*(H[i][j] + H[i][j+1])
                            : ( (H[i][j-1] > 0.0) ? 1.5*H[i][j-1] - 0.5*H[i][j] : H[i][j] );
        const PetscScalar Hs = (H[i][j-1] > 0.0) ? 0.5*(H[i][j-1] + H[i][j])
                            : ( (H[i][j+1] > 0.0) ? 1.5*H[i][j+1] - 0.5*H[i][j] : H[i][j] );
        divQ =  (uvbar[0][i][j] * He - uvbar[0][i-1][j] * Hw) / dx
              + (uvbar[1][i][j] * Hn - uvbar[1][i][j-1] * Hs) / dy;
#else
        divQ =
          (uvbar[0][i][j] * 0.5*(H[i][j] + H[i+1][j])
           - uvbar[0][i-1][j] * 0.5*(H[i-1][j] + H[i][j])) / dx
          + (uvbar[1][i][j] * 0.5*(H[i][j] + H[i][j+1])
             - uvbar[1][i][j-1] * 0.5*(H[i][j-1] + H[i][j])) / dy;
#endif
      } else { // upwinded, regular grid Div(Q), for Q = Ubar H, computed as
               //     Div(Q) = U . grad H + Div(U) H
               // note the CFL for "U . grad H" part of upwinding is checked; see
               // broadcastMacayealVelocity() and determineTimeStep()
        divQ =
          u[i][j] * (u[i][j] < 0 ? H[i+1][j]-H[i][j] : H[i][j]-H[i-1][j]) / dx
          + v[i][j] * (v[i][j] < 0 ? H[i][j+1]-H[i][j] : H[i][j]-H[i][j-1]) / dy
          + H[i][j] * ((u[i+1][j]-u[i-1][j])/(2*dx) + (v[i][j+1]-v[i][j-1])/(2*dy));
      }

      Hnew[i][j] += (accum[i][j] - divQ) * dt;
      if (includeBMRinContinuity == PETSC_TRUE) {
         Hnew[i][j] -= capBasalMeltRate(basalMeltRate[i][j]) * dt;
      }

      if (Hnew[i][j] < 0)
        // apply free boundary rule: negative thickness becomes zero
        Hnew[i][j] = 0.0;
      if ( (doOceanKill == PETSC_TRUE) 
           && (intMask(mask[i][j]) == MASK_FLOATING_OCEAN0) )
        // force zero at ocean; "accumulation-zone" b.c.
        Hnew[i][j] = 0.0;
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vbasalMeltRate, &basalMeltRate); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vHnew, &Hnew); CHKERRQ(ierr);

  // compute dH/dt (thickening rate) for viewing and for saving at end; only diagnostic
  ierr = VecWAXPY(vdHdt, -1, vH, vHnew); CHKERRQ(ierr);
  ierr = VecScale(vdHdt,1.0/dt); CHKERRQ(ierr);

  // average value of dH/dt; also d(volume)/dt
  PetscScalar gicecount;
  ierr = PetscGlobalSum(&icecount, &gicecount, grid.com); CHKERRQ(ierr);
  ierr = DALocalToGlobal(grid.da2, vdHdt, INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr = VecSum(g2, &gdHdtav); CHKERRQ(ierr);
  dvoldt = gdHdtav * grid.p->dx * grid.p->dy;  // m^3/s
  gdHdtav = gdHdtav / gicecount; // m/s

  // finally copy vHnew into vH (and communicate ghosted values at same time)
  ierr = DALocalToLocalBegin(grid.da2, vHnew, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vHnew, INSERT_VALUES, vH); CHKERRQ(ierr);

  // update h and mask
  ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::run() {
  PetscErrorCode  ierr;

  ierr = verbPrintf(2,grid.com, "$$$$$"); CHKERRQ(ierr);
  ierr = summaryPrintLine(PETSC_TRUE,PETSC_TRUE, 0.0, 0.0, 0, ' ', 0.0, 0.0, 0.0, 0.0, 0.0); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, "$$$$$"); CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep
  tempskipCountDown = 0;
  ierr = summary(true,true); CHKERRQ(ierr);  // report starting state

  dtTempAge = 0.0;
  // main loop for time evolution
  for (PetscScalar year = startYear; year < endYear; year += dt/secpera) {
    dt_force = -1.0;
    maxdt_temporary = -1.0;
    ierr = additionalAtStartTimestep(); CHKERRQ(ierr);  // might set dt_force,maxdt_temporary
    
    // update basal till yield stress if appropriate; will modify and communicate mask
    if (doPlasticTill == PETSC_TRUE) {
      ierr = updateYieldStressFromHmelt();  CHKERRQ(ierr);
    }
    
    // compute PDD; generates net accumulation, with possible ablation area, using snow accumulation
    // might set maxdt_temporary 
    if ((doPDD == PETSC_TRUE) && IsIntegralYearPDD()) {
      ierr = updateNetAccumFromPDD();  CHKERRQ(ierr);
    }

    // compute bed deformation, which only depends on current thickness and bed elevation
    if (doBedDef == PETSC_TRUE) {
      ierr = bedDefStepIfNeeded(); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }

    // always do vertically-average velocity calculation; only update velocities at depth if
    // needed for temp and age calculation
    bool updateAtDepth = (tempskipCountDown == 0);
    ierr = velocity(updateAtDepth); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com, updateAtDepth ? "v" : "V" ); CHKERRQ(ierr);

    // now that velocity field is up to date, compute grain size
    if (doGrainSize == PETSC_TRUE) {
      ierr = updateGrainSizeIfNeeded(); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }
    
    // adapt time step using velocities and diffusivity, ..., just computed
    bool useCFLforTempAgeEqntoGetTimestep = (doTemp == PETSC_TRUE);
    ierr = determineTimeStep(useCFLforTempAgeEqntoGetTimestep); CHKERRQ(ierr);
    dtTempAge += dt;
    grid.p->year += dt / secpera;  // adopt it
    // IceModel::dt,dtTempAge,grid.p->year are now set correctly according to
    //    mass-continuity-eqn-diffusivity criteria, CFL criteria, and other 
    //    criteria from derived class additionalAtStartTimestep(), and from 
    //    "-tempskip" mechanism

    // ierr = PetscPrintf(PETSC_COMM_SELF,
    //           "\n[rank=%d, year=%f, dt=%f, startYear=%f, endYear=%f]",
    //           grid.rank, grid.p->year, dt/secpera, startYear, endYear);
    //        CHKERRQ(ierr);
    
    bool tempAgeStep = (updateAtDepth && (doTemp == PETSC_TRUE));
    if (tempAgeStep) { // do temperature and age
      allowAboveMelting = PETSC_FALSE;
      ierr = temperatureAgeStep(); CHKERRQ(ierr);
      dtTempAge = 0.0;
      ierr = verbPrintf(2,grid.com, "t"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }
    
    if (doMassConserve == PETSC_TRUE) {
      ierr = massBalExplicitStep(); CHKERRQ(ierr);
      if ((doTempSkip == PETSC_TRUE) && (tempskipCountDown > 0))
        tempskipCountDown--;
      ierr = verbPrintf(2,grid.com, "f"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
    }
    
    ierr = summary(tempAgeStep,true); CHKERRQ(ierr);

//    ierr = verbPrintf(2,grid.com, " tempskipCountDown=%d, dt_from_cfl=%10.5e, dt_from_diffus=%10.5e, CFLmaxdt=%10.5e\n",
//                      tempskipCountDown, dt_from_cfl, dt_from_diffus, CFLmaxdt); CHKERRQ(ierr);
     
    ierr = updateViewers(); CHKERRQ(ierr);

    ierr = additionalAtEndTimestep(); CHKERRQ(ierr);

    if (endOfTimeStepHook() != 0) break;
  }
  
  return 0;
}


PetscErrorCode IceModel::diagnosticRun() {
  PetscErrorCode  ierr;

  // print out some stats about input state
  ierr = verbPrintf(2,grid.com, "$$$$$"); CHKERRQ(ierr);
  ierr = summaryPrintLine(PETSC_TRUE,PETSC_TRUE, 0.0, 0.0, 0, ' ', 0.0, 0.0, 0.0, 0.0, 0.0); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, "$$$$$"); CHKERRQ(ierr);
  adaptReasonFlag = ' '; // no reason for no timestep
  tempskipCountDown = 0;
  ierr = summary(true,true); CHKERRQ(ierr);  // report starting state

  // update basal till yield stress if appropriate; will modify and communicate mask
  if (doPlasticTill == PETSC_TRUE) {
    ierr = updateYieldStressFromHmelt();  CHKERRQ(ierr);
  }

  ierr = velocity(true); CHKERRQ(ierr);  // compute velocities (at depth); this is the point
  
  // update viewers and pause for a chance to view
  ierr = updateViewers(); CHKERRQ(ierr);
  PetscInt    pause_time = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, PETSC_NULL); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,"pausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
  ierr = PetscSleep(pause_time); CHKERRQ(ierr);
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


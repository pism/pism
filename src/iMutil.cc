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

#include <cstring>
#include <ctime>
#include "iceModel.hh"
#include "pism_signal.h"

#include <petscsys.h>
#include <stdarg.h>
#include <stdlib.h>
#include <petscfix.h>

// verbosity level version of PetscPrintf: print according to whether 
// (thresh <= IceModel::verbosityLevel), in which case print, or 
// (thresh > verbosityLevel) in which case no print
//
//   level  option        meaning
//   -----  ------        -------
//   0      -verbose 0    never print to std out AT ALL!
//
//   1      -verbose 1    less verbose than default: thresh must be 1 to print
//
//   2     [-verbose 2]   default
//
//   3      -verbose      somewhat verbose
//         [-verbose 3]   
//   4      -vverbose     fairly verbose
//         [-verbose 4]
//   5      -vvverbose    very verbose: if level this high then (thresh <= level) 
//         [-verbose 5]       always, so print everything
//
// note: 1 <= thresh <= 5  enforced in verbPrintf() below

PetscErrorCode setVerbosityLevel(PetscInt level) {
  if ((level < 0) || (level > 5)) {
     SETERRQ(1,"verbosity level invalid");
  }
  verbosityLevel = level;
  return 0;  
}


PetscErrorCode verbosityLevelFromOptions() {
  // verbosity options: more info to standard out
  PetscErrorCode ierr;
  PetscInt     myverbosityLevel;
  PetscTruth   verbose, verbosityLevelSet;
  PetscInt     DEFAULT_VERBOSITY_LEVEL = 2;
  
  ierr = setVerbosityLevel(DEFAULT_VERBOSITY_LEVEL);  
  ierr = PetscOptionsGetInt(PETSC_NULL, "-verbose", &myverbosityLevel, &verbosityLevelSet); CHKERRQ(ierr);
  if (verbosityLevelSet == PETSC_TRUE) {
    ierr = setVerbosityLevel(myverbosityLevel);
  } else {
    ierr = PetscOptionsHasName(PETSC_NULL, "-verbose", &verbose); CHKERRQ(ierr);
    if (verbose == PETSC_TRUE)   ierr = setVerbosityLevel(3);
  }
  ierr = PetscOptionsHasName(PETSC_NULL, "-vverbose", &verbose); CHKERRQ(ierr);
  if (verbose == PETSC_TRUE)   ierr = setVerbosityLevel(4);
  ierr = PetscOptionsHasName(PETSC_NULL, "-vvverbose", &verbose); CHKERRQ(ierr);
  if (verbose == PETSC_TRUE)   ierr = setVerbosityLevel(5);
  return 0;
}


PetscErrorCode verbPrintf(const int thresh, 
                          MPI_Comm comm,const char format[],...)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  size_t         len;
  char           *nformat,*sub1,*sub2;
  PetscReal      value;

  extern FILE *petsc_history;

  if ((thresh < 1) || (thresh > 5)) { SETERRQ(1,"invalid threshold in verbPrintf()"); }

  PetscFunctionBegin;
  if (!comm) comm = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  if (!rank && ((verbosityLevel >= thresh) || petsc_history) ) {
    va_list Argp;
    va_start(Argp,format);

    ierr = PetscStrstr(format,"%A",&sub1);CHKERRQ(ierr);
    if (sub1) {
      ierr = PetscStrstr(format,"%",&sub2);CHKERRQ(ierr);
      if (sub1 != sub2) SETERRQ(PETSC_ERR_ARG_WRONG,"%%A format must be first in format string");
      ierr    = PetscStrlen(format,&len);CHKERRQ(ierr);
      ierr    = PetscMalloc((len+16)*sizeof(char),&nformat);CHKERRQ(ierr);
      ierr    = PetscStrcpy(nformat,format);CHKERRQ(ierr);
      ierr    = PetscStrstr(nformat,"%",&sub2);CHKERRQ(ierr);
      sub2[0] = 0;
      value   = (double)va_arg(Argp,double);
      if (PetscAbsReal(value) < 1.e-12) {
        ierr    = PetscStrcat(nformat,"< 1.e-12");CHKERRQ(ierr);
      } else {
        ierr    = PetscStrcat(nformat,"%g");CHKERRQ(ierr);
        va_end(Argp);
        va_start(Argp,format);
      }
      ierr    = PetscStrcat(nformat,sub1+2);CHKERRQ(ierr);
    } else {
      nformat = (char*)format;
    }
    if (verbosityLevel >= thresh) {
      ierr = PetscVFPrintf(PETSC_STDOUT,nformat,Argp);CHKERRQ(ierr);
    }
    if (petsc_history) { // always print to history
      ierr = PetscVFPrintf(petsc_history,nformat,Argp);CHKERRQ(ierr);
    }
    va_end(Argp);
    if (sub1) {ierr = PetscFree(nformat);CHKERRQ(ierr);}
  }
  PetscFunctionReturn(0);
}


PetscErrorCode IceModel::computeFlowUbarStats
                      (PetscScalar *gUbarmax, PetscScalar *gUbarSIAav,
                       PetscScalar *gUbarstreamav, PetscScalar *gUbarshelfav,
                       PetscScalar *gicegridfrac, PetscScalar *gSIAgridfrac,
                       PetscScalar *gstreamgridfrac, PetscScalar *gshelfgridfrac) {
  // NOTE:  Assumes IceModel::vubar, vvbar, vu, vv holds correct 
  // and up-to-date values of velocities

  PetscErrorCode ierr;

  PetscScalar **H, **ubar, **vbar, **mask;
  PetscScalar Ubarmax = 0.0, UbarSIAsum = 0.0, Ubarstreamsum = 0.0,
              Ubarshelfsum = 0.0, icecount = 0.0, SIAcount = 0.0, shelfcount = 0.0;

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
        icecount += 1.0;
        const PetscScalar Ubarmag 
                           = sqrt(PetscSqr(ubar[i][j]) + PetscSqr(vbar[i][j]));
        Ubarmax = PetscMax(Ubarmax, Ubarmag);
        if (intMask(mask[i][j]) == MASK_SHEET) {
          SIAcount += 1.0;
          UbarSIAsum += Ubarmag;
        } else if (modMask(mask[i][j]) == MASK_FLOATING) {
          shelfcount += 1.0;
          Ubarshelfsum += Ubarmag;
        } else if (intMask(mask[i][j]) == MASK_DRAGGING) {
          // streamcount = icecount - SIAcount - shelfcount
          Ubarstreamsum += Ubarmag;
        } else {
          SETERRQ(1,"should not reach here!");
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&Ubarmax, gUbarmax, grid.com); CHKERRQ(ierr);
  
  // get global sums
  PetscScalar gicecount, gSIAcount, gshelfcount, gstreamcount;
  ierr = PetscGlobalSum(&icecount, &gicecount, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&SIAcount, &gSIAcount, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&shelfcount, &gshelfcount, grid.com); CHKERRQ(ierr);
  gstreamcount = gicecount - gSIAcount - gshelfcount;

  // really getting sums here (not yet averages)
  ierr = PetscGlobalSum(&UbarSIAsum, gUbarSIAav, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&Ubarshelfsum, gUbarshelfav, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&Ubarstreamsum, gUbarstreamav, grid.com); CHKERRQ(ierr);

  if (gSIAcount > 0.0) {
    *gUbarSIAav = *gUbarSIAav / gSIAcount;
  } else  *gUbarSIAav = 0.0;
  if (gshelfcount > 0.0) {
    *gUbarshelfav = *gUbarshelfav / gshelfcount;
  } else  *gUbarshelfav = 0.0;
  if (gstreamcount > 0.0) {
    *gUbarstreamav = *gUbarstreamav / gstreamcount;
  } else  *gUbarstreamav = 0.0;

  // finally make these actual fractions
  if (gicecount > 0.0) {
    *gSIAgridfrac = gSIAcount / gicecount;
    *gshelfgridfrac = gshelfcount / gicecount;
    *gstreamgridfrac = gstreamcount / gicecount;
  } else {
    *gSIAgridfrac = 0.0;
    *gshelfgridfrac = 0.0;
    *gstreamgridfrac = 0.0;
  }
 
  *gicegridfrac = gicecount / ((PetscScalar) (grid.p->Mx * grid.p->My));
  return 0;
}


PetscErrorCode IceModel::computeMaxDiffusivity(bool updateDiffusViewer) {
  // NOTE:  Assumes IceModel::vubar, vvbar holds correct 
  // and up-to-date values of velocities

  PetscErrorCode ierr;

  PetscScalar **h, **H, **ubar, **vbar, **D, **mask;
  PetscScalar Dmax = 0.0;

  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &D); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
        if (intMask(mask[i][j]) == MASK_SHEET) {
          const PetscScalar h_x=(h[i+1][j]-h[i-1][j])/(2.0*grid.p->dx);
          const PetscScalar h_y=(h[i][j+1]-h[i][j-1])/(2.0*grid.p->dy);
          const PetscScalar alpha = sqrt(PetscSqr(h_x) + PetscSqr(h_y));
          // note: when basal sliding is proportional to surface slope, as
          // it usually will be when sliding occurs in a MASK_SHEET area, then
          //    D = H Ubar / alpha
          // is the correct formula; note division by zero is avoided by
          // addition to alpha
          const PetscScalar Ubarmag
                            = sqrt(PetscSqr(ubar[i][j]) + PetscSqr(vbar[i][j]));
          const PetscScalar d =
               H[i][j] * Ubarmag/(alpha + DEFAULT_ADDED_TO_SLOPE_FOR_DIFF_IN_ADAPTIVE);
          if (d > Dmax) Dmax = d;
          D[i][j] = d;
        } else {
          D[i][j] = 0.0; // no diffusivity in non-SIA regions
        }
      } else {
        D[i][j] = 0.0; // no diffusivity if no ice
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &D); CHKERRQ(ierr);

  if (updateDiffusViewer && (diffusView != PETSC_NULL)) { // -f option: view diffusivity (m^2/s)
    ierr = DALocalToGlobal(grid.da2, vWork2d[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, diffusView); CHKERRQ(ierr);
  }

  ierr = PetscGlobalMax(&Dmax, &gDmax, grid.com); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::adaptTimeStepDiffusivity() {
  // note computeMaxDiffusivity() must be called before this to set gDmax
  // note that adapt_ratio * 2 is multiplied by dx^2/(2*maxD) so 
  // dt <= adapt_ratio * dx^2/maxD (if dx=dy)
  // reference: Morton & Mayers 2nd ed. pp 62--63
  const PetscScalar  
          gridfactor = 1.0/(grid.p->dx*grid.p->dx) + 1.0/(grid.p->dy*grid.p->dy);
  dt_from_diffus = adaptTimeStepRatio
                     * 2 / ((gDmax + DEFAULT_ADDED_TO_GDMAX_ADAPT) * gridfactor);
  if ((doTempSkip == PETSC_TRUE) && (tempskipCountDown == 0)) {
    const PetscScalar  conservativeFactor = 0.8;
    // typically "dt" in next line is from CFL, but might be from other, e.g. maxdt
    tempskipCountDown = (PetscInt) floor(conservativeFactor * (dt / dt_from_diffus));
    tempskipCountDown = (tempskipCountDown > tempskipMax) ? tempskipMax : tempskipCountDown;
  } // if tempskipCountDown > 0 then it will get decremented at the mass balance step
  if (dt_from_diffus < dt) {
    dt = dt_from_diffus;
    adaptReasonFlag = 'd';
  }
  return 0;
}


PetscErrorCode IceModel::determineTimeStep(const bool doTemperatureCFL) {
  PetscErrorCode ierr;

  if ( (diffusView != PETSC_NULL) 
       || ( (doAdaptTimeStep == PETSC_TRUE) && (doMassConserve == PETSC_TRUE) ) ) {
    ierr = computeMaxDiffusivity(true); CHKERRQ(ierr);
  }
  const PetscScalar timeToEnd = (endYear-grid.p->year) * secpera;
  if (dt_force > 0.0) {
    dt = dt_force; // override usual dt mechanism
    adaptReasonFlag = 'f';
    if (timeToEnd < dt) {
      dt = timeToEnd;
      adaptReasonFlag = 'e';
    }
  } else {
    dt = maxdt;
    adaptReasonFlag = 'm';
    if ((doAdaptTimeStep == PETSC_TRUE) && (doTemp == PETSC_TRUE)
        && doTemperatureCFL) {
      // CFLmaxdt is set by computeMax3DVelocities() in call to velocity() iMvelocity.cc
      dt_from_cfl = CFLmaxdt;
      if (dt_from_cfl < dt) {
        dt = dt_from_cfl;
        adaptReasonFlag = 'c';
      }
    } 
    if ((doAdaptTimeStep == PETSC_TRUE) && (doMassConserve == PETSC_TRUE)
        && (useMacayealVelocity)) {
      // CFLmaxdt2D is set by broadcastMacayealVelocity()
      if (CFLmaxdt2D < dt) {
        dt = CFLmaxdt2D;
        adaptReasonFlag = 'u';
      }
    }
    if ((doAdaptTimeStep == PETSC_TRUE) && (doMassConserve == PETSC_TRUE)) {
      // note: if doTempSkip then tempskipCountDown = floor(dt_from_cfl/dt_from_diffus)
      ierr = adaptTimeStepDiffusivity(); CHKERRQ(ierr); // might set adaptReasonFlag = 'd'
    }
    if ((maxdt_temporary > 0.0) && (maxdt_temporary < dt)) {
      dt = maxdt_temporary;
      adaptReasonFlag = 't';
    }
    if (timeToEnd < dt) {
      dt = timeToEnd;
      adaptReasonFlag = 'e';
    }
    if ((adaptReasonFlag == 'm') || (adaptReasonFlag == 't') || (adaptReasonFlag == 'e')) {
      if (tempskipCountDown > 1) tempskipCountDown = 1; 
    }
  }    

/*
verbPrintf(1,grid.com,
   "\n   leaving determineTimeStep();\n     endYear, grid.p->year, dt_force = %f,%f,%e\n",
   endYear, grid.p->year, dt_force);
verbPrintf(1,grid.com,
   "     maxdt, CFLmaxdt, maxdt_temporary, gDmax = %e,%e,%e,%e\n",
   maxdt, CFLmaxdt, maxdt_temporary, gDmax);
verbPrintf(1,grid.com,
   "     dt, dt/secpera = %e, %f\n",
   dt, dt/secpera);
*/
  return 0;
}


PetscErrorCode IceModel::volumeArea(PetscScalar& gvolume, PetscScalar& garea,
                                    PetscScalar& gvolSIA, PetscScalar& gvolstream, 
                                    PetscScalar& gvolshelf) {
  // returns area in units of km^2 and volume in km^3
  // though slightly less efficient when used by summaryEismint, which has to look at vH anyway,
  //   it is clearer to have this obvious ability as a procedure
  PetscErrorCode  ierr;
  PetscScalar     **H, **mask;
  PetscScalar     volume=0.0, area=0.0, volSIA=0.0, volstream=0.0, volshelf=0.0;
  
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  const PetscScalar   a = grid.p->dx * grid.p->dy * 1e-3 * 1e-3; // area unit (km^2)
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        area += a;
        const PetscScalar dv = a * H[i][j] * 1e-3;
        volume += dv;
        if (intMask(mask[i][j]) == MASK_SHEET)   volSIA += dv;
        else if (intMask(mask[i][j]) == MASK_DRAGGING)   volstream += dv;
        else if (modMask(mask[i][j]) == MASK_FLOATING)   volshelf += dv;
      }
    }
  }  
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volume, &gvolume, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&area, &garea, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volSIA, &gvolSIA, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volstream, &gvolstream, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volshelf, &gvolshelf, grid.com); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::summary(bool tempAndAge, bool useHomoTemp) {
  PetscErrorCode  ierr;
  PetscScalar     **H, ***T, ***tau;
  PetscScalar     melt, divideH, divideT, orig;
  PetscScalar     gmelt, gdivideH, gdivideT, gorig, gvolume, garea;
  PetscScalar     gvolSIA, gvolstream, gvolshelf;
  PetscScalar     meltfrac, origfrac;
  
  ierr = volumeArea(gvolume,garea,gvolSIA, gvolstream, gvolshelf); CHKERRQ(ierr);
  
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  divideH = 0; 
  if (tempAndAge || (verbosityLevel >= 3)) {
    ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vtau, &tau); CHKERRQ(ierr);
    melt = 0; divideT = 0; orig = 0;
  }
  const PetscScalar   a = grid.p->dx * grid.p->dy * 1e-3 * 1e-3; // area unit (km^2)
  const PetscScalar   v = a * grid.p->dz * 1e-3;  // volume unit (km^3)
  const PetscScalar   currtime = grid.p->year * secpera;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (i == (grid.p->Mx - 1)/2 && j == (grid.p->My - 1)/2) {
        divideH = H[i][j];
      }
      if (tempAndAge || (verbosityLevel >= 3)) {
        if (H[i][j] > 0) {
          if (useHomoTemp) {
//            if (T[i][j][0] + ice.beta_CC_grad * H[i][j] >= ice.meltingTemp)
            if (T[i][j][0] + ice.beta_CC_grad * H[i][j] >= DEFAULT_MIN_TEMP_FOR_SLIDING)
              melt += a;
          } else {
            if (T[i][j][0] >= ice.meltingTemp)
              melt += a;
          }
        }
        if (i == (grid.p->Mx - 1)/2 && j == (grid.p->My - 1)/2) {
          divideT = T[i][j][0];
        }
        const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        for (PetscInt k=1; k<=ks; k++) {
          // ice is original if it is at least one year older than current time
          if (tau[i][j][k] > currtime + secpera)
            orig += v;
        }
      }
    }
  }
  
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&divideH, &gdivideH, grid.com); CHKERRQ(ierr);
  if (tempAndAge || (verbosityLevel >= 3)) {
    ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vtau, &tau); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&melt, &gmelt, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&orig, &gorig, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalMax(&divideT, &gdivideT, grid.com); CHKERRQ(ierr);
    if (gvolume>0) origfrac=gorig/gvolume;
    else origfrac=0.0;
    if (garea>0) meltfrac=gmelt/garea;
    else meltfrac=0.0;
  }
 
  // give summary data a la EISMINT II:
  //    year (+ dt[NR]) (years),
  //       [ note 
  //            N = tempskipCountDown
  //            R = on character reason for dt: 
  //            m = maxdt applied, 
  //            e = time to end, 
  //            d = diffusive limit from mass continuity,
  //            c = CFL for temperature equation,
  //            f = forced by derived class,
  //            t = temporarily truncated (e.g. to next integer yr) by derived class,
  //            [space] = no time step taken]
  //    volume (10^6 cubic km),
  //    area (10^6 square km),
  //    basal melt ratio [=(area of base at *absolute*  or *pressure* melting);
  //                      depends on useHomoTemp],
  //    divide thickness (m),
  //    temp at base at divide (K)  (not homologous),
  // NOTE DERIVED CLASSES MAY HAVE OTHER DISPLAYED QUANTITIES
  ierr = summaryPrintLine(PETSC_FALSE,(PetscTruth)tempAndAge,grid.p->year,dt,
                          tempskipCountDown,adaptReasonFlag,
                          gvolume,garea,meltfrac,gdivideH,gdivideT); CHKERRQ(ierr);

  if (CFLviolcount > 0.5) { // report any CFL violations at end of line
    ierr = verbPrintf(2,grid.com,"  [!CFL#=%1.0f]\n",CFLviolcount); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com,"\n"); CHKERRQ(ierr);
  }
  
  if (verbosityLevel >= 3) {
    // show additional info
    PetscScalar Ubarmax, UbarSIAav, Ubarstreamav, Ubarshelfav, icegridfrac,
         SIAgridfrac, streamgridfrac, shelfgridfrac;
    ierr = computeMaxDiffusivity(false); CHKERRQ(ierr); 
    ierr = computeFlowUbarStats(&Ubarmax,
              &UbarSIAav, &Ubarstreamav, &Ubarshelfav, &icegridfrac,
              &SIAgridfrac, &streamgridfrac, &shelfgridfrac); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, 
           "    (volume of ice which is SIA, stream, shelf:  %8.3f,  %8.3f,  %8.3f)\n",
                         gvolSIA/1.0e6, gvolstream/1.0e6, gvolshelf/1.0e6);
                         CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, 
           "  d(volume)/dt of ice (km^3/a):    %11.2f\n",
                         dvoldt*secpera*1.0e-9); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, 
           "  average value of dH/dt (m/a):    %11.5f\n",
                         gdHdtav*secpera); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, 
           "  area percent covered by ice:       %9.4f\n",
                         icegridfrac*100.0); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, 
           "    (area percent ice SIA, stream, shelf:        %8.4f,  %8.4f,  %8.4f)\n",
                         SIAgridfrac*100.0, streamgridfrac*100.0, shelfgridfrac*100.0);
                         CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, 
           "  max diffusivity D on SIA (m^2/s):  %9.3f\n",
                         gDmax); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, 
           "  max |bar U| in all ice (m/a):     %10.3f\n", Ubarmax*secpera); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, 
           "    (av |bar U| in SIA, stream, shelf (m/a):    %9.3f, %9.3f, %9.3f)\n",
           UbarSIAav*secpera, Ubarstreamav*secpera, Ubarshelfav*secpera); CHKERRQ(ierr);
    if (tempAndAge) {
      ierr = PetscPrintf(grid.com, 
           "  maximum |u|,|v|,|w| in ice (m/a): "); CHKERRQ(ierr);
      if (gmaxu == DEFAULT_MAX_VEL_FOR_CFL) {
        ierr = PetscPrintf(grid.com,            "     <N/A>\n"); CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(grid.com,            "%10.3f,%10.3f, %9.3f\n",
        gmaxu*secpera, gmaxv*secpera, gmaxw*secpera); CHKERRQ(ierr);
      }
      ierr = PetscPrintf(grid.com, 
           "  fraction of ice which is original: %9.3f\n",
                         origfrac); CHKERRQ(ierr);
    }
  }

  return 0;
}


PetscErrorCode IceModel::summaryPrintLine(
    const PetscTruth printPrototype, const PetscTruth tempAndAge,
    const PetscScalar year, const PetscScalar dt, 
    const PetscInt tempskipCount, const char adaptReason,
    const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
    const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0) {

  PetscErrorCode ierr;
  if (printPrototype == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,
      "       YEAR (+     STEP[N$]):     VOL    AREA    MELTF     THICK0     TEMP0\n");
  } else {
    ierr = verbPrintf(2,grid.com, "%11.3f (+%9.5f[%d%c]):%8.3f%8.3f%9.3f%11.3f%10.3f",
                       year, dt/secpera, tempskipCount, adaptReason, 
                       volume_kmcube/1.0e6, area_kmsquare/1.0e6, meltfrac,
                       H0,T0); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::additionalAtStartTimestep() {
  // does nothing; allows derived classes to do more per-step computation,
  // reporting and checking, etc.
  return 0;
}


PetscErrorCode IceModel::additionalAtEndTimestep() {
  // does nothing; allows derived classes to do more per-step computation,
  // reporting and checking, etc.
  return 0;
}


PetscErrorCode IceModel::getHorSliceOf3D(Vec v3D, Vec &gslice, PetscInt k) {
  PetscErrorCode ierr;
  PetscScalar    ***val, **slice_val;

  ierr = DAVecGetArray(grid.da2, gslice, &slice_val); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, v3D, &val); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      slice_val[i][j] = val[i][j][k];
    }
  }
  ierr = DAVecRestoreArray(grid.da2, gslice, &slice_val); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, v3D, &val); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::getSurfaceValuesOf3D(Vec v3D, Vec &gsurf) {
  PetscErrorCode ierr;
  PetscScalar    ***val, **H, **surf_val;

  ierr = DAVecGetArray(grid.da2, gsurf, &surf_val); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, v3D, &val); CHKERRQ(ierr);
  const PetscScalar dz = grid.p->dz;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/dz));
      if (ks == grid.p->Mz-1) {  // if at top of grid
        surf_val[i][j] = val[i][j][ks];
      } else {  // general case is to linearly interpolate to z=H[i][j]
        surf_val[i][j] = ((H[i][j] - ks * dz) * val[i][j][ks]
                         + ((ks+1) * dz - H[i][j]) * val[i][j][ks+1]) / dz;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, gsurf, &surf_val); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, v3D, &val); CHKERRQ(ierr);
  return 0;
}



PetscErrorCode IceModel::initFromOptions() {
  PetscErrorCode ierr;
  PetscTruth inFileSet, bootstrapSet, bootstrapSetLegacy;
  char inFile[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-bif", inFile,
                               PETSC_MAX_PATH_LEN, &bootstrapSet);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-bif_legacy", inFile,
                               PETSC_MAX_PATH_LEN, &bootstrapSetLegacy);
  CHKERRQ(ierr);
  
  if (bootstrapSet == PETSC_TRUE) {
    ierr = bootstrapFromFile_netCDF(inFile); CHKERRQ(ierr);
  } else if(bootstrapSetLegacy == PETSC_TRUE) {
    ierr = bootstrapFromFile_netCDF_legacyAnt(inFile); CHKERRQ(ierr);
  } else if (inFileSet == PETSC_TRUE) {
    ierr = initFromFile(inFile); CHKERRQ(ierr);
  }
  
  if (! isInitialized()) {
    SETERRQ(1,"Model has not been initialized.");
  }
  
  ierr = afterInitHook(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::afterInitHook() {
  PetscErrorCode ierr;
  PetscTruth     regridFileSet = PETSC_FALSE;
  char           regridFile[PETSC_MAX_PATH_LEN];
  PetscTruth     LxSet, LySet, LzSet;
  PetscScalar    my_Lx, my_Ly, my_Lz;

  if (yearsStartRunEndDetermined == PETSC_FALSE) {
    ierr = setStartRunEndYearsFromOptions(PETSC_FALSE);  CHKERRQ(ierr);
  }
  
  // runtime options take precedence in setting of -Lx,-Ly,-Lz *including*
  // if initialization is from an input file
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lx", &my_Lx, &LxSet); CHKERRQ(ierr);
  if (LxSet == PETSC_TRUE) {
    ierr = grid.rescale(my_Lx*1000.0, grid.p->Ly, grid.p->Lz); CHKERRQ(ierr);
  }  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Ly", &my_Ly, &LySet); CHKERRQ(ierr);
  if (LySet == PETSC_TRUE) {
    ierr = grid.rescale(grid.p->Lx, my_Ly*1000.0, grid.p->Lz); CHKERRQ(ierr);
  }  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-Lz", &my_Lz, &LzSet); CHKERRQ(ierr);
  if (LzSet == PETSC_TRUE) {
    ierr = grid.rescale(grid.p->Lx, grid.p->Ly, my_Lz); CHKERRQ(ierr);
  }

  // When we regrid, the `small' model inherits the option `-regrid' but we
  // don't want it to regrid itself, so we short-circuit that action.
  if (allowRegridding == PETSC_TRUE) {
    ierr = PetscOptionsGetString(PETSC_NULL, "-regrid", regridFile, PETSC_MAX_PATH_LEN,
                                 &regridFileSet); CHKERRQ(ierr);
    if (regridFileSet == PETSC_TRUE) {
      ierr = regrid(regridFile); CHKERRQ(ierr);
      ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);  // resolve incompatibilites immediately;
                                                              // this is called anyway at end of first time
                                                              // step
    }
  }

  ierr = stampHistoryCommand(); CHKERRQ(ierr);

  // Note that basal melt rate is not read in as part of a stored model state or
  // an initialization file.  Because basal melt rate is used in computing vertical
  // velocity on the regular grid, and because that calculation could/does precede
  // the first call to temperatureAgeStep(), we want to make sure the first vert
  // velocities do not use junk from uninitialized basal melt rate.
  ierr = VecSet(vbasalMeltRate, 0.0); CHKERRQ(ierr);

  ierr = initBasalTillModel(); CHKERRQ(ierr);
    
  ierr = bedDefSetup(); CHKERRQ(ierr);

  ierr = initPDDFromOptions(); CHKERRQ(ierr);

  ierr = setSoundingFromOptions(); CHKERRQ(ierr);
  ierr = initSounding(); CHKERRQ(ierr);
  tempskipCountDown = 0;

  // initialization should be done; report on starting state
  ierr = verbPrintf(2,grid.com, 
           "  [computational box for ice: (%8.2f km) x (%8.2f km)",
           2*grid.p->Lx/1000.0,2*grid.p->Ly/1000.0); CHKERRQ(ierr);
  if (grid.p->Mbz > 1) {
    ierr = verbPrintf(2,grid.com,
         "\n                                 x (%8.2f m + %7.2f m bedrock)]\n"
         ,grid.p->Lz,grid.p->Lbz); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com," x (%8.2f m)]\n",grid.p->Lz); CHKERRQ(ierr);
  }
  ierr = verbPrintf(2,grid.com, 
           "  [grid cell dimensions     : (%8.2f km) x (%8.2f km) x (%8.2f m)]\n",
           grid.p->dx/1000.0,grid.p->dy/1000.0,grid.p->dz); CHKERRQ(ierr);

  // if -verbose then actually list all of IceParam
  ierr = verbPrintf(3,grid.com,"  IceParam: Mx = %d, My = %d, Mz = %d, Mbz = %d,\n",
                    grid.p->Mx,grid.p->My,grid.p->Mz,grid.p->Mbz); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
           "            Lx = %6.2f km, Ly = %6.2f m, Lz = %6.2f m, Lbz = %6.2f m,\n",
           grid.p->Lx/1000.0,grid.p->Ly/1000.0,grid.p->Lz,grid.p->Lbz); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
           "            dx = %6.3f km, dy = %6.3f km, dz = %6.3f m, year = %8.4f,\n",
           grid.p->dx/1000.0,grid.p->dy/1000.0,grid.p->dz,grid.p->year); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,
     "            history = ****************\n%s            **************************\n"
     ,grid.p->history); CHKERRQ(ierr);
  
  ierr = createViewers(); CHKERRQ(ierr);

  return 0;
}


int IceModel::endOfTimeStepHook() {

  // SIGTERM makes PISM end, saving state under original -o name
  if (pism_signal == SIGTERM) {
    verbPrintf(1, grid.com, "Caught signal SIGTERM: exiting early.\n");
    return 1;
  }
  
  // SIGUSR1 makes PISM save state under filename based on the current year
  if (pism_signal == SIGUSR1) {
    char file_name[PETSC_MAX_PATH_LEN];
    snprintf(file_name, PETSC_MAX_PATH_LEN, "pism-%f.nc", grid.p->year);
    verbPrintf(1, grid.com, "Caught signal SIGUSR1: Writing intermediate file `%s'.\n",
               file_name);
    pism_signal = 0;
    dumpToFile_netCDF(file_name);
  }

  return 0;
}


PetscErrorCode  IceModel::stampHistoryCommand() {
  PetscErrorCode ierr;
  PetscInt argc;
  char **argv;
  
  ierr = PetscGetArgs(&argc, &argv); CHKERRQ(ierr);
  
  char str[HISTORY_STRING_LENGTH];
  
  strncpy(str, argv[0], sizeof(str)); // Does not null terminate on overflow
  str[sizeof(str) - 1] = '\0';
  for (PetscInt i=1; i < argc; i++) {
    PetscInt remaining_bytes = sizeof(str) - strlen(str) - 1;
    // strncat promises to null terminate, so we must only make sure that the
    // end of the buffer is not overwritten.
    strncat(str, " ", remaining_bytes--);
    strncat(str, argv[i], remaining_bytes);
  }

  // All this hooplah is the C equivalent of the Ruby expression
  // ARGV.unshift($0).join(' ') and just ARGV.join(' ') if the executable name
  // was not needed.

  ierr = stampHistory(str); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode  IceModel::stampHistoryEnd() {
  PetscErrorCode ierr;
  PetscLogDouble flops, my_flops;
  char str[HISTORY_STRING_LENGTH];
  MPI_Datatype mpi_type;

  ierr = PetscGetFlops(&my_flops); CHKERRQ(ierr);

  ierr = PetscDataTypeToMPIDataType(PETSC_DOUBLE, &mpi_type); CHKERRQ(ierr);
  MPI_Reduce(&my_flops, &flops, 1, mpi_type, MPI_SUM, 0, grid.com);
  
  snprintf(str, sizeof(str), "Done. Petsc reports %.2f MFlops on %d processes.",
           flops * 1.0e-6, grid.size);

  ierr = stampHistory(str); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode  IceModel::stampHistory(const char* string) {
  PetscErrorCode ierr;

  time_t now;
  tm tm_now;
  now = time(NULL);
  localtime_r(&now, &tm_now);

  char date_str[50];
  // Format specifiers for strftime():
  //   %F : ISO date format
  //   %T : Full 24 hour time
  //   %Z : Time Zone name
  //   %z : Time zone offset
  strftime(date_str, sizeof(date_str), "%F %T %Z", &tm_now);

  char username[50];
  ierr = PetscGetUserName(username, sizeof(username)); CHKERRQ(ierr);
  char hostname[100];
  ierr = PetscGetHostName(hostname, sizeof(hostname)); CHKERRQ(ierr);
  
  char str[HISTORY_STRING_LENGTH];
  int length = snprintf(str, sizeof(str), "%s@%s %s : %s\n",
                        username, hostname, date_str, string);
  
  if (length < 0) {
    SETERRQ(1, "Output error or snprintf() is not C99 compliant.");
    // Old implementations of snprintf() will return `-1' when the string is
    // truncated.  If this is the case on some platform, we need to change this
    // check to allow for that possibility.
  }
  if (length > (int)sizeof(str)) {
    ierr = PetscPrintf(grid.com,
                       "Warning: command line truncated by %d chars in history.\n",
                       length + 1 - sizeof(str)); CHKERRQ(ierr);
    str[sizeof(str) - 2] = '\n';
    str[sizeof(str) - 1] = '\0';
  }

  ierr = stampHistoryString(str); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode  IceModel::stampHistoryString(const char* string) {
  PetscErrorCode ierr;

  PetscInt historyLength = strlen(grid.p->history);
  PetscInt stringLength = strlen(string);

  if (stringLength + historyLength > (int)sizeof(grid.p->history) - 1) {
    ierr = PetscPrintf(grid.com, "Warning: History string overflow.  Truncating history.\n");
    CHKERRQ(ierr);

    // Don't overflow the buffer and null terminate.
    strncpy(grid.p->history, string, sizeof(grid.p->history));
    grid.p->history[sizeof(grid.p->history) - 1] = '\0';
  } else { // We are safe, so we can just write it.
    strcat(grid.p->history, string);
  }
  
  return 0;
}

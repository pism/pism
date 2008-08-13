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

#include <cstring>
#include "iceModel.hh"

#include <petscsys.h>
#include <stdarg.h>
#include <stdlib.h>

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
 
  *gicegridfrac = gicecount / ((PetscScalar) (grid.Mx * grid.My));
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
  const PetscScalar   a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
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
  PetscScalar     **H, **Tbase;
  PetscScalar     melt, divideH, divideT, orig;
  PetscScalar     gmelt, gdivideH, gdivideT, gorig, gvolume, garea;
  PetscScalar     gvolSIA, gvolstream, gvolshelf;
  PetscScalar     meltfrac = 0.0, origfrac = 0.0;
  PetscScalar     *tau;

  ierr = volumeArea(gvolume,garea,gvolSIA, gvolstream, gvolshelf); CHKERRQ(ierr);
  
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  divideH = 0; 
  if (tempAndAge || (getVerbosityLevel() >= 3)) {
    ierr = T3.needAccessToVals(); CHKERRQ(ierr);
    ierr = T3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);
    ierr = T3.doneAccessToVals(); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[0], &Tbase); CHKERRQ(ierr);
    ierr = tau3.needAccessToVals(); CHKERRQ(ierr);
    melt = 0; divideT = 0; orig = 0;
  }
  const PetscScalar   a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
  const PetscScalar   currtime = grid.year * secpera;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
        divideH = H[i][j];
      }
      if (tempAndAge || (getVerbosityLevel() >= 3)) {
        if (H[i][j] > 0) {
          if (useHomoTemp) {
            if (Tbase[i][j] + ice->beta_CC_grad * H[i][j] >= min_temperature_for_SIA_sliding)
              melt += a;
          } else {
            if (Tbase[i][j] >= ice->meltingTemp)
              melt += a;
          }
        }
        if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
          divideT = Tbase[i][j];
        }
        const PetscInt  ks = grid.kBelowHeight(H[i][j]);
        ierr = tau3.getInternalColumn(i,j,&tau); CHKERRQ(ierr);
        for (PetscInt k=1; k<=ks; k++) {
          // ice is original if it is at least one year older than current time
//          if (tau[k] > currtime + secpera)
//            orig += v;
          if (0.5*(tau[k-1]+tau[k]) > currtime + secpera)
            orig += a * 1.0e-3 * (grid.zlevels[k] - grid.zlevels[k-1]);
        }
      }
    }
  }
  
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&divideH, &gdivideH, grid.com); CHKERRQ(ierr);
  if (tempAndAge || (getVerbosityLevel() >= 3)) {
    ierr = DAVecGetArray(grid.da2, vWork2d[0], &Tbase); CHKERRQ(ierr);
    ierr = tau3.doneAccessToVals(); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&melt, &gmelt, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&orig, &gorig, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalMax(&divideT, &gdivideT, grid.com); CHKERRQ(ierr);
    if (gvolume>0) origfrac=gorig/gvolume;
    else origfrac=0.0;
    if (garea>0) meltfrac=gmelt/garea;
    else meltfrac=0.0;
  }

  if (CFLviolcount > 0.5) { // report any CFL violations
    ierr = verbPrintf(2,grid.com,"  [!CFL#=%1.0f (=%8.6f%% of 3D grid)]\n",
              CFLviolcount,100.0 * CFLviolcount/(grid.Mx * grid.Mz * grid.Mz)); CHKERRQ(ierr);
  }
   
  // give summary data a la EISMINT II:
  //    year (+ dt[NR]) (years),
  //       [ note 
  //            N = skipCountDown
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
  ierr = summaryPrintLine(PETSC_FALSE,(PetscTruth)tempAndAge,grid.year,dt,
                          gvolume,garea,meltfrac,gdivideH,gdivideT); CHKERRQ(ierr);
  
  if (getVerbosityLevel() >= 3) {
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
      if ((gmaxu < 0.0) || (gmaxv < 0.0) || (gmaxw < 0.0)) {
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
    const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
    const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0) {

  PetscErrorCode ierr;
  if (printPrototype == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,
      "P         YEAR:     ivol   iarea    meltf     thick0     temp0\n");
    ierr = verbPrintf(2,grid.com,
      "U        years 10^6_km^3 10^6_km^2 (none)          m         K\n");
  } else {
    if (tempAndAge == PETSC_FALSE) {
      ierr = verbPrintf(2,grid.com, "S %12.5f: %8.5f %7.4f   <same> %10.3f    <same>\n",
                         year, volume_kmcube/1.0e6, area_kmsquare/1.0e6, H0); CHKERRQ(ierr);
    } else { // general case
      ierr = verbPrintf(2,grid.com, "S %12.5f: %8.5f %7.4f %8.4f %10.3f %9.4f\n",
                         year, volume_kmcube/1.0e6, area_kmsquare/1.0e6, meltfrac,
                         H0,T0); CHKERRQ(ierr);
    }
  }
  return 0;
}


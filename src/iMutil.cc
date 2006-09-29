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

#include <cstring>
#include <ctime>
#include "iceModel.hh"

#include "petscsys.h"
#include <stdarg.h>
#include <stdlib.h>
#include "petscfix.h"

// verbosity level version of PetscPrintf: only prints if beVerbose == PETSC_TRUE
// modification of PetscPrintf() in src/sys/fileio/mprint.c
// FIXME: want to use an actual integer value for verbosity level
PetscErrorCode IceModel::vPetscPrintf(MPI_Comm comm,const char format[],...)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  size_t         len;
  char           *nformat,*sub1,*sub2;
  PetscReal      value;

  extern FILE *petsc_history;

  PetscFunctionBegin;
  if (!comm) comm = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  if (!rank && ((beVerbose == PETSC_TRUE) || petsc_history) ) {
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
    if (beVerbose == PETSC_TRUE) { // print only if -verbose
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



PetscErrorCode IceModel::computeMaxDiffusivityAndUbar
                      (PetscScalar *gDmax, bool updateDiffusViewer,
                       PetscScalar *gUbarmax, PetscScalar *gUbarSIAav,
                       PetscScalar *gUbarstreamav, PetscScalar *gUbarshelfav,
                       PetscScalar *gicegridfrac, PetscScalar *gSIAgridfrac,
                       PetscScalar *gstreamgridfrac, PetscScalar *gshelfgridfrac) {
  // NOTE:  Assumes IceModel::vubar, vvbar, vu, vv holds correct 
  // and up-to-date values of velocities

  PetscErrorCode ierr;

  PetscScalar **h, **H, **ubar, **vbar, **D, **mask, ***u, ***v;
  PetscScalar Dmax = 0.0, Ubarmax = 0.0, UbarSIAsum = 0.0, Ubarstreamsum = 0.0,
              Ubarshelfsum = 0.0, icecount = 0.0, SIAcount = 0.0, shelfcount = 0.0;

  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &D); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
        icecount += 1.0;
        const PetscScalar Ubarmag = sqrt(PetscSqr(ubar[i][j]) + PetscSqr(vbar[i][j]));
        if (Ubarmag > Ubarmax) Ubarmax = Ubarmag;
        if (intMask(mask[i][j]) == MASK_SHEET) {
          SIAcount += 1.0;
          UbarSIAsum += Ubarmag;
          const PetscScalar h_x=(h[i+1][j]-h[i-1][j])/(2.0*grid.p->dx);
          const PetscScalar h_y=(h[i][j+1]-h[i][j-1])/(2.0*grid.p->dy);
          const PetscScalar alpha = sqrt(PetscSqr(h_x) + PetscSqr(h_y));
          const PetscScalar udef = ubar[i][j] - u[i][j][0];
          const PetscScalar vdef = vbar[i][j] - v[i][j][0];
          const PetscScalar Udefmag = sqrt(PetscSqr(udef) + PetscSqr(vdef));
          const PetscScalar d =
               H[i][j] * Udefmag/(alpha + DEFAULT_ADDED_TO_SLOPE_FOR_DIFF_IN_ADAPTIVE);
          if (d > Dmax) Dmax = d;
          D[i][j] = d;
        } else {
          D[i][j] = 0.0; // no diffusivity in non-SIA regions
          if (modMask(mask[i][j]) == MASK_FLOATING) {
            shelfcount += 1.0;
            Ubarshelfsum += Ubarmag;
          } else if (intMask(mask[i][j]) == MASK_DRAGGING) {
            // streamcount = icecount - SIAcount - shelfcount
            Ubarstreamsum += Ubarmag;
          } else {
            SETERRQ(1,"should not reach here!");
          }
        }
      } else {
        D[i][j] = 0.0;
      }
    }
  }
 
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
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

  ierr = PetscGlobalMax(&Dmax, gDmax, grid.com); CHKERRQ(ierr);
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


PetscErrorCode IceModel::computeMaxDiffusivityONLY
                      (PetscScalar *gDmax, bool updateDiffusViewer) {
  // NOTE:  Assumes IceModel::vubar, vvbar, vu, vv holds correct 
  // and up-to-date values of velocities

  PetscErrorCode ierr;

  PetscScalar **h, **H, **ubar, **vbar, **D, **mask, ***u, ***v;
  PetscScalar Dmax = 0.0;

  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &D); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
        if (intMask(mask[i][j]) == MASK_SHEET) {
          const PetscScalar h_x=(h[i+1][j]-h[i-1][j])/(2.0*grid.p->dx);
          const PetscScalar h_y=(h[i][j+1]-h[i][j-1])/(2.0*grid.p->dy);
          const PetscScalar alpha = sqrt(PetscSqr(h_x) + PetscSqr(h_y));
          const PetscScalar udef = ubar[i][j] - u[i][j][0];
          const PetscScalar vdef = vbar[i][j] - v[i][j][0];
          const PetscScalar Udefmag = sqrt(PetscSqr(udef) + PetscSqr(vdef));
          const PetscScalar d =
               H[i][j] * Udefmag/(alpha + DEFAULT_ADDED_TO_SLOPE_FOR_DIFF_IN_ADAPTIVE);
          if (d > Dmax) Dmax = d;
          D[i][j] = d;
        } else {
          D[i][j] = 0.0; // no diffusivity in non-SIA regions
        }
      } else {
        D[i][j] = 0.0;
      }
    }
  }
 
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
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

  ierr = PetscGlobalMax(&Dmax, gDmax, grid.com); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::adaptTimeStepDiffusivity() {

  PetscErrorCode ierr;
  PetscScalar gDmax;

  ierr = computeMaxDiffusivityONLY(&gDmax, true); CHKERRQ(ierr);
  // note that adapt_ratio * 2 is multiplied by dx^2/(2*maxD) so dt <= adapt_ratio * dx^2/maxD (if dx=dy)
  // reference: Morton & Mayers 2nd ed. pp 62--63
  const PetscScalar gridfactor = 1.0/(grid.p->dx*grid.p->dx) + 1.0/(grid.p->dy*grid.p->dy);
  dt = PetscMin(dt, adaptTimeStepRatio * 2 / ((gDmax + DEFAULT_ADDED_TO_GDMAX_ADAPT) * gridfactor));
  return 0;
}


PetscErrorCode IceModel::adaptTimeStepCFL() {

  // CFLmaxdt is set by computeMaxVelocities in iMvelocity.cc
  dt = PetscMin(dt, CFLmaxdt / tempskip );
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
  if (tempAndAge || (beVerbose == PETSC_TRUE)) {
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
      if (tempAndAge || (beVerbose == PETSC_TRUE)) {
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
  if (tempAndAge || (beVerbose == PETSC_TRUE)) {
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
 
  if (tempAndAge || (beVerbose == PETSC_TRUE)) {
  // give summary data a la EISMINT II:
  //    year (+ dt) (years),
  //    volume (10^6 cubic km),
  //    area (10^6 square km),
  //    basal melt ratio (=(area of base at *absolute*  or *pressure* melting); depends on useHomoTemp),
  //    divide thickness (m),
  //    temp at base at divide (K)  (not homologous),
    ierr = PetscPrintf(grid.com, "%10.2f (+%7.3f):%8.3f%8.3f%9.3f%11.3f%10.3f\n",
                       grid.p->year, dt/secpera, gvolume/1.0e6, garea/1.0e6, meltfrac, gdivideH,
                       gdivideT); CHKERRQ(ierr);
  } else {
  // give REDUCED summary data w/o referencing temp and age fields:
  //    year (+ dt) (years),
  //    volume (10^6 cubic km),
  //    area (10^6 square km),
  //    
  //    divide thickness (m),
  //    
    ierr = PetscPrintf(grid.com, 
       "%10.2f (+%7.3f):%8.3f%8.3f   <same>%11.3f    <same>\n",
       grid.p->year, dt/secpera, gvolume/1.0e6, garea/1.0e6, gdivideH); CHKERRQ(ierr);
  }
  
  if (beVerbose == PETSC_TRUE) {
    // show additional info:
    //   original ice ratio ( = (volume of ice which is original)/(totalvolume) ),
    //   global max of abs(w) in m/a  (if it is not just the default value, which 
    //       means max|w| was not calculated)
    //   the number of CFL violations
    //   the maximum of the diffusivity D on the grid
    PetscScalar Dmax, Ubarmax, UbarSIAav, Ubarstreamav, Ubarshelfav, icegridfrac,
         SIAgridfrac, streamgridfrac, shelfgridfrac;
    ierr = computeMaxDiffusivityAndUbar(&Dmax, false, &Ubarmax,
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
                         Dmax); CHKERRQ(ierr);
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
      if (CFLviolcount > 0.5) {
        ierr = PetscPrintf(grid.com, 
           "  number of temp CFL violations:        %6.0f\n",
           CFLviolcount); CHKERRQ(ierr);
      }
    }
  }
  return 0;
}


PetscErrorCode IceModel::checkForSymmetry(Vec vec, PetscReal *normx, PetscReal *normy, PetscInt stagger) {
  const PetscInt  Mx = grid.p->Mx, My = grid.p->My;
  PetscErrorCode  ierr;
  PetscScalar     **v;
  PetscInt        ind[2];
  Vec va, vb, vc;
  AO  ao;

  ierr = DACreateGlobalVector(grid.da2, &va); CHKERRQ(ierr);
  ierr = VecDuplicate(va, &vb); CHKERRQ(ierr);
  ierr = VecDuplicate(va, &vc); CHKERRQ(ierr);
  ierr = DAGetAO(grid.da2, &ao); CHKERRQ(ierr);
  ierr = DALocalToGlobal(grid.da2, vec, INSERT_VALUES, va); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, va, &v); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (stagger == 0) {
        if (i<Mx-1) ind[0] = (Mx-2-i)*My + j;
        else ind[0] = i*My + j;
      } else {
        ind[0] = (Mx-1-i)*My + j;
      }
      if (stagger == 1) {
        if (j<My-1) ind[1] = i*My + (My-2-j);
        else ind[1] = i*My + j;
      } else {
        ind[1] = i*My + (My-1-j);
      }
      ierr = AOApplicationToPetsc(ao, 2, ind); CHKERRQ(ierr);
      ierr = VecSetValue(vb, ind[0], v[i][j], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(vc, ind[1], v[i][j], INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, va, &v); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(vb); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(vc); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vb); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vc); CHKERRQ(ierr);
  ierr = VecAXPY(vb, (stagger == 0) ? 1 : -1, va); CHKERRQ(ierr);
  ierr = VecAXPY(vc, (stagger == 1) ? 1 : -1, va); CHKERRQ(ierr);
  ierr = VecNorm(vb, NORM_INFINITY, normx); CHKERRQ(ierr);
  ierr = VecNorm(vc, NORM_INFINITY, normy); CHKERRQ(ierr);

  ierr = VecDestroy(va); CHKERRQ(ierr);
  ierr = VecDestroy(vb); CHKERRQ(ierr);
  ierr = VecDestroy(vc); CHKERRQ(ierr);

  return 0;
}


// IceModel::floating(h,H,bed) removed; will use values of vMask for consistency;
// see IceModel::updateSurfaceElevation()


PetscErrorCode IceModel::initFromOptions() {
  PetscErrorCode ierr;
  PetscTruth inFileSet;
  char inFile[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  if (inFileSet == PETSC_TRUE) {
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
  PetscTruth regridFileSet = PETSC_FALSE;
  char regridFile[PETSC_MAX_PATH_LEN];

  // When we regrid, the `small' model inherits the option `-regrid' but we
  // don't want it to regrid itself, so we short-circuit that action.
  if (allowRegridding == PETSC_TRUE) {
    ierr = PetscOptionsGetString(PETSC_NULL, "-regrid", regridFile, PETSC_MAX_PATH_LEN,
                                 &regridFileSet); CHKERRQ(ierr);
    if (regridFileSet == PETSC_TRUE) {
      ierr = regrid(regridFile); CHKERRQ(ierr);
    }
  }

  ierr = stampHistoryCommand(); CHKERRQ(ierr);

  ierr = bedDefSetup(); CHKERRQ(ierr);

  // Note that basal melt rate is not read in as part of a stored model state or
  // an initialization file.  Because basal melt rate is used in computing vertical
  // velocity on the regular grid, and because that calculation could/does precede
  // the first call to temperatureAgeStep(), we want to make sure the first vert
  // velocities do not use junk from uninitialized basal melt rate.
  ierr = VecSet(vbasalMeltRate, 0.0); CHKERRQ(ierr);
  
  if (beVerbose == PETSC_TRUE) {
    ierr = PetscBagView(grid.bag, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(grid.com, 
             "  [computational box for ice: (%8.2f km) x (%8.2f km) x (%8.2f m)]\n",
             2*grid.p->Lx/1000.0,2*grid.p->Ly/1000.0,grid.p->Lz); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, 
             "  [grid cell dimensions     : (%8.2f km) x (%8.2f km) x (%8.2f m)]\n",
             grid.p->dx/1000.0,grid.p->dy/1000.0,grid.p->dz); CHKERRQ(ierr);
  }

  ierr = createViewers(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::regrid(const char *regridFile) {
  PetscErrorCode ierr;
  PetscTruth regridVarsSet;
  char regridVars[PETSC_MAX_PATH_LEN];
  InterpCtx i2, i3, i3b;

  ierr = PetscPrintf(grid.com, "regridding data from %s\n", regridFile); CHKERRQ(ierr);

  IceGrid g(grid.com, grid.rank, grid.size);
  IceModel m(g, ice);
  m.setShowViewers(PETSC_FALSE);
  m.setAllowRegridding(PETSC_FALSE);
  m.initFromFile(regridFile);
  m.afterInitHook();

  if ((grid.p->Lx > m.grid.p->Lx)
      || (grid.p->Ly > m.grid.p->Ly)
      || (grid.p->Lz > m.grid.p->Lz)
      || (grid.p->Lbz > m.grid.p->Lbz)) {
    SETERRQ(1, "New grid exceeds extent of old grid.");
  }

  ierr = getInterpCtx(m.grid.da2, grid.da2, m, i2); CHKERRQ(ierr);
  ierr = vPetscPrintf(grid.com, "  regrid interpolation context i2 created ... "); CHKERRQ(ierr);
  ierr = getInterpCtx(m.grid.da3, grid.da3, m, i3); CHKERRQ(ierr);
  ierr = vPetscPrintf(grid.com, "i3 created ... "); CHKERRQ(ierr);
  ierr = getInterpCtx(m.grid.da3b, grid.da3b, m, i3b); CHKERRQ(ierr);
  ierr = vPetscPrintf(grid.com, "i3b created.\n"); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid_vars", regridVars,
                               PETSC_MAX_PATH_LEN, &regridVarsSet); CHKERRQ(ierr);
  if (regridVarsSet == PETSC_FALSE) {
    // As a default, we only regrid the mask and 3 dimensional quantities.  This
    // is consistent with the standard purpose which is to stick with current
    // geometry through the downscaling procedure.
    strcpy(regridVars, "mTBe");
  }

  ierr = regridVar(regridVars, 'm', i2, m.vMask, vMask);
  ierr = regridVar(regridVars, 'h', i2, m.vh, vh); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'H', i2, m.vH, vH); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'b', i2, m.vbed, vbed); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'a', i2, m.vAccum, vAccum); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 's', i2, m.vTs, vTs); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'g', i2, m.vGhf, vGhf); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'T', i3, m.vT, vT); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'B', i3b, m.vTb, vTb); CHKERRQ(ierr);
  ierr = regridVar(regridVars, 'e', i3, m.vtau, vtau); CHKERRQ(ierr);

  ierr = vPetscPrintf(grid.com, "\n"); CHKERRQ(ierr);

  ierr = destroyInterpCtx(i2); CHKERRQ(ierr);
  ierr = destroyInterpCtx(i3); CHKERRQ(ierr);
  ierr = destroyInterpCtx(i3b); CHKERRQ(ierr);

  PetscScalar run_years = endYear - startYear;
  ierr = setStartYear(m.grid.p->year); CHKERRQ(ierr);
  ierr = setRunYears(run_years); CHKERRQ(ierr);
  grid.p->year = startYear;
  
  // Stamp history
  ierr = stampHistoryString("<---- ---- Begin regrid block ---- ----"); CHKERRQ(ierr);
  ierr = stampHistoryString(m.grid.p->history); CHKERRQ(ierr);
  ierr = stampHistoryString(" ---- ----  End regrid block  ---- ---->\n"); CHKERRQ(ierr);
  
  // The model that we just grabbed data from will fall out of scope and its
  // memory will be freed at the end of this method.  
  return 0;
}


int ind(double i, double j, double k, int N, int K) {
  return (int)i * N * K + (int)j * K + (int)k;
}


PetscErrorCode IceModel::getInterpCtx(const DA dac, const DA daf,
                                      const IceModel &cmodel, InterpCtx &ic) {
  PetscErrorCode ierr;
  PetscInt Mzc, Myc, Mxc, Mzf, Myf, Mxf;
  PetscInt rowl, rowh;
  PetscReal xl, xh, xxl, xxh;
  PetscReal yl, yh, yyl, yyh;
  PetscReal zl, zh, zzl, zzh;
  
  ic.dac = dac;
  ic.daf = daf;
  
  ierr = DAGetInfo(dac, PETSC_NULL, &Mzc, &Myc, &Mxc, PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DAGetInfo(daf, PETSC_NULL, &Mzf, &Myf, &Mxf, PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

  const IceGrid &cgrid = cmodel.grid;
  // The concept of x, y, z are goofy here so that we don't have to duplicate a
  // bunch of code.  We have to make sure that the span (*h - *l) is always positive
  if (Mxf == 1) { // We are working with a 2 dimensional array
    ic.gc = cmodel.g2; ic.gf = g2;
    zl = -grid.p->Ly; zh = grid.p->Ly;
    yl = -grid.p->Lx; yh = grid.p->Lx;
    xl = 0.0; xh = 1.0;
    zzl = -cgrid.p->Ly; zzh = cgrid.p->Ly;
    yyl = -cgrid.p->Lx; yyh = cgrid.p->Lx;
    xxl = 0.0; xxh = 1.0;
  } else { // It is 3 dimensional
    xl = -grid.p->Lx; xh = grid.p->Lx;
    yl = -grid.p->Ly; yh = grid.p->Ly;
    xxl = -cgrid.p->Lx; xxh = cgrid.p->Lx;
    yyl = -cgrid.p->Ly; yyh = cgrid.p->Ly;
    if (Mzf == grid.p->Mz) {
      ic.gc = cmodel.g3; ic.gf = g3;
      zl = 0.0; zh = grid.p->Lz;
      zzl = 0.0; zzh = cgrid.p->Lz;
    } else if (Mzf == grid.p->Mbz || Mzf == 1) {
      ic.gc = cmodel.g3b; ic.gf = g3b;
      // If there is no basal component (Mzf == 1) then this actually doesn't
      // matter since it gets overwritten by the boundary condition in the
      // temperature equation, so we will just do something.
      zl = -grid.p->Lbz; zh = 0.0;
      zzl = -cgrid.p->Lbz; zzh = 0.0;
      if (zl == 0.0) zl = -1.0;
      if (zzl == 0.0) zzl = -1.0;
    } else {
      SETERRQ(1, "Something is broken here.");
    }
  }

  ierr = MatCreateMPIAIJ(grid.com, PETSC_DECIDE, PETSC_DECIDE,
                         Mxf * Myf * Mzf, Mxc * Myc * Mzc,
                         8, PETSC_NULL, 8, PETSC_NULL, &(ic.A)); CHKERRQ(ierr);
  ierr = MatGetVecs(ic.A, &(ic.coarse), &(ic.fine)); CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(ic.A, &rowl, &rowh); CHKERRQ(ierr);

  for (PetscInt row = rowl; row < rowh; row++) {
    const PetscInt i = row / (Myf * Mzf);
    const PetscInt j = (row - i * (Myf * Mzf)) / Mzf;
    const PetscInt k = row % Mzf;

    // We are multiplexing the 2D and 3D cases and here we do not want to divide
    // by zero in the 2D case.
    const PetscReal x = (Mxf == 1) ? 0 : xl + i * (xh - xl) / (Mxf - 1);
    const PetscReal y = yl + j * (yh - yl) / (Myf - 1);
    const PetscReal z = zl + k * (zh - zl) / (Mzf - 1);

    const PetscReal ii = (Mxc - 1) * (x - xxl) / (xxh - xxl);
    const PetscReal jj = (Myc - 1) * (y - yyl) / (yyh - yyl);
    const PetscReal kk = (Mzc - 1) * (z - zzl) / (zzh - zzl);

    // These are in [0,1)
    const PetscReal xx = ii - floor(ii);
    const PetscReal yy = jj - floor(jj);
    const PetscReal zz = kk - floor(kk);
    
    const PetscInt col[] = {
      ind(floor(ii), floor(jj), floor(kk), Myc, Mzc),
      ind(floor(ii), floor(jj), ceil(kk), Myc, Mzc),
      ind(floor(ii), ceil(jj), floor(kk), Myc, Mzc),
      ind(floor(ii), ceil(jj), ceil(kk), Myc, Mzc),
      ind(ceil(ii), floor(jj), floor(kk), Myc, Mzc),
      ind(ceil(ii), floor(jj), ceil(kk), Myc, Mzc),
      ind(ceil(ii), ceil(jj), floor(kk), Myc, Mzc),
      ind(ceil(ii), ceil(jj), ceil(kk), Myc, Mzc)
    };
    const PetscScalar weight[] = {
      (1 - xx) * (1 - yy) * (1 - zz),
      (1 - xx) * (1 - yy) * zz,
      (1 - xx) * yy * (1 - zz),
      (1 - xx) * yy * zz,
      xx * (1 - yy) * (1 - zz),
      xx * (1 - yy) * zz,
      xx * yy * (1 - zz),
      xx * yy * zz
    };

    ierr = MatSetValues(ic.A, 1, &row, 8, col, weight, ADD_VALUES); CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(ic.A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(ic.A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::destroyInterpCtx(InterpCtx &i) {
  PetscErrorCode ierr;

  ierr = MatDestroy(i.A); CHKERRQ(ierr);
  ierr = VecDestroy(i.coarse); CHKERRQ(ierr);
  ierr = VecDestroy(i.fine); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::regridVar(const char *vars, char c, const InterpCtx &ic,
                                   const Vec src, Vec dest) {
  PetscErrorCode ierr;

  if (! strchr(vars, c)) {
    return 0;
  }
  
  // This is a bit ugly.  We need to move the data on local vectors to a
  // global vector which is compatible with the interpolation matrix `ic.A'.
  // Then we move this onto a local vector.  We should be able to do this
  // cleanly with scatters, but my first effort was repelled.  If performance
  // is an issue, then we should improve this part.  Since regridding is only
  // done once per model run, it should not be a problem considering that we
  // do it in parallel.
  ierr = vPetscPrintf(grid.com, "  regridding %c ... ",c); CHKERRQ(ierr);

  ierr = DALocalToGlobal(ic.dac, src, INSERT_VALUES, ic.gc); CHKERRQ(ierr);

  PetscInt low, high, n;
  PetscInt *ida;
  PetscScalar *a;
  PetscInt xs, ys, zs, xm, ym, zm, Mz, My, Mx;
  ierr = DAGetCorners(ic.dac, &zs, &ys, &xs, &zm, &ym, &xm); CHKERRQ(ierr);
  ierr = DAGetInfo(ic.dac, PETSC_NULL, &Mz, &My, &Mx, PETSC_NULL, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

  ida = new PetscInt[xm * ym * zm];
  n = 0;
  for (PetscInt i = xs; i < xs + xm; i++) {
    for (PetscInt j = ys; j < ys + ym; j++) {
      for (PetscInt k = zs; k < zs + zm; k++) {
        ida[n++] = i * (My * Mz) + j * Mz + k;
      }
    }
  }
    
  ierr = VecGetArray(ic.gc, &a); CHKERRQ(ierr);
  ierr = VecSetValues(ic.coarse, n, ida, a, INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecRestoreArray(ic.gc, &a); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ic.coarse); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ic.coarse); CHKERRQ(ierr);
  delete [] ida;
    
  // ierr = VecView(gsrc, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  // ierr = PetscPrintf(grid.com, "#######################################\n"); CHKERRQ(ierr);
  // ierr = VecView(interpCtx.coarse, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = MatMult(ic.A, ic.coarse, ic.fine); CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(ic.fine, &low, &high); CHKERRQ(ierr);
  ierr = VecGetArray(ic.fine, &a); CHKERRQ(ierr);
  ida = new PetscInt[high - low];
  for (PetscInt i = 0; i < high - low; i++) {
    ida[i] = low + i;
  }
  ierr = VecSetValues(ic.gf, high - low, ida, a, INSERT_VALUES);
  ierr = VecAssemblyBegin(ic.gf); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ic.gf); CHKERRQ(ierr);
  delete [] ida;
  ierr = VecRestoreArray(ic.fine, &a); CHKERRQ(ierr);

  ierr = DAGlobalToLocalBegin(ic.daf, ic.gf, INSERT_VALUES, dest); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(ic.daf, ic.gf, INSERT_VALUES, dest); CHKERRQ(ierr);

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
  
  snprintf(str, sizeof(str), "Done. Petsc reports %.2f MFlops.", flops * 1.0e-6);

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


//  PetscErrorCode  IceModel::initCircularFlatBed()  REMOVED!
//  if you want simple circular ice caps to test things use 
//     ./verify -test [BCDFG]   or   ./run_eis -exper [AF]


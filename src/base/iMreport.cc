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

  ierr =    vH.get_array(H);    CHKERRQ(ierr);
  ierr = vubar.get_array(ubar); CHKERRQ(ierr);
  ierr = vvbar.get_array(vbar); CHKERRQ(ierr);
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
        icecount += 1.0;
        const PetscScalar Ubarmag 
                           = sqrt(PetscSqr(ubar[i][j]) + PetscSqr(vbar[i][j]));
        Ubarmax = PetscMax(Ubarmax, Ubarmag);
        if (PismIntMask(mask[i][j]) == MASK_SHEET) {
          SIAcount += 1.0;
          UbarSIAsum += Ubarmag;
        } else if (PismModMask(mask[i][j]) == MASK_FLOATING) {
          shelfcount += 1.0;
          Ubarshelfsum += Ubarmag;
        } else if (PismIntMask(mask[i][j]) == MASK_DRAGGING) {
          // streamcount = icecount - SIAcount - shelfcount
          Ubarstreamsum += Ubarmag;
        } else {
          SETERRQ(1,"should not reach here!");
        }
      }
    }
  }
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr =    vH.end_access(); CHKERRQ(ierr);
  ierr = vubar.end_access(); CHKERRQ(ierr);
  ierr = vvbar.end_access(); CHKERRQ(ierr);

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
  
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  const PetscScalar   a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        area += a;
        const PetscScalar dv = a * H[i][j] * 1e-3;
        volume += dv;
        if (PismIntMask(mask[i][j]) == MASK_SHEET)   volSIA += dv;
        else if (PismIntMask(mask[i][j]) == MASK_DRAGGING)   volstream += dv;
        else if (PismModMask(mask[i][j]) == MASK_FLOATING)   volshelf += dv;
      }
    }
  }  
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volume, &gvolume, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&area, &garea, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volSIA, &gvolSIA, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volstream, &gvolstream, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volshelf, &gvolshelf, grid.com); CHKERRQ(ierr);
  return 0;
}


/*!
Computes fraction of the base which is melted, fraction of the ice which is as 
old as the start of the run (original), and the ice basal temperature at the
center of the ice sheet.

Communication occurs here.
 */
PetscErrorCode IceModel::energyAgeStats(
                    PetscScalar ivol, PetscScalar iarea, bool useHomoTemp, 
                    PetscScalar &gmeltfrac, PetscScalar &gtemp0, PetscScalar &gorigfrac) {
  PetscErrorCode  ierr;
  PetscScalar     **H, **Tbase, *tau;
  PetscScalar     meltarea, temp0, origvol;
  
  // put basal ice temperature in vWork2d[0]
  ierr = T3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);  // z=0 slice
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(Tbase); CHKERRQ(ierr);
  ierr = tau3.begin_access(); CHKERRQ(ierr);

  const PetscScalar   a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
  const PetscScalar   currtime = grid.year * secpera;

  double min_temperature_for_SIA_sliding = config.get("minimum_temperature_for_sliding");
  meltarea = 0.0; temp0 = 0.0; origvol = 0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        // accumulate area of base which is at melt point
        if (useHomoTemp) {
          if (Tbase[i][j] + ice->beta_CC_grad * H[i][j] >= min_temperature_for_SIA_sliding)
            meltarea += a;
        } else {
          if (Tbase[i][j] >= ice->meltingTemp)
            meltarea += a;
        }
        // accumulate volume of ice which is original
        ierr = tau3.getInternalColumn(i,j,&tau); CHKERRQ(ierr);
        const PetscInt  ks = grid.kBelowHeight(H[i][j]);
        for (PetscInt k=1; k<=ks; k++) {
          // ice is original if it is at least one year older than current time
          if (0.5*(tau[k-1]+tau[k]) > currtime + secpera)
            origvol += a * 1.0e-3 * (grid.zlevels[k] - grid.zlevels[k-1]);
        }
      }
      // if you happen to be at center, record basal temp
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
        temp0 = Tbase[i][j];
      }
    }
  }
  
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = tau3.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&meltarea, &gmeltfrac, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&origvol,  &gorigfrac, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&temp0,    &gtemp0,    grid.com); CHKERRQ(ierr);

  // normalize fractions correctly
  if (ivol > 0.0)    gorigfrac = gorigfrac / ivol;
  else gorigfrac = 0.0;
  if (iarea > 0.0)   gmeltfrac = gmeltfrac / iarea;
  else gmeltfrac = 0.0;

  return 0;
}


PetscErrorCode IceModel::summary(bool tempAndAge, bool useHomoTemp) {
  PetscErrorCode  ierr;
  PetscScalar     **H;
  PetscScalar     divideH;
  PetscScalar     gdivideH, gdivideT, gvolume, garea;
  PetscScalar     gvolSIA, gvolstream, gvolshelf;
  PetscScalar     meltfrac = 0.0, origfrac = 0.0;

  ierr = volumeArea(gvolume, garea, gvolSIA, gvolstream, gvolshelf); CHKERRQ(ierr);
  
  // get thick0 = gdivideH
  ierr = vH.get_array(H); CHKERRQ(ierr);
  divideH = 0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
        divideH = H[i][j];
      }
    }
  }  
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&divideH, &gdivideH, grid.com); CHKERRQ(ierr);

  if (tempAndAge || (getVerbosityLevel() >= 3)) {
    ierr = energyAgeStats(gvolume, garea, useHomoTemp, 
                               meltfrac, gdivideT, origfrac); CHKERRQ(ierr);
  }

  // report CFL violations is there are enough
  if (CFLviolcount > 0.0) {
    const PetscScalar CFLviolpercent = 100.0 * CFLviolcount / (grid.Mx * grid.Mz * grid.Mz);
    const PetscScalar CFLVIOL_REPORT_VERB2_PERCENT = 0.1; // only report (verbosity=2) if above 0.1%
    if (CFLviolpercent > CFLVIOL_REPORT_VERB2_PERCENT) {
      ierr = verbPrintf(2,grid.com,"  [!CFL#=%1.0f (=%8.6f%% of 3D grid)]\n",
                CFLviolcount,CFLviolpercent); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(3,grid.com,"  [!CFL#=%1.0f (=%8.6f%% of 3D grid)]\n",
                CFLviolcount,CFLviolpercent); CHKERRQ(ierr);
    }
  }
   
  // main report: 'S' line
  ierr = summaryPrintLine(PETSC_FALSE,(PetscTruth)tempAndAge,grid.year,dt,
                          gvolume,garea,meltfrac,gdivideH,gdivideT); CHKERRQ(ierr);

  // extra verbose report  
  const PetscScalar EXTRAS_VERB_LEVEL = 4;
  if (getVerbosityLevel() >= EXTRAS_VERB_LEVEL) {
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


//! Print a line to stdout which summarizes the state of the modeled ice sheet at the end of the time step.
/*!
Generally, a single line is printed to stdout, starting with the character 'S' in the 
left-most column.

If IceModel::printPrototype is TRUE then alternate lines with
different left-most characters are printed:
  - 'P' line gives names of the quantities reported in the 'S' line, the "prototype", while
  - 'U' line gives units of these quantities.

The left-most character convention allows automatic tools to read PISM stdout
and produce time-series.  The 'P' and 'U' lines are intended to appear once at the
beginning of the run, while an 'S' line appears at every time step.  Some 'S'
lines report "<same>" when there is no change to a quantity.

This base class version gives a report based on the information included in the 
EISMINT II intercomparison of ice sheet models[\ref EISMINT00].  The 'P' and 'U' lines are \code
  P         YEAR:     ivol   iarea    meltf     thick0     temp0
  U        years 10^6_km^3 10^6_km^2 (none)          m         K
\endcode
The 'S' line gives the corresponding numbers.

Note that \c ivol is the ice sheet volume and \c iarea is the area occupied 
by positive thickness ice.  The pure number \c meltf is the fraction of \c iarea for which
the base is at the melting temperature.  (xactly what this melting temperature refers to
is determined in IceModel::summary().)  The final two quantities are
simply values of two variable at the center of the computational domain, namely 
the thickness (\c thick0) and basal ice absolute temperature (\c temp0) at that
center point.

For more description and examples, see the PISM User's Manual.
  
Derived classes of IceModel are encouraged to redefine this method and add 
alternate information.
 */
PetscErrorCode IceModel::summaryPrintLine(
     PetscTruth printPrototype,  bool tempAndAge,
     PetscScalar year,  PetscScalar /* delta_t */,
     PetscScalar volume_kmcube,  PetscScalar area_kmsquare,
     PetscScalar meltfrac,  PetscScalar H0,  PetscScalar T0) {

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

//! \brief Computes cbar, the magnitude of vertically-integrated horizontal
//! velocity of ice and masks out ice-free areas.
PetscErrorCode IceModel::compute_cbar(IceModelVec2 &result) {
  PetscErrorCode ierr;

  ierr = result.set_to_magnitude(vubar, vvbar); CHKERRQ(ierr);
  ierr = result.mask_by(vH); CHKERRQ(ierr); // mask out ice-free areas

  ierr = result.set_name("cbar"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", 
			  "magnitude of vertically-integrated horizontal velocity of ice",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes cflx, the magnitude of vertically-integrated horizontal
//! flux of ice.
PetscErrorCode IceModel::compute_cflx(IceModelVec2 &result, IceModelVec2 &cbar) {
  PetscErrorCode ierr;

  ierr = cbar.multiply_by(vH, result); CHKERRQ(ierr);

  ierr = result.set_name("cflx"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", 
			  "magnitude of vertically-integrated horizontal flux of ice",
			  "m2 s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m2 year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes cbase, the magnitude of horizontal velocity of ice at base
//! of ice and masks out ice-free areas.
//! Uses \c tmp as a preallocated temporary storage.
PetscErrorCode IceModel::compute_cbase(IceModelVec2 &result, IceModelVec2 &tmp) {
  PetscErrorCode ierr;

  ierr = u3.getHorSlice(result, 0.0); CHKERRQ(ierr); // result = u_{z=0}
  ierr = v3.getHorSlice(tmp, 0.0); CHKERRQ(ierr);    // tmp = v_{z=0}

  ierr = result.set_to_magnitude(result,tmp); CHKERRQ(ierr);

  ierr = result.set_name("cbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", 
			  "magnitude of horizontal velocity of ice at base of ice",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes csurf, the magnitude of horizontal velocity of ice at ice
//! surface and masks out ice-free areas. Uses \c tmp as a preallocated
//! temporary storage.
PetscErrorCode IceModel::compute_csurf(IceModelVec2 &result, IceModelVec2 &tmp) {
  PetscErrorCode ierr;

  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getSurfaceValues(result, vH); CHKERRQ(ierr);
  ierr = v3.getSurfaceValues(tmp,    vH); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);

  ierr = result.set_to_magnitude(result, tmp); CHKERRQ(ierr);

  ierr = result.set_name("csurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", 
			  "magnitude of horizontal velocity of ice at ice surface",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes uvelsurf, the x component of velocity of ice at ice surface.
PetscErrorCode IceModel::compute_uvelsurf(IceModelVec2 &result) {
  PetscErrorCode ierr;

  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getSurfaceValues(result, vH); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("uvelsurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "x component of velocity of ice at ice surface",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", -1e6); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes vvelsurf, the y component of velocity of ice at ice surface.
PetscErrorCode IceModel::compute_vvelsurf(IceModelVec2 &result) {
  PetscErrorCode ierr;

  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = v3.getSurfaceValues(result, vH); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("vvelsurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "y component of velocity of ice at ice surface",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", -1e6); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes wvelsurf, the vertical velocity of ice at ice surface.
/*! Note that there is no need to mask out ice-free areas here, because
  wvelsurf is zero at those locations.
 */
PetscErrorCode IceModel::compute_wvelsurf(IceModelVec2 &result) {
  PetscErrorCode ierr;
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = w3.getSurfaceValues(result, vH); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("wvelsurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "vertical velocity of ice at ice surface",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", -1e6); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes taud, the magnitude of driving shear stress at base of ice.
//! Uses tmp as a preallocated temporary storage.
PetscErrorCode IceModel::compute_taud(IceModelVec2 &result, IceModelVec2 &tmp) {
  PetscErrorCode ierr;

  ierr = computeDrivingStress(result, tmp); CHKERRQ(ierr);

  ierr = result.set_to_magnitude(result, tmp); CHKERRQ(ierr);
  ierr = result.set_name("taud"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic",
			  "magnitude of driving shear stress at base of ice",
			  "Pa", ""); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  return 0;
}

//! \brief Sets entrues of result to corresponding processor ranks.
PetscErrorCode IceModel::compute_rank(IceModelVec2 &result) {
  PetscErrorCode ierr;
  PetscScalar **a;

  ierr = result.get_array(a); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j)
      a[i][j] = grid.rank;
  ierr = result.end_access();

  ierr = result.set_name("rank"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "processor rank", "", ""); CHKERRQ(ierr);
  result.time_independent = true;
  return 0;
}

//! Compute the pressure-adjusted temperature in degrees C corresponding to T3, and put in a global IceModelVec3 provided by user.
/*!
This procedure is put here in IceModel to facilitate comparison of IceModel and IceEnthalpyModel
results.  It is called by giving option -temp_pa.
 */
PetscErrorCode IceModel::compute_temp_pa(IceModelVec3 &useForPATemp) {
  PetscErrorCode ierr;

  PetscScalar **thickness;
  PetscScalar *Tpaij, *Tij; // columns of these values
  ierr = useForPATemp.begin_access(); CHKERRQ(ierr);
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = vH.get_array(thickness); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = useForPATemp.getInternalColumn(i,j,&Tpaij); CHKERRQ(ierr);
      ierr = T3.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        Tpaij[k] = Tij[k] - ice->meltingTemp;  // un-adjusted, but in deg_C
        const PetscScalar depth = thickness[i][j] - grid.zlevels[k];
        if (depth > 0.0)  Tpaij[k] += ice->beta_CC_grad * depth;
      }
    }
  }
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = useForPATemp.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = useForPATemp.set_name("temp_pa"); CHKERRQ(ierr);
  ierr = useForPATemp.set_attrs("diagnostic",
       "pressure-adjusted ice temperature (degrees C)", "", ""); CHKERRQ(ierr);

  // communication not done; we allow global IceModelVec3s as useForPATemp
  return 0;
}

//! Computes the vertically-averaged ice hardness.
PetscErrorCode IceModel::compute_hardav(IceModelVec2 &result) {
  PetscErrorCode ierr;
  
  const PetscScalar fillval = -1.0;
  
  PetscScalar **hardav, **thk;
  PetscScalar *Tij; // columns of temperature values
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = vH.get_array(thk); CHKERRQ(ierr);
  ierr = result.get_array(hardav); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = T3.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      const PetscScalar H = thk[i][j];
      if (H > 0.0) {
        hardav[i][j] = ice->averagedHarness(H, grid.kBelowHeight(H), grid.zlevels, Tij);
      } else { // put negative value below valid range
        hardav[i][j] = fillval;
      }
    }
  }
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("hardav"); CHKERRQ(ierr);
  const PetscScalar power = 1.0 / ice->exponent();
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);
  ierr = result.set_attrs("diagnostic", "vertical average of ice hardness",
			  unitstr, ""); CHKERRQ(ierr);

  ierr = result.set_attr("valid_min",0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue",fillval); CHKERRQ(ierr);
  return 0;
}

//! Computes the basal water pressure
/*!
  \f[p_w = \alpha\, \frac{w}{w_{\text{max}}}\, \rho\, g\, H,\f]
  where 

  - \f$\alpha\f$ is the till pore water fraction (till_pw_fraction),
  - \f$w\f$ is the effective thickness of subglacial melt water (bwat)
  - \f$w_{\text{max}}\f$ is the maximum allowed value for \f$w\f$ (max_hmelt),
  - \f$\rho\f$ is the ice density (ice_density)
  - \f$H\f$ is the ice thickness (thk)
 */
PetscErrorCode IceModel::compute_bwp(IceModelVec2 &result) {
  PetscErrorCode ierr;
  PetscScalar
    alpha     = config.get("till_pw_fraction"),
    wmax      = config.get("max_hmelt"),
    g         = config.get("standard_gravity"),
    rho       = config.get("ice_density");

  ierr = vH.multiply_by(vHmelt, result); CHKERRQ(ierr);
  ierr = result.scale((alpha / wmax) * rho * g); CHKERRQ(ierr);

  ierr = result.set_name("bwp"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "basal pore water pressure",
			  "Pa", ""); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes the rate of change of ice surface elevation as a sum of the
//! bedrock uplift rate and the thickness rate of change.
PetscErrorCode IceModel::compute_dhdt(IceModelVec2 &result) {
  PetscErrorCode ierr;

  ierr = vuplift.add(1.0, vdHdt, result); CHKERRQ(ierr);

  ierr = result.set_name("dhdt"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "rate of change of surface elevation",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1");
  result.write_in_glaciological_units = true;

  return 0;
}

PetscErrorCode IceModel::compute_uvelbase(IceModelVec2 &result) {
  PetscErrorCode ierr;

  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getHorSlice(result, 0.0); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("uvelbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "x component of ice velocity at the base of ice",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", -1e6); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::compute_vvelbase(IceModelVec2 &result) {
  PetscErrorCode ierr;
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = v3.getHorSlice(result, 0.0); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("vvelbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "y component of ice velocity at the base of ice",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", -1e6); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::compute_wvelbase(IceModelVec2 &result) {
  PetscErrorCode ierr;
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = w3.getHorSlice(result, 0.0); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("wvelbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "vertical velocity of ice at the base of ice",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", -1e6); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes ice temperature at the base of ice.
PetscErrorCode IceModel::compute_tempbase(IceModelVec2 &result) {
  PetscErrorCode ierr;

  // put basal ice temperature in vWork2d[0]
  ierr = T3.getHorSlice(result, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = result.set_name("tempbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "ice temperature at the base of ice",
			  "Kelvin", ""); CHKERRQ(ierr);

  PetscScalar fill_value = GSL_NAN;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes ice temperature at the 1 m below the surface.
PetscErrorCode IceModel::compute_tempsurf(IceModelVec2 &result) {
  PetscErrorCode ierr;
  PetscScalar fill_value = GSL_NAN;

  // compute levels corresponding to 1 m below the ice surface:

  PetscScalar **thk, **res;
  ierr = vH.get_array(thk); CHKERRQ(ierr);
  ierr = result.get_array(res); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	res[i][j] = PetscMax(thk[i][j] - 1.0, 0.0);
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = T3.getSurfaceValues(result, result); CHKERRQ(ierr);  // z=0 slice

  ierr = vH.get_array(thk); CHKERRQ(ierr);
  ierr = result.get_array(res); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (thk[i][j] <= 1.0)
	res[i][j] = fill_value;
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("tempsurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "ice temperature at 1m below the ice surface",
			  "Kelvin", ""); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);
  return 0;
}


//! \brief Computes a diagnostic quantity given by \c name and returns a
//! pointer to a pre-allocated work vector containing it.
/*! For 2D quantities, result will point to vWork2d[0].

  For 3D -- to T3new, because we don't have a general-purpose 3D work vector.

  Note that (depending on the quantity requested) vWork2d[1] might get used as
  a temporary storage.
 */
PetscErrorCode IceModel::compute_by_name(string name, IceModelVec* &result) {
  PetscErrorCode ierr;

  result = NULL;		// if clauses can override this

  if (name == "bwp") {
    ierr = compute_bwp(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "cbar") {
    ierr = compute_cbar(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "cbase") {
    ierr = compute_cbase(vWork2d[0], vWork2d[1]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "cflx") {
    ierr = compute_cbar(vWork2d[1]); CHKERRQ(ierr);
    ierr = compute_cflx(vWork2d[0], vWork2d[1]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "csurf") {
    ierr = compute_csurf(vWork2d[0], vWork2d[1]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "dhdt") {
    ierr = compute_dhdt(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "hardav") {
    ierr = compute_hardav(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "taud") {
    ierr = compute_taud(vWork2d[0], vWork2d[1]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "temp_pa") {
    ierr = compute_temp_pa(Tnew3); CHKERRQ(ierr);
    result = &Tnew3;
    return 0;
  }

  if (name == "tempsurf") {
    ierr = compute_tempsurf(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "tempbase") {
    ierr = compute_tempbase(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "uvelbase") {
    ierr = compute_uvelbase(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "vvelbase") {
    ierr = compute_vvelbase(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "wvelbase") {
    ierr = compute_wvelbase(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "uvelsurf") {
    ierr = compute_uvelsurf(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "vvelsurf") {
    ierr = compute_vvelsurf(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "wvelsurf") {
    ierr = compute_wvelsurf(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "rank") {
    ierr = compute_rank(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  return 0;
}

//! Computes the ice volume, in m^3.
PetscErrorCode IceModel::compute_ice_volume(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     **H;
  PetscScalar     volume=0.0;
  const PetscScalar a = grid.dx * grid.dy; // cell area
  
  ierr = vH.get_array(H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        volume += a * H[i][j];
      }
    }
  }  
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     **H;
  PetscScalar     area=0.0;
  const PetscScalar a = grid.dx * grid.dy; // cell area
  
  ierr = vH.get_array(H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0)
        area += a;
    }
  }  
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes grounded ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_grounded(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     **H, **mask;
  PetscScalar     area=0.0;
  const PetscScalar a = grid.dx * grid.dy; // cell area
  
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ((H[i][j] > 0) &&
	  ((PismIntMask(mask[i][j]) == MASK_SHEET) || (PismIntMask(mask[i][j]) == MASK_DRAGGING)))
        area += a;
    }
  }  
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes floating ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_floating(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     **H, **mask;
  PetscScalar     area=0.0;
  const PetscScalar a = grid.dx * grid.dy; // cell area
  
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ((H[i][j] > 0) && (PismModMask(mask[i][j]) == MASK_FLOATING))
        area += a;
    }
  }  
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Compute a scalar diagnostic quantity by name.
PetscErrorCode IceModel::compute_by_name(string name, PetscScalar &result) {
  PetscErrorCode ierr, errcode = 1;

  if (name == "ivol") {
    errcode = 0;
    ierr = compute_ice_volume(result); CHKERRQ(ierr);
  }

  if (name == "imass") {
    errcode = 0;
    PetscScalar ice_density = config.get("ice_density");
    ierr = compute_ice_volume(result); CHKERRQ(ierr);
    result *= ice_density;
  }

  if (name == "iarea") {
    errcode = 0;
    ierr = compute_ice_area(result); CHKERRQ(ierr);
  }

  if (name == "iareag") {
    errcode = 0;
    ierr = compute_ice_area_grounded(result); CHKERRQ(ierr);
  }

  if (name == "iareaf") {
    errcode = 0;
    ierr = compute_ice_area_floating(result); CHKERRQ(ierr);
  }

  if (name == "dt") {
    errcode = 0;
    result = dt;
  }

  if (name == "divoldt") {
    errcode = 0;
    result = dvoldt;
  }

  if (name == "dimassdt") {
    errcode = 0;
    PetscScalar ice_density = config.get("ice_density");
    result = dvoldt * ice_density;
  }

  if (name == "total_basal_ice_flux") {
    errcode = 0;
    result = total_basal_ice_flux;
  }

  if (name == "total_surface_ice_flux") {
    errcode = 0;
    result = total_surface_ice_flux;
  }

  if (name == "total_sub_shelf_ice_flux") {
    errcode = 0;
    result = total_sub_shelf_ice_flux;
  }

  return errcode;
}

/*! Computes total ice fluxes (in kg per second) at 3 interfaces:

  \li "ice-atmoshpere" interface using the accumulation/ablation rate
  \li "ice-bedrock" interface using the basal melt rate
  \li bottom of ice shelves using the sub-ice-shelf melt rate

  These quantities should be understood as <i>instantaneous at the beginning of
  the time-step.</i> Multiplying by dt will \b not necessarily give the
  corresponding change from the beginning to the end of the time-step.
 */
PetscErrorCode IceModel::ice_mass_bookkeeping() {
  PetscErrorCode ierr;
  IceModelVec2 *surface_mf, *subshelf_mf;
  PetscScalar **thk, **acab, **bmr_grounded, **bmr_floating, **mask;
  PetscReal cell_area = grid.dx * grid.dy;

  // call sets pccsmf to point to IceModelVec2 with current surface mass flux
  if (atmosPCC != PETSC_NULL) {
    ierr = atmosPCC->updateSurfMassFluxAndProvide(grid.year, dt / secpera,
						  surface_mf); CHKERRQ(ierr);
  } else { SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL"); }
  // call sets pccsbmf to point to IceModelVec2 with current mass flux under shelf base
  if (oceanPCC != PETSC_NULL) {
    ierr = oceanPCC->updateShelfBaseMassFluxAndProvide(grid.year, dt / secpera,
						       subshelf_mf); CHKERRQ(ierr);
  } else { SETERRQ(2,"PISM ERROR: oceanPCC == PETSC_NULL"); }

  total_basal_ice_flux = 0;
  total_surface_ice_flux = 0;
  total_sub_shelf_ice_flux = 0;

  ierr = surface_mf->get_array(acab); CHKERRQ(ierr);
  ierr = subshelf_mf->get_array(bmr_floating); CHKERRQ(ierr);
  ierr = vbasalMeltRate.get_array(bmr_grounded); CHKERRQ(ierr);
  ierr = vH.get_array(thk); CHKERRQ(ierr);
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // ignore ice-free cells:
      if (thk[i][j] <= 0.0)
	continue;

      // add the accumulation/ablation rate:
      total_surface_ice_flux += acab[i][j]; // note the "+="!

      // add the sub-shelf melt rate;
      if (PismIntMask(mask[i][j]) == MASK_FLOATING) {
	total_sub_shelf_ice_flux -= bmr_floating[i][j]; // note the "-="!
      }

      // add the basal melt rate:
      if ( (PismIntMask(mask[i][j]) == MASK_SHEET) ||
	   (PismIntMask(mask[i][j]) == MASK_DRAGGING) ) {
	total_basal_ice_flux -= bmr_grounded[i][j]; // note the "-="!
      }

    }
  }  
  ierr = surface_mf->end_access(); CHKERRQ(ierr);
  ierr = subshelf_mf->end_access(); CHKERRQ(ierr);
  ierr = vbasalMeltRate.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  PetscScalar ice_density = config.get("ice_density");

  total_surface_ice_flux     *= ice_density * cell_area;
  total_sub_shelf_ice_flux   *= ice_density * cell_area;
  total_basal_ice_flux       *= ice_density * cell_area;

  ierr = PetscGlobalSum(&total_surface_ice_flux,   &total_surface_ice_flux,   grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&total_sub_shelf_ice_flux, &total_sub_shelf_ice_flux, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&total_basal_ice_flux,     &total_basal_ice_flux,     grid.com); CHKERRQ(ierr);

  return 0;
}

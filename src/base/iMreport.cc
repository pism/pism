// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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
  // NOTE:  Assumes IceModel::vel_bar, vu, vv holds correct 
  // and up-to-date values of velocities

  PetscErrorCode ierr;

  PetscScalar **H;
  PetscScalar Ubarmax = 0.0, UbarSIAsum = 0.0, Ubarstreamsum = 0.0,
              Ubarshelfsum = 0.0, icecount = 0.0, SIAcount = 0.0, shelfcount = 0.0;

  ierr =    vH.get_array(H);    CHKERRQ(ierr);
  ierr = vel_bar.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
        icecount += 1.0;
        const PetscScalar Ubarmag 
                           = sqrt(PetscSqr(vel_bar(i,j).u) + PetscSqr(vel_bar(i,j).v));
        Ubarmax = PetscMax(Ubarmax, Ubarmag);
        if (vMask.value(i,j) == MASK_SHEET) {
          SIAcount += 1.0;
          UbarSIAsum += Ubarmag;
        } else if (vMask.is_floating(i,j)) {
          shelfcount += 1.0;
          Ubarshelfsum += Ubarmag;
        } else if (vMask.value(i,j) == MASK_DRAGGING_SHEET) {
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
  ierr = vel_bar.end_access(); CHKERRQ(ierr);

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


//!  Computes volume and area of ice sheet, for reporting purposes.
/*!
Also computes ice volumes in the three mask regions SHEET, DRAGGING, FLOATING.

Communication done for global max and global sum.

Returns area in units of m^2 and volume in m^3.
 */
PetscErrorCode IceModel::volumeArea(PetscScalar& gvolume, PetscScalar& garea,
                                    PetscScalar& gvolSIA, PetscScalar& gvolstream, 
                                    PetscScalar& gvolshelf) {

  PetscErrorCode  ierr;
  PetscScalar     **H;
  PetscScalar     volume=0.0, area=0.0, volSIA=0.0, volstream=0.0, volshelf=0.0;
  
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  const PetscScalar   a = grid.dx * grid.dy; // area unit (m^2)
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        area += a;
        const PetscScalar dv = a * H[i][j];
        volume += dv;
        if (vMask.value(i,j) == MASK_SHEET)   volSIA += dv;
        else if (vMask.value(i,j) == MASK_DRAGGING_SHEET)   volstream += dv;
        else if (vMask.is_floating(i,j))   volshelf += dv;
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
Computes fraction of the base which is melted and the ice basal temperature at the
center of the ice sheet.

Communication occurs here.
 */
PetscErrorCode IceModel::energyStats(PetscScalar iarea, PetscScalar &gmeltfrac,
                                     PetscScalar &gtemp0) {
  PetscErrorCode    ierr;
  PetscScalar       meltarea = 0.0, temp0 = 0.0;
  const PetscScalar a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
  
  ierr = vH.begin_access(); CHKERRQ(ierr);

  // use Enth3 to get stats
  ierr = Enth3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);  // z=0 slice
  PetscScalar **Enthbase;
  ierr = vWork2d[0].get_array(Enthbase); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
	// accumulate area of base which is at melt point
	if (EC->isTemperate(Enthbase[i][j], EC->getPressureFromDepth(vH(i,j)) ))  
	  meltarea += a;
      }
      // if you happen to be at center, record absolute basal temp there
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
	ierr = EC->getAbsTemp(Enthbase[i][j],EC->getPressureFromDepth(vH(i,j)), temp0);
	CHKERRQ(ierr);
      }
    }
  }
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);

  ierr = vH.end_access(); CHKERRQ(ierr);

  // communication
  ierr = PetscGlobalSum(&meltarea, &gmeltfrac, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&temp0,    &gtemp0,    grid.com); CHKERRQ(ierr);

  // normalize fraction correctly
  if (iarea > 0.0)   gmeltfrac = gmeltfrac / iarea;
  else gmeltfrac = 0.0;

  return 0;
}


/*!
Computes fraction of the ice which is as old as the start of the run (original).
Communication occurs here.
 */
PetscErrorCode IceModel::ageStats(PetscScalar ivol, PetscScalar &gorigfrac) {
  PetscErrorCode  ierr;

  gorigfrac = -1.0;  // result value if not do_age

  if (!config.get_flag("do_age")) 
    return 0;  // leave now

  const PetscScalar  a = grid.dx * grid.dy * 1e-3 * 1e-3, // area unit (km^2)
                     currtime = grid.year * secpera; // seconds

  PetscScalar **H, *tau, origvol = 0.0;
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = tau3.begin_access(); CHKERRQ(ierr);

  // compute local original volume
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        // accumulate volume of ice which is original
        ierr = tau3.getInternalColumn(i,j,&tau); CHKERRQ(ierr);
        const PetscInt  ks = grid.kBelowHeight(H[i][j]);
        for (PetscInt k=1; k<=ks; k++) {
          // ice in segment is original if it is as old as one year less than current time
          if (0.5*(tau[k-1]+tau[k]) > currtime - secpera)
            origvol += a * 1.0e-3 * (grid.zlevels[k] - grid.zlevels[k-1]);
        }
      }
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = tau3.end_access(); CHKERRQ(ierr);

  // communicate to turn into global original fraction
  ierr = PetscGlobalSum(&origvol,  &gorigfrac, grid.com); CHKERRQ(ierr);

  // normalize fraction correctly
  if (ivol > 0.0)    gorigfrac = gorigfrac / ivol;
  else gorigfrac = 0.0;

  return 0;
}


PetscErrorCode IceModel::summary(bool tempAndAge) {
  PetscErrorCode  ierr;
  PetscScalar     **H;
  PetscScalar     divideH;
  PetscScalar     gdivideH, gdivideT, gvolume, garea;
  PetscScalar     gvolSIA, gvolstream, gvolshelf;
  PetscScalar     meltfrac = 0.0, origfrac = 0.0;

  // get volumes in m^3 and areas in m^2
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
    ierr = energyStats(garea, meltfrac, gdivideT); CHKERRQ(ierr);
  }

  if ((tempAndAge || (getVerbosityLevel() >= 3)) && (config.get_flag("do_age"))) {
    ierr = ageStats(gvolume, origfrac); CHKERRQ(ierr);
  }

  // report CFL violations
  if (CFLviolcount > 0.0) {
    const PetscScalar CFLviolpercent = 100.0 * CFLviolcount / (grid.Mx * grid.Mz * grid.Mz);
    // at default verbosity level, only report CFL viols if above:
    const PetscScalar CFLVIOL_REPORT_VERB2_PERCENT = 0.1;
    if (   (CFLviolpercent > CFLVIOL_REPORT_VERB2_PERCENT) 
        || (getVerbosityLevel() > 2) ) {
      char tempstr[90] = "";
      snprintf(tempstr,90,
              "  [!CFL#=%1.0f (=%5.2f%% of 3D grid)] ",
              CFLviolcount,CFLviolpercent);
      stdout_flags = tempstr + stdout_flags;
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
    ierr = computeFlowUbarStats(&Ubarmax,
              &UbarSIAav, &Ubarstreamav, &Ubarshelfav, &icegridfrac,
              &SIAgridfrac, &streamgridfrac, &shelfgridfrac); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, 
           "    (volume of ice which is SIA, stream, shelf:  %8.3f,  %8.3f,  %8.3f)\n",
           gvolSIA/(1.0e6*1.0e9), gvolstream/(1.0e6*1.0e9), gvolshelf/(1.0e6*1.0e9));
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
      if (config.get_flag("do_age")) {
        ierr = PetscPrintf(grid.com, 
           "  fraction of ice which is original: %9.3f\n",
                         origfrac); CHKERRQ(ierr);
      }
    }
  }
  return 0;
}


//! Print a line to stdout which summarizes the state of the modeled ice sheet at the end of the time step.
/*!
Generally, two lines are printed to stdout, the first starting with a space 
and the second starting with the character 'S' in the left-most column (column 1).

The first line shows flags for which processes executed, and the length of the
time step (and/or substeps under option -skip).  See IceModel::run()
for meaning of these flags.

If IceModel::printPrototype is TRUE then the first line does not appear and
the second line has alternate appearance.  Specifically, different column 1
characters are printed:
  - 'P' line gives names of the quantities reported in the 'S' line, the
    "prototype", while
  - 'U' line gives units of these quantities.
This column 1 convention allows automatic tools to read PISM stdout
and produce time-series.  The 'P' and 'U' lines are intended to appear once at
the beginning of the run, while an 'S' line appears at every time step.
This base class version gives a report based on the information included in the 
EISMINT II intercomparison of ice sheet models[\ref EISMINT00].

Note that the inputs \c volume and \c area to this method are in m^3 and m^2,
respectively.  Thus all inputs to this method are in MKS except for \c year.

The resulting numbers on an 'S' line have the following meaning in this base
class version:
  - \c ivol is the total ice sheet volume
  - \c iarea is the total area occupied by positive thickness ice
  - \c thick0 is the ice thickness at the center of the computational domain
  - \c temp0 is the ice basal temperature at the center of the computational domain
The last two can be interpreted as "sanity checks", because they give
information about a location which may or may not be "typical".

For more description and examples, see the PISM User's Manual.
Derived classes of IceModel may redefine this method and print alternate
information.  Use of DiagnosticTimeseries may be superior, however.
 */
PetscErrorCode IceModel::summaryPrintLine(
     PetscTruth printPrototype,  bool tempAndAge,
     PetscScalar year,  PetscScalar delta_t,
     PetscScalar volume,  PetscScalar area,
     PetscScalar /* meltfrac */,  PetscScalar H0,  PetscScalar T0) {

  PetscErrorCode ierr;
  const bool do_temp = config.get_flag("do_temp");

  const int log10scale = static_cast<int>(config.get("summary_volarea_scale_factor_log10"));
  const double scale = pow(10.0, static_cast<double>(log10scale));
  char  volscalestr[10] = "     ", // for special case: blank when 10^0 = 1 scaling
        areascalestr[10] = "   ";  // ditto
  if (log10scale != 0) {
    snprintf(volscalestr, sizeof(volscalestr), "10^%1d_", log10scale);
    strcpy(areascalestr,volscalestr);
  }

  // this version keeps track of what has been done so as to minimize stdout:
  static string stdout_flags_count0;
  static int    mass_cont_sub_counter = 0;  
  static double mass_cont_sub_dtsum = 0.0;
  
  if (printPrototype == PETSC_TRUE) {
    if (do_temp) {
      ierr = verbPrintf(2,grid.com,
          "P         YEAR:     ivol   iarea     thick0     temp0\n");
      ierr = verbPrintf(2,grid.com,
          "U        years %skm^3 %skm^2        m         K\n",
          volscalestr,areascalestr);
    } else {
      ierr = verbPrintf(2,grid.com,
          "P         YEAR:     ivol   iarea     thick0\n");
      ierr = verbPrintf(2,grid.com,
          "U        years %skm^3 %skm^2        m\n",
          volscalestr,areascalestr);
    }
  } else {
    if (mass_cont_sub_counter == 0)
      stdout_flags_count0 = stdout_flags;
    if (delta_t > 0.0) {
      mass_cont_sub_counter++;      
      mass_cont_sub_dtsum += delta_t;
    }
    if ((tempAndAge == PETSC_TRUE) || (!do_temp) || (getVerbosityLevel() > 2)) {
      char tempstr[90] = "";
      const PetscScalar major_dt_years = mass_cont_sub_dtsum / secpera;
      if (mass_cont_sub_counter == 1) {
        snprintf(tempstr,90, " (dt=%.5f)", major_dt_years);
      } else {
        snprintf(tempstr,90, " (dt=%.5f in %d substeps; av dt_sub_mass_cont=%.5f)",
          major_dt_years, mass_cont_sub_counter, major_dt_years / mass_cont_sub_counter);
      }
      stdout_flags_count0 += tempstr;
      if (delta_t > 0.0) { // avoids printing an empty line if we have not done anything
        stdout_flags_count0 += "\n";
        ierr = verbPrintf(2,grid.com, stdout_flags_count0.c_str()); CHKERRQ(ierr);
      }
      if (stdout_ssa.length() > 0) {
        stdout_ssa += "\n";
        ierr = verbPrintf(2,grid.com, stdout_ssa.c_str()); CHKERRQ(ierr);
      }
      if (do_temp) {
        ierr = verbPrintf(2,grid.com, 
          "S %12.5f: %8.5f %7.4f %10.3f %9.4f\n",
          year, volume/(scale*1.0e9), area/(scale*1.0e6), H0, T0); CHKERRQ(ierr);
      } else {
        ierr = verbPrintf(2,grid.com, 
          "S %12.5f: %8.5f %7.4f %10.3f\n",
          year, volume/(scale*1.0e9), area/(scale*1.0e6), H0); CHKERRQ(ierr);
      }
      mass_cont_sub_counter = 0;      
      mass_cont_sub_dtsum = 0.0;
    }
  }
  return 0;
}


//! \brief Computes cbar, the magnitude of vertically-integrated horizontal
//! velocity of ice and masks out ice-free areas.
PetscErrorCode IceModel::compute_cbar(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscScalar fill_value = -0.01/secpera;

  ierr = vel_bar.magnitude(result); CHKERRQ(ierr);
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas

  ierr = result.set_name("cbar"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", 
			  "magnitude of vertically-integrated horizontal velocity of ice",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes cflx, the magnitude of vertically-integrated horizontal
//! flux of ice.
PetscErrorCode IceModel::compute_cflx(IceModelVec2S &result, IceModelVec2S &cbar) {
  PetscErrorCode ierr;
  PetscScalar fill_value = -0.01/secpera;

  ierr = cbar.multiply_by(vH, result); CHKERRQ(ierr);

  ierr = result.set_name("cflx"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", 
			  "magnitude of vertically-integrated horizontal flux of ice",
			  "m2 s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m2 year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes cbase, the magnitude of horizontal velocity of ice at base
//! of ice and masks out ice-free areas.
//! Uses \c tmp as a preallocated temporary storage.
PetscErrorCode IceModel::compute_cbase(IceModelVec2S &result, IceModelVec2S &tmp) {
  PetscErrorCode ierr;
  PetscScalar fill_value = -0.01/secpera;

  ierr = u3.getHorSlice(result, 0.0); CHKERRQ(ierr); // result = u_{z=0}
  ierr = v3.getHorSlice(tmp, 0.0); CHKERRQ(ierr);    // tmp = v_{z=0}

  ierr = result.set_to_magnitude(result,tmp); CHKERRQ(ierr);

  ierr = result.set_name("cbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", 
			  "magnitude of horizontal velocity of ice at base of ice",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes csurf, the magnitude of horizontal velocity of ice at ice
//! surface and masks out ice-free areas. Uses \c tmp as a preallocated
//! temporary storage.
PetscErrorCode IceModel::compute_csurf(IceModelVec2S &result, IceModelVec2S &tmp) {
  PetscErrorCode ierr;
  PetscScalar fill_value = -0.01/secpera;

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

  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes uvelsurf, the x component of velocity of ice at ice surface.
PetscErrorCode IceModel::compute_uvelsurf(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getSurfaceValues(result, vH); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("uvelsurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "x component of velocity of ice at ice surface",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = 0.0;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", -1e6/secpera); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6/secpera); CHKERRQ(ierr);
//   ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes vvelsurf, the y component of velocity of ice at ice surface.
PetscErrorCode IceModel::compute_vvelsurf(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = v3.getSurfaceValues(result, vH); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("vvelsurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "y component of velocity of ice at ice surface",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = 0.0;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", -1e6/secpera); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6/secpera); CHKERRQ(ierr);
//   ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes vertical ice velocity (relative to the geoid).
/*!
  \f[
  w(s) = \tilde w(s) + \frac{\partial b}{\partial t} + U(s) \cdot \nabla b
  \f]
 */
PetscErrorCode IceModel::compute_wvel(IceModelVec3 &result) {
  PetscErrorCode ierr;
  PetscScalar *res, *w, *u, *v;
  
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = vuplift.begin_access(); CHKERRQ(ierr); 
  ierr = result.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = u3.getInternalColumn(i, j, &u); CHKERRQ(ierr);
      ierr = v3.getInternalColumn(i, j, &v); CHKERRQ(ierr);
      ierr = w3.getInternalColumn(i, j, &w); CHKERRQ(ierr);
      ierr = result.getInternalColumn(i, j, &res); CHKERRQ(ierr);

      for (PetscInt k = 0; k < grid.Mz; ++k)
	res[k] = w[k] + vuplift(i,j) + u[k] * vbed.diff_x_p(i,j) + v[k] * vbed.diff_y_p(i,j);
    }
  }

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vuplift.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr); 

  ierr = result.set_name("wvel"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", 
			  "vertical velocity of ice, relative to geoid",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  ierr = result.set_attr("valid_min", -1e6/secpera); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6/secpera); CHKERRQ(ierr);

  return 0;
}
//! Computes vertical velocity at the base of ice.
PetscErrorCode IceModel::compute_wvelbase(IceModelVec3 &wvel, IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = wvel.begin_access(); CHKERRQ(ierr);
  ierr = wvel.getHorSlice(result, 0.0); CHKERRQ(ierr);
  ierr = wvel.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("wvelbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "vertical velocity of ice at the base of ice, relative to the geoid",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = 0.0;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", -1e6/secpera); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6/secpera); CHKERRQ(ierr);
//   ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes wvelsurf, the vertical velocity of ice at ice surface.
/*! Note that there is no need to mask out ice-free areas here, because
  wvelsurf is zero at those locations.
 */
PetscErrorCode IceModel::compute_wvelsurf(IceModelVec3 &wvel, IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = wvel.begin_access(); CHKERRQ(ierr);
  ierr = wvel.getSurfaceValues(result, vH); CHKERRQ(ierr);
  ierr = wvel.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("wvelsurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "vertical velocity of ice at ice surface, relative to the geoid",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = 0.0;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", -1e6/secpera); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6/secpera); CHKERRQ(ierr);
//   ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes taud, the magnitude of driving shear stress at base of ice.
//! Uses tmp as a preallocated temporary storage.
PetscErrorCode IceModel::compute_taud(IceModelVec2S &result, IceModelVec2S &tmp) {
  PetscErrorCode ierr;

  ierr = computeDrivingStress(result, tmp); CHKERRQ(ierr);

  ierr = result.set_to_magnitude(result, tmp); CHKERRQ(ierr);
  ierr = result.set_name("taud"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic",
			  "magnitude of driving shear stress at base of ice",
			  "Pa", ""); CHKERRQ(ierr);

  PetscScalar fill_value = -0.01;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr); // mask out ice-free areas
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! \brief Sets entrues of result to corresponding processor ranks.
PetscErrorCode IceModel::compute_rank(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j)
      result(i,j) = grid.rank;
  ierr = result.end_access();

  ierr = result.set_name("rank"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "processor rank", "", ""); CHKERRQ(ierr);
  result.time_independent = true;
  return 0;
}

//! \brief Sets entrues of result to corresponding processor ranks.
PetscErrorCode IceModel::compute_proc_ice_area(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscInt ice_filled_cells = 0;

  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j)
      if (vH(i,j) > 0) {
	ice_filled_cells += 1;
      }
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j)
      result(i,j) = ice_filled_cells;
  ierr = result.end_access();

  ierr = result.set_name("proc_ice_area"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic",
			  "number of cells containing ice in a processor's domain",
			  "", ""); CHKERRQ(ierr);
  result.time_independent = true;
  return 0;
}

//! \brief Computes the map of f(|v|) (see [\ref BBssasliding], equation 22)
PetscErrorCode IceModel::compute_bueler_brown_f(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal fill_value = -0.01;

  ierr = vel_ssa.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i,j) = bueler_brown_f(vel_ssa(i,j).magnitude_squared());
    }
  }

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vel_ssa.end_access(); CHKERRQ(ierr);

  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);

  ierr = result.set_name("bueler_brown_f"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "f(|v|) in Bueler and Brown (2009), equation 22", "", ""); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  return 0;
}


//! Compute the CTS field, CTS = E/E_s(p) and put in a global IceModelVec3 provided by user.
PetscErrorCode IceModel::compute_cts(IceModelVec3 &useForCTS) {
  PetscErrorCode ierr;

  ierr = setCTSFromEnthalpy(useForCTS); CHKERRQ(ierr);
  ierr = useForCTS.set_name("cts"); CHKERRQ(ierr);
  ierr = useForCTS.set_attrs(
     "diagnostic",
     "cts = E/E_s(p), so cold-temperate transition surface is at cts = 1",
     "", ""); CHKERRQ(ierr);
  return 0;
}


//! Compute the liquid fraction, and put in a global IceModelVec3 provided by user.
PetscErrorCode IceModel::compute_liqfrac(IceModelVec3 &useForLiqfrac) {
  PetscErrorCode ierr;

  if (config.get_flag("do_cold_ice_methods")) {
    ierr = useForLiqfrac.set(0.0); CHKERRQ(ierr);
  } else {
    ierr = setLiquidFracFromEnthalpy(useForLiqfrac); CHKERRQ(ierr);
  }
  ierr = useForLiqfrac.set_name("liqfrac"); CHKERRQ(ierr);
  ierr = useForLiqfrac.set_attrs(
     "diagnostic","liquid water fraction in ice (between 0 and 1)",
     "", ""); CHKERRQ(ierr);
  return 0;
}


//! \brief Compute the pressure-adjusted temperature in degrees C corresponding
//! to ice temperature, and put in a global IceModelVec3 provided by user.
PetscErrorCode IceModel::compute_temp_pa(IceModelVec3 &result) {
  PetscErrorCode ierr;

  ierr = result.set_name("temp_pa"); CHKERRQ(ierr);
  ierr = result.set_attrs(
     "diagnostic",
     "pressure-adjusted ice temperature (degrees above pressure-melting point)",
     "deg_C", ""); CHKERRQ(ierr);

  PetscScalar *Tpaij, *Enthij; // columns of these values
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result.getInternalColumn(i,j,&Tpaij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k];
	const PetscScalar p = EC->getPressureFromDepth(depth);
        ierr = EC->getPATemp(Enthij[k],p,Tpaij[k]);
          CHKERRQ(ierr);
	  if (config.get_flag("do_cold_ice_methods")) {
	    if ( EC->isTemperate(Enthij[k],p) && (vH(i,j) > 0)) {
	      Tpaij[k] = config.get("water_melting_temperature");
	    }
	  }
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  // make deg C:
  ierr = result.shift(-config.get("water_melting_temperature")); CHKERRQ(ierr);

  // communication not done; we allow global IceModelVec3s as result
  return 0;
}


//! Computes the vertically-averaged ice hardness.
PetscErrorCode IceModel::compute_hardav(IceModelVec2S &result) {
  PetscErrorCode ierr;
  
  const PetscScalar fillval = -0.01;
  
  PetscScalar *Eij; // columns of temperature values
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = Enth3.getInternalColumn(i,j,&Eij); CHKERRQ(ierr);
      const PetscScalar H = vH(i,j);
      if (H > 0.0) {
        result(i,j) = ice->averagedHardness_from_enth(H, grid.kBelowHeight(H),
                                                      grid.zlevels, Eij);
      } else { // put negative value below valid range
        result(i,j) = fillval;
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
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

//! Computes the subglacial (basal) water pressure
/*!
  \f[p_w = \alpha\, \frac{w}{w_{\text{max}}}\, \rho\, g\, H,\f]
  where 

  - \f$\alpha\f$ is the till pore water fraction (till_pw_fraction),
  - \f$w\f$ is the effective thickness of subglacial melt water (bwat)
  - \f$w_{\text{max}}\f$ is the maximum allowed value for \f$w\f$ (hmelt_max),
  - \f$\rho\f$ is the ice density (ice_density)
  - \f$H\f$ is the ice thickness (thk)

Result is set to invalid (_FillValue) where the ice is floating, there being
no meaning to the above calculation.
 */
PetscErrorCode IceModel::compute_bwp(IceModelVec2S &result) {
  PetscErrorCode ierr;
  const PetscScalar
    alpha     = config.get("till_pw_fraction"),
    wmax      = config.get("hmelt_max"),
    fillval   = -0.01;

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHmelt.begin_access(); CHKERRQ(ierr);
  ierr = vbmr.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0.0) {
        result(i,j) = getBasalWaterPressure(
                        vH(i,j), vHmelt(i,j), vbmr(i,j), alpha, wmax);
      } else { // put negative value below valid range
        result(i,j) = fillval;
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vHmelt.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = vMask.fill_where_floating(result, fillval); CHKERRQ(ierr);

  ierr = result.set_name("bwp"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "subglacial (pore) water pressure",
			  "Pa", ""); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min",0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue",fillval); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes the rate of change of ice surface elevation as a sum of the
//! bedrock uplift rate and the thickness rate of change.
PetscErrorCode IceModel::compute_dhdt(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = result.copy_from(vdHdt); CHKERRQ(ierr); // result = dHdt
  ierr = result.mask_by(vH,0.0); CHKERRQ(ierr);	// set _FillValue areas to 0.0
  ierr = result.add(1.0, vuplift); CHKERRQ(ierr); // result += dbdt = uplift

  ierr = result.set_name("dhdt"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "rate of change of surface elevation",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1");
  result.write_in_glaciological_units = true;

  return 0;
}

PetscErrorCode IceModel::compute_uvelbase(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getHorSlice(result, 0.0); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("uvelbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "x component of ice velocity at the base of ice",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = 0.0;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", -1e6/secpera); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6/secpera); CHKERRQ(ierr);
//   ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::compute_vvelbase(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = v3.getHorSlice(result, 0.0); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("vvelbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "y component of ice velocity at the base of ice",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  result.write_in_glaciological_units = true;

  PetscScalar fill_value = 0.0;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", -1e6/secpera); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 1e6/secpera); CHKERRQ(ierr);
//   ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes the thickness of the basal layer of the temperate ice.
/*!
 * Uses linear interpolation to go beyond vertical grid resolution.
 */
PetscErrorCode IceModel::compute_tempicethk_basal(IceModelVec2S& result) {
  PetscErrorCode ierr;
  PetscScalar *Enth, fill_value = -0.01;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {

      // if we have no ice, go on to the next grid point (this cell will be
      // marked as "missing" later)
      if (vH(i,j) < 0.1)
        continue;

      ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
      PetscReal pressure;
      PetscInt ks = grid.kBelowHeight(vH(i,j)),
        k = 0;

      while (k <= ks) {
        pressure = EC->getPressureFromDepth(vH(i,j) - grid.zlevels[k]);

        if (EC->isTemperate(Enth[k],pressure))
          k++;
        else
          break;
      }
      // after this loop 'pressure' is equal to the pressure at the first level
      // that is cold

      // no temperate ice at all; go to the next grid point
      if (k == 0) {
        result(i,j) = 0;
        continue;
      }
      
      // the whole column is temperate (except, possibly, some ice between
      // zlevels[ks] and the total thickness; we ignore it)
      if (k == ks + 1) {
        result(i,j) = grid.zlevels[ks];
        continue;
      }

      PetscReal 
        pressure_0 = EC->getPressureFromDepth(vH(i,j) - grid.zlevels[k-1]),
        dz = grid.zlevels[k] - grid.zlevels[k-1],
        slope1 = (Enth[k] - Enth[k-1]) / dz,
        slope2 = (EC->getEnthalpyCTS(pressure) - EC->getEnthalpyCTS(pressure_0)) / dz;
      
      if (slope1 != slope2) {
        result(i,j) = grid.zlevels[k-1] +
          (EC->getEnthalpyCTS(pressure_0) - Enth[k-1]) / (slope1 - slope2);

        // check if the resulting thickness is valid:
        result(i,j) = PetscMax(result(i,j), grid.zlevels[k-1]);
        result(i,j) = PetscMin(result(i,j), grid.zlevels[k]);
      } else {
        SETERRQ4(1, "This should never happen: (i=%d, j=%d, k=%d, ks=%d)\n",
                    i, j, k, ks);
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("tempicethk_basal"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "thickness of the basal layer of temperate ice",
			  "m", ""); CHKERRQ(ierr);

  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes thickness of the temperate ice layer (if there is any).
PetscErrorCode IceModel::compute_tempicethk(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscScalar *Enth;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0.) {
	ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
	PetscScalar tithk = 0.;
	const PetscInt ks = grid.kBelowHeight(vH(i,j));
        
	for (PetscInt k=0; k<ks; ++k) {
          PetscReal pressure = EC->getPressureFromDepth(vH(i,j) - grid.zlevels[k]);

	  if (EC->isTemperate(Enth[k], pressure)) {
	    tithk += grid.zlevels[k+1] - grid.zlevels[k];
	  }
	}

	if (EC->isTemperate(Enth[ks],EC->getPressureFromDepth(vH(i,j) - grid.zlevels[ks]))) {
	  tithk += vH(i,j) - grid.zlevels[ks];
	}

	result(i,j) = tithk;
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("tempicethk"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "temperate ice thickness (total column content)",
			  "m", ""); CHKERRQ(ierr);

  PetscScalar fill_value = -0.01;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes ice enthalpy at the base of ice.
PetscErrorCode IceModel::compute_enthalpybase(IceModelVec2S &result) {
  PetscErrorCode ierr;

  // put basal ice temperature in vWork2d[0]
  ierr = Enth3.getHorSlice(result, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = result.set_name("enthalpybase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "ice enthalpy at the base of ice",
			  "J kg-1", ""); CHKERRQ(ierr);

  PetscScalar fill_value = -0.01;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes ice temperature at the base of ice.
PetscErrorCode IceModel::compute_tempbase(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = compute_temp(vWork3d); CHKERRQ(ierr);

  // put basal ice temperature in vWork2d[0]
  ierr = vWork3d.getHorSlice(result, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = result.set_name("tempbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "ice temperature at the base of ice",
			  "Kelvin", ""); CHKERRQ(ierr);

  PetscScalar fill_value = -0.01;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes pressure-adjusted ice temperature at the base of ice.
PetscErrorCode IceModel::compute_temppabase(IceModelVec3 &hasPATemp,
                                            IceModelVec2S &result) {
  PetscErrorCode ierr;

  // put basal pressure-adjusted ice temperature in 2d result
  ierr = hasPATemp.getHorSlice(result, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = result.set_name("temppabase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", 
                          "pressure-adjusted ice temperature at the base of ice",
			  "degrees Celsius", ""); CHKERRQ(ierr);

  PetscScalar fill_value = 0.01;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes ice temperature at the 1 m below the surface.
PetscErrorCode IceModel::compute_tempsurf(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscScalar fill_value = -0.01;

  // compute levels corresponding to 1 m below the ice surface:

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      result(i,j) = PetscMax(vH(i,j) - 1.0, 0.0);
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = compute_temp(vWork3d); CHKERRQ(ierr);

  ierr = vWork3d.getSurfaceValues(result, result); CHKERRQ(ierr);  // z=H-1 slice

  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) <= 1.0)
	result(i,j) = fill_value;
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

//! Computes ice enthalpy at the 1 m below the surface.
PetscErrorCode IceModel::compute_enthalpysurf(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscScalar fill_value = -0.01;

  // compute levels corresponding to 1 m below the ice surface:

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      result(i,j) = PetscMax(vH(i,j) - 1.0, 0.0);
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = Enth3.getSurfaceValues(result, result); CHKERRQ(ierr);  // z=0 slice

  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) <= 1.0)
	result(i,j) = fill_value;
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("tempsurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "ice enthalpy at 1m below the ice surface",
			  "J kg-1", ""); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);
  return 0;
}


//! \brief Computes the multiplier \f$\theta\f$ in Schoof's (2003) theory of the
//! effect of bed roughness on the diffusivity of the SIA.
/*!
See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
 */
PetscErrorCode IceModel::compute_schoofs_theta(IceModelVec2S &result) {
  PetscErrorCode ierr;

  if (sia_bed_smoother==NULL) {
    SETERRQ(1,"PISM ERROR: sia_bed_smoother==NULL in compute_schoofs_theta()");
  }
  
  ierr = sia_bed_smoother->preprocess_bed(vbed,
               config.get("Glen_exponent"), config.get("bed_smoother_range") );
               CHKERRQ(ierr);
  ierr = sia_bed_smoother->get_theta(
               vh, config.get("Glen_exponent"), 0, &result); CHKERRQ(ierr);

  ierr = result.set_name("schoofs_theta"); CHKERRQ(ierr);
  ierr = result.set_attrs(
        "diagnostic", 
        "multiplier 'theta' in Schoof's (2003) theory of bed roughness in SIA",
        "", ""); CHKERRQ(ierr);

  ierr = result.set_attr("valid_max", 1.0); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  const PetscScalar fill_value = -1.0;
  ierr = result.mask_by(vH,fill_value); CHKERRQ(ierr);	// set no-ice to _FillValue
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);
  return 0;
}


//! \brief Computes the smoothed bed elevation from Schoof's (2003) theory of the
//! effect of bed roughness on the SIA.
/*!
See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
 */
PetscErrorCode IceModel::compute_topgsmooth(IceModelVec2S &result) {
  PetscErrorCode ierr;

  if (sia_bed_smoother==NULL) {
    SETERRQ(1,"PISM ERROR: sia_bed_smoother==NULL in compute_topgsmooth()");
  }
  
  ierr = sia_bed_smoother->preprocess_bed(vbed,
               config.get("Glen_exponent"), config.get("bed_smoother_range") );
               CHKERRQ(ierr);
  ierr = result.copy_from(sia_bed_smoother->topgsmooth); CHKERRQ(ierr);

  ierr = result.set_name("topgsmooth"); CHKERRQ(ierr);
  ierr = result.set_attrs(
        "diagnostic", 
        "smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
        "m", ""); CHKERRQ(ierr);
  return 0;
}


//! \brief Computes the thickness relative to the smoothed bed elevation in
//! Schoof's (2003) theory of the effect of bed roughness on the SIA.
/*!
See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
 */
PetscErrorCode IceModel::compute_thksmooth(IceModelVec2S &result) {
  PetscErrorCode ierr;

  if (sia_bed_smoother==NULL) {
    SETERRQ(1,"PISM ERROR: sia_bed_smoother==NULL in compute_thksmooth()");
  }
  
  ierr = sia_bed_smoother->preprocess_bed(vbed,
               config.get("Glen_exponent"), config.get("bed_smoother_range") );
               CHKERRQ(ierr);
  ierr = sia_bed_smoother->get_smoothed_thk(vh, vH, 0, &result); CHKERRQ(ierr);

  ierr = result.set_name("thksmooth"); CHKERRQ(ierr);
  ierr = result.set_attrs(
        "diagnostic", 
        "thickness relative to smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
        "m", ""); CHKERRQ(ierr);
  //ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  return 0;
}


//! \brief Computes a diagnostic quantity given by \c name and returns a
//! pointer to a pre-allocated work vector containing it.
/*! 
For 2D quantities, result will point to vWork2d[0].  For 3D -- to vWork3d.
Depending on the quantity requested, vWork2d[1] might get used as temporary
storage.
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

  if (name == "cts") {
    ierr = compute_cts(vWork3d); CHKERRQ(ierr);
    result = &vWork3d;
    return 0;
  }

  if (name == "dhdt") {
    ierr = compute_dhdt(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "enthalpybase") {
    ierr = compute_enthalpybase(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "enthalpysurf") {
    ierr = compute_enthalpysurf(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "hardav") {
    ierr = compute_hardav(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "liqfrac") {
    ierr = compute_liqfrac(vWork3d); CHKERRQ(ierr);
    result = &vWork3d;
    return 0;
  }

  if (name == "schoofs_theta") {
    ierr = compute_schoofs_theta(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "taud") {
    ierr = compute_taud(vWork2d[0], vWork2d[1]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "temp") {
    ierr = compute_temp(vWork3d); CHKERRQ(ierr);
    result = &vWork3d;
    return 0;
  }

  if (name == "temp_pa") {
    ierr = compute_temp_pa(vWork3d); CHKERRQ(ierr);
    result = &vWork3d;
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

  if (name == "temppabase") {
    ierr = compute_temp_pa(vWork3d); CHKERRQ(ierr);
    ierr = compute_temppabase(vWork3d,vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "tempicethk") {
    ierr = compute_tempicethk(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "tempicethk_basal") {
    ierr = compute_tempicethk_basal(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "thksmooth") {
    ierr = compute_thksmooth(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "topgsmooth") {
    ierr = compute_topgsmooth(vWork2d[0]); CHKERRQ(ierr);
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

  if (name == "wvel") {
    ierr = compute_wvel(vWork3d); CHKERRQ(ierr);
    result = &vWork3d;
    return 0;
  }

  if (name == "wvelbase") {
    ierr = compute_wvel(vWork3d); CHKERRQ(ierr);
    ierr = compute_wvelbase(vWork3d, vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "wvelsurf") {
    ierr = compute_wvel(vWork3d); CHKERRQ(ierr);
    ierr = compute_wvelsurf(vWork3d, vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "rank") {
    ierr = compute_rank(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "proc_ice_area") {
    ierr = compute_proc_ice_area(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "diffusivity") {
    ierr = compute_diffusivity(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "bueler_brown_f") {
    ierr = compute_bueler_brown_f(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  return 0;
}


//! Computes the ice volume, in m^3.
PetscErrorCode IceModel::compute_ice_volume(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     volume=0.0;
  
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0)
        volume += vH(i,j) * cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
  return 0;
}


//! Computes the temperate ice volume, in m^3.
PetscErrorCode IceModel::compute_ice_volume_temperate(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     volume=0.0;
  
  PetscScalar *Enth;  // do NOT delete this pointer: space returned by
                      //   getInternalColumn() is allocated already
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
        const PetscInt ks = grid.kBelowHeight(vH(i,j));
        ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
        for (PetscInt k=0; k<ks; ++k) {	  
          if (EC->isTemperate(Enth[k],EC->getPressureFromDepth(vH(i,j)))) {
            volume += (grid.zlevels[k+1] - grid.zlevels[k]) * cell_area(i,j);
          }
        }
        if (EC->isTemperate(Enth[ks],EC->getPressureFromDepth(vH(i,j)))) {
          volume += (vH(i,j) - grid.zlevels[ks]) * cell_area(i,j);
        }
      }
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes the cold ice volume, in m^3.
PetscErrorCode IceModel::compute_ice_volume_cold(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     volume=0.0;
  
  PetscScalar *Enth;  // do NOT delete this pointer: space returned by
                      //   getInternalColumn() is allocated already
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
        const PetscInt ks = grid.kBelowHeight(vH(i,j));
        ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
        for (PetscInt k=0; k<ks; ++k) {	  
          if (!EC->isTemperate(Enth[k],EC->getPressureFromDepth(vH(i,j)))) {
            volume += (grid.zlevels[k+1] - grid.zlevels[k]) * cell_area(i,j);
          }
        }
        if (!EC->isTemperate(Enth[ks],EC->getPressureFromDepth(vH(i,j)))) {
          volume += (vH(i,j) - grid.zlevels[ks]) * cell_area(i,j);
        }
      }
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0)
        area += cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes temperate ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_temperate(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  
  ierr = Enth3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);  // z=0 slice
  PetscScalar **Enthbase;
  ierr = vWork2d[0].get_array(Enthbase); CHKERRQ(ierr);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ( (vH(i,j) > 0) && (EC->isTemperate(Enthbase[i][j],EC->getPressureFromDepth(vH(i,j)))) )
        area += cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes cold ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_cold(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  
  ierr = Enth3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);  // z=0 slice
  PetscScalar **Enthbase;
  ierr = vWork2d[0].get_array(Enthbase); CHKERRQ(ierr);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ( (vH(i,j) > 0) && (!EC->isTemperate(Enthbase[i][j],EC->getPressureFromDepth(vH(i,j)))) )
        area += cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes grounded ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_grounded(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ( (vH(i,j) > 0) && vMask.is_grounded(i,j) )
        area += cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes floating ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_floating(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ( (vH(i,j) > 0) && vMask.is_floating(i,j) )
        area += cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}


//! Computes the total ice enthalpy in J.
/*!
Units of the specific enthalpy field \f$E=\f$(IceModelVec3::Enth3) are J kg-1.  We integrate
\f$E(t,x,y,z)\f$ over the entire ice fluid region \f$\Omega(t)\f$, multiplying
by the density to get units of energy:
   \f[ E_{\text{total}}(t) = \int_{\Omega(t)} E(t,x,y,z) \rho_i \,dx\,dy\,dz. \f]
 */
PetscErrorCode IceModel::compute_ice_enthalpy(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar enthalpysum = 0.0;

  PetscScalar *Enth;  // do NOT delete this pointer: space returned by
                      //   getInternalColumn() is allocated already
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
        const PetscInt ks = grid.kBelowHeight(vH(i,j));
        ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
        for (PetscInt k=0; k<ks; ++k) {
          enthalpysum += Enth[k] * (grid.zlevels[k+1] - grid.zlevels[k]);
        }
        enthalpysum += Enth[ks] * (vH(i,j) - grid.zlevels[ks]);
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);  
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  
  enthalpysum *= config.get("ice_density") * (grid.dx * grid.dy);
  
  ierr = PetscGlobalSum(&enthalpysum, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Compute diffusivity associated to the SIA mass continuity equation (for diagnostic purposes).
/*! 
Though inefficient (perhaps), we simply re-compute the SIA velocity update.
This at least means that all quantities are updated using the same settings.

Note that the SIA update puts the diffusivity on the staggered grid, so we do
have to average it onto the regular grid.  But at least velocitySIAStaggered()
fills in ghosted points so no communication is needed.

Uses  vWork2d[0,1,2,3,4,5] indirectly, in the above calls.
 */
PetscErrorCode IceModel::compute_diffusivity(IceModelVec2S &result) {
  PetscErrorCode ierr;

  // compare to calls in IceModel::velocity()
  ierr = surfaceGradientSIA(); CHKERRQ(ierr);
  ierr = velocitySIAStaggered(); CHKERRQ(ierr); // ends with staggered
                                                //   diffusivity in vWork2dStag
  ierr = vWork2dStag.staggered_to_regular(result); CHKERRQ(ierr); 

  ierr = result.set_name("diffusivity"); CHKERRQ(ierr);
  ierr = result.set_attrs(
        "diagnostic", 
        "diffusivity of SIA mass continuity equation",
        "m2 s-1", ""); CHKERRQ(ierr);

  return 0;
}

//! Compute a scalar diagnostic quantity by name.
PetscErrorCode IceModel::compute_by_name(string name, PetscScalar &result) {
  PetscErrorCode ierr, errcode = 1;

  if (name == "ivol") {
    errcode = 0;
    ierr = compute_ice_volume(result); CHKERRQ(ierr);
  }

  if (name == "ivoltemp") {
    errcode = 0;
    ierr = compute_ice_volume_temperate(result); CHKERRQ(ierr);
  }

  if (name == "ivoltempf") {
    errcode = 0;
    PetscScalar ivol;
    ierr = compute_ice_volume(ivol); CHKERRQ(ierr);
    ierr = compute_ice_volume_temperate(result); CHKERRQ(ierr);
    result /= ivol;
  }

  if (name == "ivolcold") {
    errcode = 0;
    ierr = compute_ice_volume_cold(result); CHKERRQ(ierr);
  }

  if (name == "ivolcoldf") {
    errcode = 0;
    PetscScalar ivol;
    ierr = compute_ice_volume(ivol); CHKERRQ(ierr);
    ierr = compute_ice_volume_cold(result); CHKERRQ(ierr);
    result /= ivol;
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

  if (name == "iareatemp") {
    errcode = 0;
    ierr = compute_ice_area_temperate(result); CHKERRQ(ierr);
  }

  if (name == "iareatempf") {
    errcode = 0;
    PetscScalar iarea;
    ierr = compute_ice_area(iarea); CHKERRQ(ierr);
    ierr = compute_ice_area_temperate(result); CHKERRQ(ierr);
    result /= iarea;
  }

  if (name == "iareacold") {
    errcode = 0;
    ierr = compute_ice_area_cold(result); CHKERRQ(ierr);
  }

  if (name == "iareacoldf") {
    errcode = 0;
    PetscScalar iarea;
    ierr = compute_ice_area(iarea); CHKERRQ(ierr);
    ierr = compute_ice_area_cold(result); CHKERRQ(ierr);
    result /= iarea;
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

  if (name == "ienthalpy") {
    errcode = 0;
    ierr = compute_ice_enthalpy(result); CHKERRQ(ierr);
  }

  if (name == "basal_ice_flux") {
    errcode = 0;
    result = total_basal_ice_flux;
  }

  if (name == "surface_ice_flux") {
    errcode = 0;
    result = total_surface_ice_flux;
  }

  if (name == "sub_shelf_ice_flux") {
    errcode = 0;
    result = total_sub_shelf_ice_flux;
  }

  if (name == "nonneg_rule_flux") {
    errcode = 0;
    result = nonneg_rule_flux;
  }

  if (name == "ocean_kill_flux") {
    errcode = 0;
    result = ocean_kill_flux;
  }

  if (name == "float_kill_flux") {
    errcode = 0;
    result = float_kill_flux;
  }

  return errcode;
}


/*! Computes total ice fluxes in kg s-1 at 3 interfaces:

  \li the ice-atmosphere interface: gets surface mass balance rate from
      PISMSurfaceModel *surface,
  \li the ice-ocean interface at the bottom of ice shelves: gets ocean-imposed
      basal melt rate from PISMOceanModel *ocean, and
  \li the ice-bedrock interface: gets basal melt rate from IceModelVec2S vbmr.

A unit-conversion occurs for all three quantities, from ice-equivalent m s-1
to kg s-1.  The sign convention about these fluxes is that positve flux means
ice is being \e added to the ice fluid volume at that interface.

These quantities should be understood as <i>instantaneous at the beginning of
the time-step.</i>  Multiplying by dt will \b not necessarily give the
corresponding change from the beginning to the end of the time-step.

FIXME:  The calving rate can be computed by post-processing:
divoldt = surface_ice_flux * iarea + basal_ice_flux * iareag + sub_shelf_ice_flux * iareaf + calving_flux_vol_rate
 */
PetscErrorCode IceModel::ice_mass_bookkeeping() {
  PetscErrorCode ierr;

  bool include_bmr_in_continuity = config.get_flag("include_bmr_in_continuity");

  // note acab and shelfbmassflux are IceModelVec2S owned by IceModel
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_mass_flux(grid.year, dt / secpera, acab);
    CHKERRQ(ierr);
  } else { SETERRQ(2,"PISM ERROR: surface == PETSC_NULL"); }

  if (ocean != PETSC_NULL) {
    ierr = ocean->shelf_base_mass_flux(grid.year, dt / secpera, shelfbmassflux);
    CHKERRQ(ierr);
  } else { SETERRQ(2,"PISM ERROR: ocean == PETSC_NULL"); }

  PetscScalar my_total_surface_ice_flux = 0.0, my_total_basal_ice_flux = 0.0,
    my_total_sub_shelf_ice_flux = 0.0;

  ierr = acab.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vbmr.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // ignore ice-free cells:
      if (vH(i,j) <= 0.0)
        continue;

      my_total_surface_ice_flux += acab(i,j) * cell_area(i,j); // note the "+="!

      if ((vMask.value(i,j) == MASK_FLOATING) && include_bmr_in_continuity) {
        // note: we are deliberately *not* including fluxes in
        //   MASK_ICE_FREE_OCEAN and MASK_OCEAN_AT_TIME_0 areas
        my_total_sub_shelf_ice_flux -= shelfbmassflux(i,j) * cell_area(i,j); // note the "-="!
      }

      if (vMask.is_grounded(i,j) && include_bmr_in_continuity) {
        my_total_basal_ice_flux -= vbmr(i,j) * cell_area(i,j); // note the "-="!
      }
    }	// j
  } // i

  ierr = acab.end_access(); CHKERRQ(ierr);
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);

  PetscScalar ice_density = config.get("ice_density");
  my_total_surface_ice_flux     *= ice_density;
  my_total_sub_shelf_ice_flux   *= ice_density;
  my_total_basal_ice_flux       *= ice_density;

  ierr = PetscGlobalSum(&my_total_surface_ice_flux,   &total_surface_ice_flux,   grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&my_total_sub_shelf_ice_flux, &total_sub_shelf_ice_flux, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&my_total_basal_ice_flux,     &total_basal_ice_flux,     grid.com); CHKERRQ(ierr);

  return 0;
}


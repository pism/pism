// Copyright (C) 2004-2012 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <petscsys.h>
#include <stdarg.h>
#include <stdlib.h>

#include "iceModel.hh"
#include "Mask.hh"
#include "PISMStressBalance.hh"
#include "PISMOcean.hh"
#include "enthalpyConverter.hh"
#include "PISMTime.hh"

//!  Computes volume and area of ice sheet, for reporting purposes.
/*!
  Communication done for global max and global sum.

  Returns area in units of m^2 and volume in m^3.
 */
PetscErrorCode IceModel::volumeArea(PetscScalar& gvolume, PetscScalar& garea) {

  PetscErrorCode  ierr;
  PetscScalar     volume=0.0, area=0.0;

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  MaskQuery mask(vMask);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
        area += cell_area(i,j);
        const PetscScalar dv = cell_area(i,j) * vH(i,j);
        volume += dv;
      }
    }
  }

  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&volume, &gvolume, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&area, &garea, grid.com); CHKERRQ(ierr);
  return 0;
}


/*!
  Computes fraction of the base which is melted.

  Communication occurs here.

  FIXME: energyStats should use cell_area(i,j).
 */
PetscErrorCode IceModel::energyStats(PetscScalar iarea, PetscScalar &gmeltfrac) {
  PetscErrorCode    ierr;
  PetscScalar       meltarea = 0.0, temp0 = 0.0;
  const PetscScalar a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
  IceModelVec2S &Enthbase = vWork2d[0];

  // use Enth3 to get stats
  ierr = Enth3.getHorSlice(Enthbase, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Enthbase.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
	// accumulate area of base which is at melt point
	if (EC->isTemperate(Enthbase(i,j), EC->getPressureFromDepth(vH(i,j)) )) // FIXME issue #15
	  meltarea += a;
      }
      // if you happen to be at center, record absolute basal temp there
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
	ierr = EC->getAbsTemp(Enthbase(i,j),EC->getPressureFromDepth(vH(i,j)), temp0); // FIXME issue #15
	CHKERRQ(ierr);
      }
    }
  }
  ierr = Enthbase.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  // communication
  ierr = PISMGlobalSum(&meltarea, &gmeltfrac, grid.com); CHKERRQ(ierr);

  // normalize fraction correctly
  if (iarea > 0.0)   gmeltfrac = gmeltfrac / iarea;
  else gmeltfrac = 0.0;

  return 0;
}


/*!
  Computes fraction of the ice which is as old as the start of the run (original).
  Communication occurs here.

  FIXME: ageStats should use cell_area(i,j).
 */
PetscErrorCode IceModel::ageStats(PetscScalar ivol, PetscScalar &gorigfrac) {
  PetscErrorCode  ierr;

  gorigfrac = -1.0;  // result value if not do_age

  if (!config.get_flag("do_age"))
    return 0;  // leave now

  const PetscScalar  a = grid.dx * grid.dy * 1e-3 * 1e-3, // area unit (km^2)
    currtime = grid.time->current(); // seconds

  PetscScalar *tau, origvol = 0.0;
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = tau3.begin_access(); CHKERRQ(ierr);

  // compute local original volume
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
        // accumulate volume of ice which is original
        ierr = tau3.getInternalColumn(i,j,&tau); CHKERRQ(ierr);
        const PetscInt  ks = grid.kBelowHeight(vH(i,j));
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
  ierr = PISMGlobalSum(&origvol,  &gorigfrac, grid.com); CHKERRQ(ierr);

  // normalize fraction correctly
  if (ivol > 0.0)    gorigfrac = gorigfrac / ivol;
  else gorigfrac = 0.0;

  return 0;
}


PetscErrorCode IceModel::summary(bool tempAndAge) {
  PetscErrorCode  ierr;
  PetscScalar     gvolume, garea;
  PetscScalar     meltfrac = 0.0, origfrac = 0.0;
  PetscScalar     max_diffusivity;

  // get volumes in m^3 and areas in m^2
  ierr = volumeArea(gvolume, garea); CHKERRQ(ierr);

  if (tempAndAge || (getVerbosityLevel() >= 3)) {
    ierr = energyStats(garea, meltfrac); CHKERRQ(ierr);
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

  // get maximum diffusivity
  ierr = stress_balance->get_max_diffusivity(max_diffusivity); CHKERRQ(ierr);

  // main report: 'S' line
  ierr = summaryPrintLine(PETSC_FALSE, tempAndAge, grid.time->date(), dt,
                          gvolume,garea,meltfrac,max_diffusivity); CHKERRQ(ierr);

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
  - \c max_diff is the maximum diffusivity
thick0 and temp0 can be interpreted as "sanity checks", because they give
information about a location which may or may not be "typical".

For more description and examples, see the PISM User's Manual.
Derived classes of IceModel may redefine this method and print alternate
information.  Use of DiagnosticTimeseries may be superior, however.
 */
PetscErrorCode IceModel::summaryPrintLine(PetscBool printPrototype,  bool tempAndAge,
                                          string date,  PetscScalar delta_t,
                                          PetscScalar volume,  PetscScalar area,
                                          PetscScalar /* meltfrac */,  PetscScalar max_diffusivity) {

  PetscErrorCode ierr;
  const bool do_energy = config.get_flag("do_energy");

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
    ierr = verbPrintf(2,grid.com,
                      "P       YEAR:       ivol      iarea  max_diffusivity  max_hor_vel\n");
    ierr = verbPrintf(2,grid.com,
                      "U      years   %skm^3  %skm^2         m^2 s^-1       m/year\n",
                      volscalestr,areascalestr);
    return 0;
  }

  if (mass_cont_sub_counter == 0)
    stdout_flags_count0 = stdout_flags;
  if (delta_t > 0.0) {
    mass_cont_sub_counter++;
    mass_cont_sub_dtsum += delta_t;
  }

  if ((tempAndAge == PETSC_TRUE) || (!do_energy) || (getVerbosityLevel() > 2)) {
    char tempstr[90] = "";
    const PetscScalar major_dt_years = convert(mass_cont_sub_dtsum, "seconds", "years");

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

    if (stdout_ssa.empty() == false) {
      ierr = verbPrintf(2, grid.com, "%s\n", stdout_ssa.c_str()); CHKERRQ(ierr);
    }

    ierr = verbPrintf(2,grid.com,
                      "S %s: %8.5f %9.5f %12.5f %14.5f\n",
                      date.c_str(), volume/(scale*1.0e9), area/(scale*1.0e6), max_diffusivity,
                      convert(gmaxu > gmaxv ? gmaxu : gmaxv, "m/s", "m/year")); CHKERRQ(ierr);

    mass_cont_sub_counter = 0;
    mass_cont_sub_dtsum = 0.0;
  }

  return 0;
}

  //! Computes the ice volume, in m^3.
PetscErrorCode IceModel::compute_ice_volume(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     volume=0.0;

  ierr = cell_area.begin_access(); CHKERRQ(ierr);

  {
    ierr = vH.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (vH(i,j) > 0)
          volume += vH(i,j) * cell_area(i,j);
      }
    }
    ierr = vH.end_access(); CHKERRQ(ierr);
  }

  // Add the volume of the ice in Href:
  if (config.get_flag("part_grid")) {
    ierr = vHref.begin_access(); CHKERRQ(ierr);
    for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
        volume += vHref(i,j) * cell_area(i,j);
      }
    }
    ierr = vHref.end_access(); CHKERRQ(ierr);
  }

  ierr = cell_area.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes the ice volume, which is relevant for sea-level rise in m^3 in SEA-WATER EQUIVALENT.
PetscErrorCode IceModel::compute_sealevel_volume(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     volume=0.0;
  MaskQuery mask(vMask);
  double ocean_rho = config.get("sea_water_density");
  double ice_rho = config.get("ice_density");

  if (ocean == PETSC_NULL) {  SETERRQ(grid.com, 1, "PISM ERROR: ocean == PETSC_NULL");  }
  PetscReal sea_level;
  ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);

  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access();  CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (mask.grounded_ice(i,j)){
        if (vH(i,j) > 0) {
          if(vbed(i, j) > sea_level){
            volume += vH(i,j) * cell_area(i,j) * ice_rho/ocean_rho ;
          } else {
            volume += vH(i,j) * cell_area(i,j) * ice_rho/ocean_rho - cell_area(i,j) * ( sea_level - vbed(i, j) );
          }
        }
      }
    }
  }
  const PetscScalar oceanarea=3.61e14;//in square meters
  volume /= oceanarea;
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
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
          if (EC->isTemperate(Enth[k],EC->getPressureFromDepth(vH(i,j)))) { // FIXME issue #15
            volume += (grid.zlevels[k+1] - grid.zlevels[k]) * cell_area(i,j);
          }
        }
        if (EC->isTemperate(Enth[ks],EC->getPressureFromDepth(vH(i,j)))) { // FIXME issue #15
          volume += (vH(i,j) - grid.zlevels[ks]) * cell_area(i,j);
        }
      }
    }
  }
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
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
          if (!EC->isTemperate(Enth[k],EC->getPressureFromDepth(vH(i,j)))) { // FIXME issue #15
            volume += (grid.zlevels[k+1] - grid.zlevels[k]) * cell_area(i,j);
          }
        }
        if (!EC->isTemperate(Enth[ks],EC->getPressureFromDepth(vH(i,j)))) { // FIXME issue #15
          volume += (vH(i,j) - grid.zlevels[ks]) * cell_area(i,j);
        }
      }
    }
  }
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
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

  ierr = PISMGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes area of basal ice which is temperate, in m^2.
PetscErrorCode IceModel::compute_ice_area_temperate(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  IceModelVec2S &Enthbase = vWork2d[0];

  ierr = Enth3.getHorSlice(Enthbase, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = Enthbase.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ( (vH(i,j) > 0) && (EC->isTemperate(Enthbase(i,j),EC->getPressureFromDepth(vH(i,j)))) ) // FIXME issue #15
        area += cell_area(i,j);
    }
  }
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = Enthbase.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes area of basal ice which is cold, in m^2.
PetscErrorCode IceModel::compute_ice_area_cold(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  IceModelVec2S &Enthbase = vWork2d[0];

  ierr = Enth3.getHorSlice(Enthbase, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = Enthbase.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ( (vH(i,j) > 0) && (!EC->isTemperate(Enthbase(i,j),EC->getPressureFromDepth(vH(i,j)))) ) // FIXME issue #15
        area += cell_area(i,j);
    }
  }
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = Enthbase.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes grounded ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_grounded(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;

  MaskQuery mask(vMask);

  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (mask.grounded_ice(i,j))
        area += cell_area(i,j);
    }
  }
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes floating ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_floating(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;

  MaskQuery mask(vMask);

  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (mask.floating_ice(i,j))
        area += cell_area(i,j);
    }
  }
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
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

  ierr = PISMGlobalSum(&enthalpysum, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

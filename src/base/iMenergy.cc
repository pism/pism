// Copyright (C) 2004-2011, 2013, 2014 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "iceModel.hh"
#include "Mask.hh"
#include "bedrockThermalUnit.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "enthalpyConverter.hh"
#include <assert.h>

namespace pism {

//! \file iMenergy.cc Methods of IceModel which address conservation of energy.
//! Common to enthalpy (polythermal) and temperature (cold-ice) methods.

//! Manage the solution of the energy equation, and related parallel communication.
/*!
  This method updates three fields:
  - IceModelVec3 Enth3
  - IceModelVec2 basal_melt_rate
  - IceModelVec2 vHmelt
  That is, energyStep() is in charge of calling other methods that actually update, and
  then it is in charge of doing the ghost communication as needed.  If
  do_cold_ice_methods == true, then energyStep() must also update this field
  - IceModelVec3 T3

  Normally calls the method enthalpyAndDrainageStep().  Calls temperatureStep() if
  do_cold_ice_methods == true.
*/
PetscErrorCode IceModel::energyStep() {
  PetscErrorCode  ierr;

  double  myCFLviolcount = 0.0,   // these are counts but they are type "double"
    myVertSacrCount = 0.0,  //   because that type works with GlobalSum()
    myBulgeCount = 0.0;
  double gVertSacrCount, gBulgeCount;

  // always count CFL violations for sanity check (but can occur only if -skip N with N>1)
  ierr = countCFLViolations(&myCFLviolcount); CHKERRQ(ierr);

  // operator-splitting occurs here (ice and bedrock energy updates are split):
  //   tell BedThermalUnit* btu that we have an ice base temp; it will return
  //   the z=0 value of geothermal flux when called inside temperatureStep() or
  //   enthalpyAndDrainageStep()
  ierr = get_bed_top_temp(bedtoptemp); CHKERRQ(ierr);
  ierr = btu->update(t_TempAge, dt_TempAge); CHKERRQ(ierr);  // has ptr to bedtoptemp

  if (config.get_flag("do_cold_ice_methods")) {
    // new temperature values go in vTnew; also updates Hmelt:
    ierr = temperatureStep(&myVertSacrCount,&myBulgeCount); CHKERRQ(ierr);

    ierr = vWork3d.update_ghosts(T3); CHKERRQ(ierr);

    // compute_enthalpy_cold() updates ghosts of Enth3 using
    // update_ghosts(). Is not optimized because this
    // (do_cold_ice_methods) is a rare case.
    ierr = compute_enthalpy_cold(T3, Enth3);  CHKERRQ(ierr);

  } else {
    // new enthalpy values go in vWork3d; also updates (and communicates) Hmelt
    double myLiquifiedVol = 0.0, gLiquifiedVol;

    ierr = enthalpyAndDrainageStep(&myVertSacrCount,&myLiquifiedVol,&myBulgeCount);
    CHKERRQ(ierr);

    ierr = vWork3d.update_ghosts(Enth3); CHKERRQ(ierr);

    ierr = GlobalSum(&myLiquifiedVol, &gLiquifiedVol, grid.com); CHKERRQ(ierr);
    if (gLiquifiedVol > 0.0) {
      ierr = verbPrintf(1,grid.com,
                        "\n PISM WARNING: fully-liquified cells detected: volume liquified = %.3f km^3\n\n",
                        gLiquifiedVol / 1.0e9); CHKERRQ(ierr);
    }
  }

  ierr = GlobalSum(&myCFLviolcount, &CFLviolcount, grid.com); CHKERRQ(ierr);

  ierr = GlobalSum(&myVertSacrCount, &gVertSacrCount, grid.com); CHKERRQ(ierr);
  if (gVertSacrCount > 0.0) { // count of when BOMBPROOF switches to lower accuracy
    const double bfsacrPRCNT = 100.0 * (gVertSacrCount / (grid.Mx * grid.My));
    const double BPSACR_REPORT_VERB2_PERCENT = 5.0; // only report if above 5%
    if (bfsacrPRCNT > BPSACR_REPORT_VERB2_PERCENT &&
        getVerbosityLevel() > 2) {
      char tempstr[50] = "";
      snprintf(tempstr,50, "  [BPsacr=%.4f%%] ", bfsacrPRCNT);
      stdout_flags = tempstr + stdout_flags;
    }
  }

  ierr = GlobalSum(&myBulgeCount, &gBulgeCount, grid.com); CHKERRQ(ierr);
  if (gBulgeCount > 0.0) {   // count of when advection bulges are limited;
                             //    frequently it is identically zero
    char tempstr[50] = "";
    snprintf(tempstr,50, " BULGE=%d ", static_cast<int>(ceil(gBulgeCount)));
    stdout_flags = tempstr + stdout_flags;
  }

  return 0;
}

//! @brief Combine basal melt rate in grounded and floating areas.
/**
 * Grounded basal melt rate is computed as a part of the energy
 * (enthalpy or temperature) step; floating basal melt rate is
 * provided by an ocean model.
 *
 * This method updates IceModel::basal_melt_rate (in meters per second
 * ice-equivalent).
 *
 * The sub shelf mass flux provided by an ocean model is in [kg m-2
 * s-1], so we divide by the ice density to convert to [m/s].
 */
PetscErrorCode IceModel::combine_basal_melt_rate() {
  PetscErrorCode ierr;

  assert(ocean != NULL);
  ierr = ocean->shelf_base_mass_flux(shelfbmassflux); CHKERRQ(ierr);

  const bool sub_gl = config.get_flag("sub_groundingline");
  if (sub_gl == true) {
    ierr = gl_mask.begin_access(); CHKERRQ(ierr);
  }

  MaskQuery mask(vMask);

  double ice_density = config.get("ice_density");

  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = basal_melt_rate.begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);

  for (int i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      double lambda = 1.0;      // 1.0 corresponds to the grounded case
      // Note: here we convert shelf base mass flux from [kg m-2 s-1] to [m s-1]:
      const double
        M_grounded   = basal_melt_rate(i,j),
        M_shelf_base = shelfbmassflux(i,j) / ice_density;

      // Use the fractional floatation mask to adjust the basal melt
      // rate near the grounding line:
      if (sub_gl == true) {
        lambda = gl_mask(i,j);
      } else if (mask.ocean(i,j)) {
        lambda = 0.0;
      }
      basal_melt_rate(i,j) = lambda * M_grounded + (1.0 - lambda) * M_shelf_base;
    }
  }

  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);
  ierr = basal_melt_rate.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  if (sub_gl == true) {
    ierr = gl_mask.end_access(); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Extract from enthalpy field (Enth3) the temperature which the top of
//! the bedrock thermal layer will see.
PetscErrorCode IceModel::get_bed_top_temp(IceModelVec2S &result) {
  PetscErrorCode  ierr;
  double sea_level = 0,
    T0 = config.get("water_melting_point_temperature"),
    beta_CC_grad_sea_water = (config.get("beta_CC") * config.get("sea_water_density") *
                              config.get("standard_gravity")); // K m-1

  // will need coupler fields in ice-free land and
  assert(surface != NULL);
  ierr = surface->ice_surface_temperature(ice_surface_temp); CHKERRQ(ierr);

  assert(ocean != NULL);
  ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);

  // start by grabbing 2D basal enthalpy field at z=0; converted to temperature if needed, below
  ierr = Enth3.getHorSlice(result, 0.0); CHKERRQ(ierr);

  MaskQuery mask(vMask);

  ierr = bed_topography.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = ice_surface_temp.begin_access(); CHKERRQ(ierr);
  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (mask.grounded(i,j)) {
        if (mask.ice_free(i,j)) { // no ice: sees air temp
          result(i,j) = ice_surface_temp(i,j);
        } else { // ice: sees temp of base of ice
          const double pressure = EC->getPressureFromDepth(ice_thickness(i,j));
          double temp;
          // ignore return code when getting temperature: we are committed to
          //   this enthalpy field; getAbsTemp() only returns temperatures at or
          //   below pressure melting
          EC->getAbsTemp(result(i,j), pressure, temp);
          result(i,j) = temp;
        }
      } else { // floating: apply pressure melting temp as top of bedrock temp
        result(i,j) = T0 - (sea_level - bed_topography(i,j)) * beta_CC_grad_sea_water;
      }
    }
  }
  ierr = ice_thickness.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = ice_surface_temp.end_access(); CHKERRQ(ierr);
  ierr = bed_topography.end_access(); CHKERRQ(ierr);

  return 0;
}


//! \brief Is one of my neighbors below a critical thickness to apply advection
//! in enthalpy or temperature equation?
bool IceModel::checkThinNeigh(IceModelVec2S &thickness, int i, int j, const double threshold) {
  const double
    N  = thickness(i, j + 1),
    E  = thickness(i + 1, j),
    S  = thickness(i, j - 1),
    W  = thickness(i - 1, j),
    NW = thickness(i - 1, j + 1),
    SW = thickness(i - 1, j - 1),
    NE = thickness(i + 1, j + 1),
    SE = thickness(i + 1, j - 1);

  return ((E < threshold) || (NE < threshold) || (N < threshold) || (NW < threshold) ||
          (W < threshold) || (SW < threshold) || (S < threshold) || (SE < threshold));
}

} // end of namespace pism

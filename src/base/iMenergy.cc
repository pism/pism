// Copyright (C) 2004-2011, 2013 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "iceModel.hh"
#include "Mask.hh"
#include "bedrockThermalUnit.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "enthalpyConverter.hh"

//! \file iMenergy.cc Methods of IceModel which address conservation of energy.
//! Common to enthalpy (polythermal) and temperature (cold-ice) methods.

//! Manage the solution of the energy equation, and related parallel communication.
/*!
This method updates three fields:
  - IceModelVec3 Enth3
  - IceModelVec2 vbmr
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

  PetscScalar  myCFLviolcount = 0.0,   // these are counts but they are type "PetscScalar"
               myVertSacrCount = 0.0,  //   because that type works with PISMGlobalSum()
               myBulgeCount = 0.0;
  PetscScalar gVertSacrCount, gBulgeCount;

  // always count CFL violations for sanity check (but can occur only if -skip N with N>1)
  ierr = countCFLViolations(&myCFLviolcount); CHKERRQ(ierr);

  // operator-splitting occurs here (ice and bedrock energy updates are split):
  //   tell PISMBedThermalUnit* btu that we have an ice base temp; it will return
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
    PetscScalar myLiquifiedVol = 0.0, gLiquifiedVol;

    ierr = enthalpyAndDrainageStep(&myVertSacrCount,&myLiquifiedVol,&myBulgeCount);
    CHKERRQ(ierr);

    ierr = vWork3d.update_ghosts(Enth3); CHKERRQ(ierr);

    ierr = PISMGlobalSum(&myLiquifiedVol, &gLiquifiedVol, grid.com); CHKERRQ(ierr);
    if (gLiquifiedVol > 0.0) {
      ierr = verbPrintf(1,grid.com,
        "\n PISM WARNING: fully-liquified cells detected: volume liquified = %.3f km^3\n\n",
        gLiquifiedVol / 1.0e9); CHKERRQ(ierr);
    }
  }

  // Both cases above update the basal melt rate field; here we update its
  // ghosts, which are needed to compute tauc locally
  ierr = vbmr.update_ghosts(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&myCFLviolcount, &CFLviolcount, grid.com); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&myVertSacrCount, &gVertSacrCount, grid.com); CHKERRQ(ierr);
  if (gVertSacrCount > 0.0) { // count of when BOMBPROOF switches to lower accuracy
    const PetscScalar bfsacrPRCNT = 100.0 * (gVertSacrCount / (grid.Mx * grid.My));
    const PetscScalar BPSACR_REPORT_VERB2_PERCENT = 5.0; // only report if above 5%
    if (   (bfsacrPRCNT > BPSACR_REPORT_VERB2_PERCENT) 
        && (getVerbosityLevel() > 2)                    ) {
      char tempstr[50] = "";
      snprintf(tempstr,50, "  [BPsacr=%.4f%%] ", bfsacrPRCNT);
      stdout_flags = tempstr + stdout_flags;
    }
  }

  ierr = PISMGlobalSum(&myBulgeCount, &gBulgeCount, grid.com); CHKERRQ(ierr);
  if (gBulgeCount > 0.0) {   // count of when advection bulges are limited;
                             //    frequently it is identically zero
    char tempstr[50] = "";
    snprintf(tempstr,50, " BULGE=%d ", static_cast<int>(ceil(gBulgeCount)) );
    stdout_flags = tempstr + stdout_flags;
  }

  return 0;
}


//! \brief Extract from enthalpy field (Enth3) the temperature which the top of
//! the bedrock thermal layer will see.
PetscErrorCode IceModel::get_bed_top_temp(IceModelVec2S &result) {
  PetscErrorCode  ierr;
  PetscReal sea_level = 0,
    T0 = config.get("water_melting_point_temperature"),
    beta_CC_grad_sea_water = (config.get("beta_CC") * config.get("sea_water_density") *
                              config.get("standard_gravity")); // K m-1

  // will need coupler fields in ice-free land and 
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_temperature(artm); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 1,"PISM ERROR: surface == PETSC_NULL");
  }
  if (ocean != PETSC_NULL) {
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 5,"PISM ERROR: ocean == PETSC_NULL");
  }

  // start by grabbing 2D basal enthalpy field at z=0; converted to temperature if needed, below
  ierr = Enth3.getHorSlice(result, 0.0); CHKERRQ(ierr);

  MaskQuery mask(vMask);

  ierr = vbed.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = artm.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (mask.grounded(i,j)) {
        if (mask.ice_free(i,j)) { // no ice: sees air temp
          result(i,j) = artm(i,j);
        } else { // ice: sees temp of base of ice
          const PetscReal pressure = EC->getPressureFromDepth(vH(i,j));
          PetscReal temp;
          // ignore return code when getting temperature: we are committed to
          //   this enthalpy field; getAbsTemp() only returns temperatures at or
          //   below pressure melting
          EC->getAbsTemp(result(i,j), pressure, temp);
          result(i,j) = temp;
        }
      } else { // floating: apply pressure melting temp as top of bedrock temp
        result(i,j) = T0 - (sea_level - vbed(i,j)) * beta_CC_grad_sea_water;
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = artm.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);

  return 0;
}


//! \brief Is one of my neighbors below a critical thickness to apply advection
//! in enthalpy or temperature equation?
bool IceModel::checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                              PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE) {
  // FIXME: silly hard-wired critical level, but we want to avoid config.get() in loops.
  const PetscScalar THIN = 100.0;  // thin = (at most 100m thick)
  return (   (E < THIN) || (NE < THIN) || (N < THIN) || (NW < THIN)
          || (W < THIN) || (SW < THIN) || (S < THIN) || (SE < THIN) );
}


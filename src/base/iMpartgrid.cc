// Copyright (C) 2011 Torsten Albrecht and Ed Bueler
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


// methods implementing PIK logic for -part_grid; see Albrecht et al 2011


//! Compute staggered grid velocities according to mask and regular grid velocities.
/*!
  In the finite volume interpretation, these are normal velocities at the faces
  of the cell.  The method avoids differencing velocities from ice free ocean locations.

  The outputs are vel_E, vel_W, vel_N, vel_S.
*/
PetscErrorCode IceModel::velsPartGrid(int M_ij, 
                                      int M_e, int M_w, 
                                      int M_n, int M_s, 
                                      PISMVector2 vreg_ij, 
                                      PISMVector2 vreg_e, PISMVector2 vreg_w, 
                                      PISMVector2 vreg_n, PISMVector2 vreg_s, 
                                      PetscReal &vel_E, PetscReal &vel_W, 
                                      PetscReal &vel_N, PetscReal &vel_S) {
  const bool oneneighboricefree = (M_e > MASK_FLOATING ||
                                   M_w > MASK_FLOATING ||
                                   M_n > MASK_FLOATING ||
                                   M_s > MASK_FLOATING), 
    oneneighboricefilled = (M_e <= MASK_FLOATING ||
                            M_w <= MASK_FLOATING ||
                            M_n <= MASK_FLOATING ||
                            M_s <= MASK_FLOATING);

  //case1: [i][j] in the middle of ice or bedrock: default scheme
  if (M_ij <= MASK_FLOATING && (!oneneighboricefree)) {
    // compute (i, j) - centered "face" velocity components by average
    vel_E = 0.5 * (vreg_ij.u + vreg_e.u);
    vel_W = 0.5 * (vreg_w.u + vreg_ij.u);
    vel_N = 0.5 * (vreg_ij.v + vreg_n.v);
    vel_S = 0.5 * (vreg_s.v + vreg_ij.v);
    //case2: [i][j] on floating or grounded ice, but next to a ice - free ocean grid cell
  } else if (M_ij <= MASK_FLOATING && (oneneighboricefree)) {
    vel_E = (M_e > MASK_FLOATING ? vreg_ij.u : 0.5 * (vreg_ij.u + vreg_e.u));
    vel_W = (M_w > MASK_FLOATING ? vreg_ij.u : 0.5 * (vreg_w.u + vreg_ij.u));
    vel_N = (M_n > MASK_FLOATING ? vreg_ij.v : 0.5 * (vreg_ij.v + vreg_n.v));
    vel_S = (M_s > MASK_FLOATING ? vreg_ij.v : 0.5 * (vreg_s.v + vreg_ij.v));
    //case3: [i][j] on ice - free ocean (or partially filled), but next to ice grid cell
  } else if (M_ij > MASK_FLOATING && (oneneighboricefilled)) {
    vel_E = (M_e <= MASK_FLOATING ? vreg_e.u : 0.0);
    vel_W = (M_w <= MASK_FLOATING ? vreg_w.u : 0.0);
    vel_N = (M_n <= MASK_FLOATING ? vreg_n.v : 0.0);
    vel_S = (M_s <= MASK_FLOATING ? vreg_s.v : 0.0);
    //case4: [i][j] on ice - free ocean, and no ice neighbors, and else
  } else {
    vel_E = 0.0;
    vel_W = 0.0;
    vel_N = 0.0;
    vel_S = 0.0;
  }

  return 0;
}


//! For ice-free (or partially-filled) cells adjacent to "full" floating cells, update Href.
/*!
  Should only be called if one of the neighbors is floating, i.e. only if at
  least one of M_e, M_w, M_n, M_s is MASK_FLOATING.

  FIXME: does not account for grounded tributaries: thin ice shelves may
  evolve from grounded tongue
*/
PetscReal IceModel::getHav(bool do_redist, 
                           int M_e, int M_w, int M_n, int M_s, 
                           PetscReal H_e, PetscReal H_w, PetscReal H_n, PetscReal H_s) {

  // get mean ice thickness over adjacent floating ice shelf neighbors
  PetscReal H_average = 0.0;
  PetscInt N = 0;
  if (M_e == MASK_FLOATING) { H_average += H_e; N++; }
  if (M_w == MASK_FLOATING) { H_average += H_w; N++; }
  if (M_n == MASK_FLOATING) { H_average += H_n; N++; }
  if (M_s == MASK_FLOATING) { H_average += H_s; N++; }

  if (N == 0) {
    SETERRQ(1, "N == 0;  call me only if a neighbor is floating!\n");
  }

  H_average = H_average / N;

  // reduces the guess at the front
  if (do_redist) {
    const PetscReal  mslope = 2.4511e-18*grid.dx / (300*600 / secpera);
    // for declining front C / Q0 according to analytical flowline profile in
    //   vandeveen with v0 = 300m / yr and H0 = 600m
    H_average -= 0.8*mslope*pow(H_average, 5);
  }

  return H_average;
}


//! Redistribute residual ice mass from subgrid-scale parameterization, when using -part_redist option.
/*!
  See [\ref Albrechtetal2011].  Manages the loop.

  FIXME: Reporting!

  FIXME: repeatRedist should be config flag?

  FIXME: resolve fixed number (=3) of loops issue
*/
PetscErrorCode IceModel::redistResiduals() {
  PetscErrorCode ierr;
  const PetscInt max_loopcount = 3;
  ierr = calculateRedistResiduals(); CHKERRQ(ierr); //while loop?

  for (int i = 0; i < max_loopcount && repeatRedist == PETSC_TRUE; ++i) {
    ierr = calculateRedistResiduals(); CHKERRQ(ierr); // sets repeatRedist
    ierr = verbPrintf(4, grid.com, "redistribution loopcount = %d\n", i); CHKERRQ(ierr);
  }
  return 0;
}


// This routine carries-over the ice mass when using -part_redist option, one step in the loop.
PetscErrorCode IceModel::calculateRedistResiduals() {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com, "calculateRedistResiduals() is called\n"); CHKERRQ(ierr);

  IceModelVec2S vHnew = vWork2d[0];
  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);

  ierr = vHref.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);

  IceModelVec2S vHresidualnew = vWork2d[1];
  ierr = vHresidual.copy_to(vHresidualnew); CHKERRQ(ierr);
  ierr = vHresidual.begin_access(); CHKERRQ(ierr);
  ierr = vHresidualnew.begin_access(); CHKERRQ(ierr);


  if (ocean == PETSC_NULL) { SETERRQ(1, "PISM ERROR: ocean == PETSC_NULL");  }
  PetscReal sea_level = 0.0; //FIXME
  ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);

  PetscScalar minHRedist = 0.0; // to avoid the propagation of thin ice shelf tongues

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      // first step: distributing residual ice masses
      if (vHresidual(i, j) > 0.0) {

        PetscInt N = 0; // counting empty / partially filled neighbors
        PetscTruth exEast = PETSC_FALSE, exWest = PETSC_FALSE, exNorth = PETSC_FALSE, exSouth = PETSC_FALSE;

        // check for partially filled / empty grid cell neighbors (mask not updated yet, but vH is)
        if (vH(i + 1, j) == 0.0 && vbed(i + 1, j) < sea_level) {N += 1; exEast = PETSC_TRUE;}
        if (vH(i - 1, j) == 0.0 && vbed(i - 1, j) < sea_level) {N += 1; exWest = PETSC_TRUE;}
        if (vH(i, j + 1) == 0.0 && vbed(i, j + 1) < sea_level) {N += 1; exNorth = PETSC_TRUE;}
        if (vH(i, j - 1) == 0.0 && vbed(i, j - 1) < sea_level) {N += 1; exSouth = PETSC_TRUE;}

        if (N > 0 && vH(i, j) > minHRedist)  {
          //remainder ice mass will be redistributed equally to all adjacent imfrac boxes (is there a more physical way?)
          if (exEast) vHref(i + 1, j) += vHresidual(i, j) / N;
          if (exWest) vHref(i - 1, j) += vHresidual(i, j) / N;
          if (exNorth) vHref(i, j + 1) += vHresidual(i, j) / N;
          if (exSouth) vHref(i, j - 1) += vHresidual(i, j) / N;

          vHresidualnew(i, j) = 0.0;
        } else {
          vHnew(i, j) += vHresidual(i, j); // mass conservation, but thick ice at one grid cell possible
          vHresidualnew(i, j) = 0.0;
          ierr = verbPrintf(4, grid.com, 
                            "!!! PISM WARNING: Hresidual has %d partially filled neighbors, "
                            " set ice thickness to vHnew = %.2e at %d, %d \n", 
                            N, vHnew(i, j), i, j ); CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);


  double  ocean_rho = config.get("sea_water_density");
  double  ice_rho = config.get("ice_density");
  PetscScalar     H_average;
  PetscScalar     Hcut = 0.0;
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      // second step: if neighbors which gained redistributed ice also become
      // full, this needs to be redistributed in a repeated loop
      if (vHref(i, j) > 0.0) {
        H_average = 0.0;
        PetscInt N = 0; // number of full floating ice neighbors (mask not yet updated)

        if (vH(i + 1, j) > 0.0 && (vbed(i + 1, j) < (sea_level - ice_rho / ocean_rho*vH(i + 1, j))))
          { H_average += vH(i + 1, j); N += 1;}
        if (vH(i - 1, j) > 0.0 && (vbed(i - 1, j) < (sea_level - ice_rho / ocean_rho*vH(i - 1, j))))
          { H_average += vH(i - 1, j); N += 1;}
        if (vH(i, j + 1) > 0.0 && (vbed(i, j + 1) < (sea_level - ice_rho / ocean_rho*vH(i, j + 1))))
          { H_average += vH(i, j + 1); N += 1;}
        if (vH(i, j - 1) > 0.0 && (vbed(i, j - 1) < (sea_level - ice_rho / ocean_rho*vH(i, j - 1))))
          { H_average += vH(i, j - 1); N += 1;}

        if (N > 0){
          H_average = H_average / N;

          PetscScalar coverageRatio = vHref(i, j) / H_average;
          if (coverageRatio > 1.0) { // partially filled grid cell is considered to be full
            vHresidualnew(i, j) = vHref(i, j) - H_average;
            Hcut += vHresidualnew(i, j); // summed up to decide, if methods needs to be run once more
            vHnew(i, j) += H_average;
            vHref(i, j) = 0.0;
          }
        } else { // no full floating ice neighbor
          vHnew(i, j) += vHref(i, j); // mass conservation, but thick ice at one grid cell possible
          vHref(i, j) = 0.0;
          vHresidualnew(i, j) = 0.0;
          ierr = verbPrintf(4, grid.com, 
                            "!!! PISM_WARNING: No floating ice neighbors to calculate H_average, "
                            " set ice thickness to vHnew = %.2e at %d,%d \n", 
                            vHnew(i, j), i, j); CHKERRQ(ierr);
        }
      }
    }
  }

  PetscScalar gHcut; //check, if redistribution should be run once more
  ierr = PetscGlobalSum(&Hcut, &gHcut, grid.com); CHKERRQ(ierr);
  if (gHcut > 0.0) { repeatRedist = PETSC_TRUE;}
  else { repeatRedist = PETSC_FALSE;}
  //ierr = verbPrintf(3, grid.com, "!!! Hcut = %f \n", gHcut); CHKERRQ(ierr);

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vHnew.end_access(); CHKERRQ(ierr);
  // finally copy vHnew into vH and communicate ghosted values
  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);

  ierr = vHref.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);

  ierr = vHresidual.end_access(); CHKERRQ(ierr);
  ierr = vHresidualnew.end_access(); CHKERRQ(ierr);
  ierr = vHresidualnew.beginGhostComm(vHresidual); CHKERRQ(ierr);
  ierr = vHresidualnew.endGhostComm(vHresidual); CHKERRQ(ierr);


  return 0;
}

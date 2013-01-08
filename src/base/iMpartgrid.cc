// Copyright (C) 2011, 2012, 2013 Torsten Albrecht and Ed Bueler and Constantine Khroulev
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
#include <petscdmda.h>

#include "iceModel.hh"
#include "Mask.hh"
#include "PISMStressBalance.hh"
#include "PISMOcean.hh"


//! \file iMpartgrid.cc Methods implementing PIK option -part_grid [\ref Albrechtetal2011].

//! For ice-free (or partially-filled) cells adjacent to "full" floating cells, update Href.
/*!
  Should only be called if one of the neighbors is floating.

  FIXME: does not account for grounded tributaries: thin ice shelves may
  evolve from grounded tongue
*/
PetscReal IceModel::get_average_thickness(bool do_redist, planeStar<int> M, planeStar<PetscScalar> H) {

  // get mean ice thickness over adjacent floating ice shelf neighbors
  PetscReal H_average = 0.0;
  PetscInt N = 0;
  Mask m;

  if (m.floating_ice(M.e)) { H_average += H.e; N++; }
  if (m.floating_ice(M.w)) { H_average += H.w; N++; }
  if (m.floating_ice(M.n)) { H_average += H.n; N++; }
  if (m.floating_ice(M.s)) { H_average += H.s; N++; }

  if (N == 0) {
    SETERRQ(grid.com, 1, "N == 0;  call this only if a neighbor is floating!\n");
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

  IceModelVec2S vHnew = vWork2d[0];
  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);

  IceModelVec2S vHresidualnew = vWork2d[1];
  ierr = vHresidual.copy_to(vHresidualnew); CHKERRQ(ierr);

  if (ocean == PETSC_NULL) { SETERRQ(grid.com, 1, "PISM ERROR: ocean == PETSC_NULL");  }
  PetscReal sea_level = 0.0; //FIXME
  ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);

  PetscScalar minHRedist = 50.0; // to avoid the propagation of thin ice shelf tongues

  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHref.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  ierr = vHresidual.begin_access(); CHKERRQ(ierr);
  ierr = vHresidualnew.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      // first step: distributing residual ice masses
      if (vHresidual(i, j) > 0.0 && putOnTop==PETSC_FALSE) {
        
        planeStar<PetscScalar> thk = vH.star(i, j),
          bed = vbed.star(i, j);

        PetscInt N = 0; // counting empty / partially filled neighbors
        planeStar<bool> neighbors;
        neighbors.e = neighbors.w = neighbors.n = neighbors.s = false;

        // check for partially filled / empty grid cell neighbors (mask not updated yet, but vH is)
        if (thk.e == 0.0 && bed.e < sea_level) {N++; neighbors.e = true;}
        if (thk.w == 0.0 && bed.w < sea_level) {N++; neighbors.w = true;}
        if (thk.n == 0.0 && bed.n < sea_level) {N++; neighbors.n = true;}
        if (thk.s == 0.0 && bed.s < sea_level) {N++; neighbors.s = true;}

        if (N > 0)  {
          //remainder ice mass will be redistributed equally to all adjacent
          //imfrac boxes (is there a more physical way?)
          if (neighbors.e) vHref(i + 1, j) += vHresidual(i, j) / N;
          if (neighbors.w) vHref(i - 1, j) += vHresidual(i, j) / N;
          if (neighbors.n) vHref(i, j + 1) += vHresidual(i, j) / N;
          if (neighbors.s) vHref(i, j - 1) += vHresidual(i, j) / N;
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
  ierr = vHnew.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = vHnew.update_ghosts(vH); CHKERRQ(ierr);

  double  ocean_rho = config.get("sea_water_density"),
    ice_rho = config.get("ice_density"),
    C = ice_rho / ocean_rho;
  PetscScalar     H_average;
  PetscScalar     Hcut = 0.0;

  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      // second step: if neighbors which gained redistributed ice also become
      // full, this needs to be redistributed in a repeated loop
      if (vHref(i, j) > 0.0) {
        H_average = 0.0;
        PetscInt N = 0; // number of full floating ice neighbors (mask not yet updated)

        planeStar<PetscScalar> thk = vH.star(i, j),
          bed = vbed.star(i, j);

        if (thk.e > 0.0 && bed.e < sea_level - C * thk.e) { N++; H_average += thk.e; }
        if (thk.w > 0.0 && bed.w < sea_level - C * thk.w) { N++; H_average += thk.w; }
        if (thk.n > 0.0 && bed.n < sea_level - C * thk.n) { N++; H_average += thk.n; }
        if (thk.s > 0.0 && bed.s < sea_level - C * thk.s) { N++; H_average += thk.s; }

        if (N > 0){
          H_average = H_average / N;

          PetscScalar coverageRatio = vHref(i, j) / H_average;
          if (coverageRatio > 1.0) { // partially filled grid cell is considered to be full
            vHresidualnew(i, j) = vHref(i, j) - H_average;
            Hcut += vHresidualnew(i, j); // summed up to decide, if methods needs to be run once more
            vHnew(i, j) += H_average; //SMB?
            vHref(i, j) = 0.0;
          }
        } else { // no full floating ice neighbor
          vHnew(i, j) += vHref(i, j); // mass conservation, but thick ice at one grid cell possible
          vHref(i, j) = 0.0;
          vHresidualnew(i, j) = 0.0;
	      ierr = verbPrintf(4, grid.com, 
                            "!!! PISM WARNING: Hresidual=%.2f with %d partially filled neighbors, "
                            " set ice thickness to vHnew = %.2f at %d, %d \n", 
                            vHresidual(i, j), N , vHnew(i, j), i, j ); CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vHnew.end_access(); CHKERRQ(ierr);
  ierr = vHref.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  ierr = vHresidual.end_access(); CHKERRQ(ierr);
  ierr = vHresidualnew.end_access(); CHKERRQ(ierr);
  
  PetscScalar gHcut; //check, if redistribution should be run once more
  ierr = PISMGlobalSum(&Hcut, &gHcut, grid.com); CHKERRQ(ierr);
  putOnTop = PETSC_FALSE;
  if (gHcut > 0.0) {
    repeatRedist = PETSC_TRUE;
    // avoid repetition for the redistribution of very thin vHresiduals
    if (gHcut < minHRedist) { putOnTop = PETSC_TRUE; }
  } else {
    repeatRedist = PETSC_FALSE;
  }

  // finally copy vHnew into vH and communicate ghosted values
  ierr = vHnew.update_ghosts(vH); CHKERRQ(ierr);
  ierr = vHresidualnew.update_ghosts(vHresidual); CHKERRQ(ierr);

  return 0;
}


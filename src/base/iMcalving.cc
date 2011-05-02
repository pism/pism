// Copyright (C) 2004--2011 Torsten Albrecht
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
#include <petscda.h>
#include "iceModel.hh"
#include "pism_signal.h"
#include "Mask.hh"

//! \brief Uses principal strain rates to apply "eigencalving" with constant K.
/*!
  See equation (26) in [\ref Winkelmannetal2010TCD].
*/
PetscErrorCode IceModel::eigenCalving() {
  const PetscScalar   dx = grid.dx, dy = grid.dy;
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com, "######### applyCalvingRate is called \n");    CHKERRQ(ierr);

  // is ghost communication really needed here?
  ierr = vH.beginGhostComm(); CHKERRQ(ierr);
  ierr = vH.endGhostComm(); CHKERRQ(ierr);
  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);

  double ocean_rho = config.get("sea_water_density");
  double ice_rho = config.get("ice_density");

  const PetscScalar eigenCalvFactor = config.get("eigen_calving");

  PetscReal sea_level = 0;

  if (ocean != NULL) {
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
  } else { SETERRQ(2, "PISM ERROR: ocean == NULL"); }

  IceModelVec2S vHnew = vWork2d[0];
  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);

  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  ierr = vHref.begin_access(); CHKERRQ(ierr);
  ierr = vPrinStrain1.begin_access(); CHKERRQ(ierr);
  ierr = vPrinStrain2.begin_access(); CHKERRQ(ierr);

  IceModelVec2S vDiffCalvRate = vWork2d[1];
  ierr = vDiffCalvRate.set(0.0); CHKERRQ(ierr);
  ierr = vDiffCalvRate.begin_access(); CHKERRQ(ierr);

  if(dx != dy) {
    ierr = PetscPrintf(grid.com,
                       "PISMPIK_ERROR: Calvingrate using a non-quadratic grid does not work (yet), "
                       " since it has no direction!!!\n");
    PISMEnd();
  }

  // Distance from calving front where straine rate is evaluated
  PetscInt offset = 2;

  MaskQuery mask(vMask);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      // Average of strain-rate eigenvalues in adjacent floating gird cells to
      // be used for eigencalving
      PetscScalar eigen1 = 0.0, eigen2 = 0.0;

      // Number of directly adjacent floating boxes
      PetscInt N = 0;

      // Neighbor-averaged ice thickness
      PetscScalar H_average = 0.0; // is calculated here as average over direct neighbors

      // Counting adjacent floating boxes (with distance "offset")
      PetscInt M = 0;

      // find partially filled or empty grid boxes on the icefree ocean, which
      // have floating ice neighbors after massContExplicitStep (mask not
      // updated)
      bool next_to_floating =
        ((vH(i + 1, j) > 0.0 && (vbed(i + 1, j) < (sea_level - ice_rho / ocean_rho*vH(i + 1, j)))) ||
         (vH(i - 1, j) > 0.0 && (vbed(i - 1, j) < (sea_level - ice_rho / ocean_rho*vH(i - 1, j)))) ||
         (vH(i, j + 1) > 0.0 && (vbed(i, j + 1) < (sea_level - ice_rho / ocean_rho*vH(i, j + 1)))) ||
         (vH(i, j - 1) > 0.0 && (vbed(i, j - 1) < (sea_level - ice_rho / ocean_rho*vH(i, j - 1)))));

      bool ice_free_ocean = ( vH(i, j) == 0.0 && vbed(i, j) < sea_level );

      H_average = 0.0; eigen1 = 0.0; eigen2 = 0.0; M = 0, N = 0;

      if (ice_free_ocean && next_to_floating) {

        if ( mask.floating_ice(i + 1, j) ) { N += 1; H_average += vH(i + 1, j); }
        if ( mask.floating_ice(i - 1, j) ) { N += 1; H_average += vH(i - 1, j); }
        if ( mask.floating_ice(i, j + 1) ) { N += 1; H_average += vH(i, j + 1); }
        if ( mask.floating_ice(i, j - 1) ) { N += 1; H_average += vH(i, j - 1); }
        if (N > 0)
          H_average /= N;

        if ( mask.floating_ice(i + offset, j) ){
          eigen1 += vPrinStrain1(i + offset, j);
          eigen2 += vPrinStrain2(i + offset, j);
          M += 1;
        }

        if ( mask.floating_ice(i - offset, j) ){
          eigen1 += vPrinStrain1(i - offset, j);
          eigen2 += vPrinStrain2(i - offset, j);
          M += 1;
        }

        if ( mask.floating_ice(i, j + offset) ){
          eigen1 += vPrinStrain1(i, j + offset);
          eigen2 += vPrinStrain2(i, j + offset);
          M += 1;
        }

        if ( mask.floating_ice(i, j - offset) ){
          eigen1 += vPrinStrain1(i, j - offset);
          eigen2 += vPrinStrain2(i, j - offset);
          M += 1;
        }

        if (M > 0) {
          eigen1 /= M;
          eigen2 /= M;
        }


        PetscScalar calvrateHorizontal = 0.0,
          eigenCalvOffset = 0.0; // if it's not exactly the zero line of
                               // transition from compressive to extensive flow
                               // regime

        // calving law
        if ( eigen2 > eigenCalvOffset && eigen1 > 0.0) { // if spreading in all directions
          calvrateHorizontal = eigenCalvFactor * eigen1 * (eigen2 - eigenCalvOffset);
          // eigen1 * eigen2 has units [s^ - 2] and calvrateHorizontal [m*s^1]
          // hence, eigenCalvFactor has units [m*s]
        } else calvrateHorizontal = 0.0;

        // calculate mass loss with respect to the associated ice thickness and the grid size:
        PetscScalar calvrate = calvrateHorizontal * H_average / dx; // in m/s

        // apply calving rate at partially filled or empty grid cells
        if (calvrate > 0.0) {
          //PetscScalar Href_old = vHref(i, j);
          //vDiffCalvRate(i, j) = 0.0;
          vHref(i, j) -= calvrate * dt; // in m
          if(vHref(i, j) < 0.0) { // i.e. partially filled grid cell has completely calved off
            vDiffCalvRate(i, j) =  - vHref(i, j) / dt;// in m/s, means additional ice loss
          }
          vHref(i, j) = 0.0;
          if(N > 0){
            vDiffCalvRate(i, j) = vDiffCalvRate(i, j) / N;
          }
        }
      }
    }
  }
  ierr = vDiffCalvRate.end_access(); CHKERRQ(ierr);

  ierr = vDiffCalvRate.beginGhostComm(); CHKERRQ(ierr);
  ierr = vDiffCalvRate.endGhostComm(); CHKERRQ(ierr);

  ierr = vDiffCalvRate.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      PetscScalar restCalvRate = 0.0;
      bool hereFloating = (vH(i, j) > 0.0 && (vbed(i, j) < (sea_level - ice_rho / ocean_rho*vH(i, j))));

      if (hereFloating &&
          (vDiffCalvRate(i + 1, j) > 0.0 || vDiffCalvRate(i - 1, j) > 0.0 ||
           vDiffCalvRate(i, j + 1) > 0.0 || vDiffCalvRate(i, j - 1) > 0.0 )) {

        restCalvRate = (vDiffCalvRate(i + 1, j) +
                        vDiffCalvRate(i - 1, j) +
                        vDiffCalvRate(i, j + 1) +
                        vDiffCalvRate(i, j - 1));     // in m/s

        vHref(i, j) = vH(i, j) - (restCalvRate * dt); // in m

        //Hav = vH(i, j);
        vHnew(i, j) = 0.0;

        if(vHref(i, j) < 0.0) { // i.e. terminal floating ice grid cell has calved off completely.
          // We do not account for further calving ice - inwards!
          // Alternatively CFL criterion for time stepping could be adjusted to maximum of calving rate.
          //Hav = 0.0;
          vHref(i, j) = 0.0;
        }
      }
    }
  }
  ierr = vHnew.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);

  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  ierr = vHref.end_access(); CHKERRQ(ierr);
  ierr = vPrinStrain1.end_access(); CHKERRQ(ierr);
  ierr = vPrinStrain2.end_access(); CHKERRQ(ierr);
  ierr = vDiffCalvRate.end_access(); CHKERRQ(ierr);

  return 0;
}


/*!
This calving condition applies for terminal floating ice shelf grid cells when their thickness is less than a threshold
*/
PetscErrorCode IceModel::calvingAtThickness() {
  //const PetscScalar   dx = grid.dx, dy = grid.dy;
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com, "######### callvingAtThickness is called \n");    CHKERRQ(ierr);

  ierr = vH.beginGhostComm(); CHKERRQ(ierr);
  ierr = vH.endGhostComm(); CHKERRQ(ierr);

  IceModelVec2S vHnew = vWork2d[0];
  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);

  //ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);

  double ocean_rho = config.get("sea_water_density"),
    ice_rho = config.get("ice_density");

  PetscReal sea_level;
  if (ocean != NULL) {
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
  } else { SETERRQ(2, "PISM ERROR: ocean == NULL"); }

  const PetscScalar Hcalving = config.get("calving_at_thickness");
  //ierr = verbPrintf(3, grid.com, "!!! Hcalving=%f \n", Hcalving);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      bool hereFloating = (vH(i, j) > 0.0 && (vbed(i, j) < (sea_level - ice_rho / ocean_rho * vH(i, j))));
      bool icefreeOceanNeighbor = ( (vH(i + 1, j) == 0.0 && vbed(i + 1, j) < sea_level) ||
                                    (vH(i - 1, j) == 0.0 && vbed(i - 1, j) < sea_level) ||
                                    (vH(i, j + 1) == 0.0 && vbed(i, j + 1) < sea_level) ||
                                    (vH(i, j - 1) == 0.0 && vbed(i, j - 1) < sea_level));


      if (hereFloating && vH(i, j) <= Hcalving && icefreeOceanNeighbor) {
        vHnew(i, j) = 0.0;
      }
    }
  }
  ierr = vHnew.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);

  //ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);

  return 0;
}


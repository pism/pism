// Copyright (C) 2004--2013 Torsten Albrecht and Constantine Khroulev
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


#include <cmath>
#include <petscdmda.h>
#include "iceModel.hh"
#include "pism_signal.h"
#include "Mask.hh"
#include "PISMOcean.hh"
#include "PISMOceanKill.hh"
#include "PISMCalvingAtThickness.hh"
#include "PISMStressBalance.hh"
#include "PISMIcebergRemover.hh"

PetscErrorCode IceModel::do_calving() {
  PetscErrorCode ierr;

  if (config.get_flag("do_eigen_calving") && config.get_flag("use_ssa_velocity")) {
    bool dteigencalving = config.get_flag("cfl_eigencalving");
    if (dteigencalving == false){
      // otherwise calculation of strain rates has been done in iMadaptive.cc already
      IceModelVec2V *ssa_velocity;
      ierr = stress_balance->get_2D_advective_velocity(ssa_velocity); CHKERRQ(ierr);
      ierr = stress_balance->compute_2D_principal_strain_rates(*ssa_velocity,
                                                               vMask, strain_rates); CHKERRQ(ierr);
    }
    ierr = eigenCalving(); CHKERRQ(ierr);
  }
  
  if (ocean_kill_calving != NULL) {
    ierr = ocean_kill_calving->update(vMask, vH); CHKERRQ(ierr);
  }

  if (thickness_threshold_calving != NULL) {
    ierr = thickness_threshold_calving->update(vMask, vH); CHKERRQ(ierr);
  }

  if (config.get_flag("kill_icebergs") || iceberg_remover != NULL) {
    ierr = iceberg_remover->update(vMask, vH); CHKERRQ(ierr);
    // the call above modifies ice thickness and updates the mask
    // accordingly
    ierr = update_surface_elevation(vbed, vH, vh); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Uses principal strain rates to apply "eigencalving" with constant K.
/*!
  See equation (26) in [\ref Winkelmannetal2011].
*/
PetscErrorCode IceModel::eigenCalving() {
  PetscErrorCode ierr;
  const PetscScalar eigen_calving_K = config.get("eigen_calving_K");

  // Distance (grid cells) from calving front where strain rate is evaluated
  // Note: the strain_rates field has to have stencil width of at least 'offset'.
  PetscInt offset = 2;

  IceModelVec2S thk_loss = vWork2d[0];
  ierr = thk_loss.set(0.0); CHKERRQ(ierr);

  if (PetscAbs(grid.dx - grid.dy) / PetscMin(grid.dx, grid.dy) > 1e-2) {
    PetscPrintf(grid.com,
                "PISM ERROR: -eigen_calving using a non-square grid cell is not implemented (yet);\n"
                "            dx = %f, dy = %f, rel. diff = %f",
                grid.dx, grid.dy, PetscAbs(grid.dx - grid.dy) / PetscMax(grid.dx, grid.dy));
    PISMEnd();
  }

  MaskQuery mask(vMask);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vHref.begin_access(); CHKERRQ(ierr);
  ierr = strain_rates.begin_access(); CHKERRQ(ierr);
  ierr = thk_loss.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      PetscScalar
        eigen1    = 0.0,        // Average of strain-rate eigenvalues
                                // in adjacent floating grid cells to
                                // be used for eigencalving
        eigen2    = 0.0,
        H_average = 0.0;        // Neighbor-averaged ice thickness


      PetscInt N_floating_neighbors = 0, M = 0;
      // M is the number of adjacent floating cells (with distance
      // "offset")

      // find partially filled or empty grid boxes on the icefree ocean, which
      // have floating ice neighbors after the mass continuity step
      if (mask.ice_free_ocean(i, j) && mask.next_to_floating_ice(i, j)) {

        if (mask.floating_ice(i + 1, j)) { N_floating_neighbors += 1; H_average += vH(i + 1, j); }
        if (mask.floating_ice(i - 1, j)) { N_floating_neighbors += 1; H_average += vH(i - 1, j); }
        if (mask.floating_ice(i, j + 1)) { N_floating_neighbors += 1; H_average += vH(i, j + 1); }
        if (mask.floating_ice(i, j - 1)) { N_floating_neighbors += 1; H_average += vH(i, j - 1); }

        if (N_floating_neighbors > 0)
          H_average /= N_floating_neighbors;

        for (int p = -1; p < 2; p += 2) {
          int i_offset = p * offset;
          if (mask.floating_ice(i + i_offset, j) &&
              mask.ice_margin(i + i_offset, j) == false) {
            eigen1 += strain_rates(i + i_offset, j, 0);
            eigen2 += strain_rates(i + i_offset, j, 1);
            M += 1;
          }
        }

        for (int q = -1; q < 2; q += 2) {
          int j_offset = q * offset;
          if (mask.floating_ice(i, j + j_offset) &&
              mask.ice_margin(i , j + j_offset) == false) {
            eigen1 += strain_rates(i, j + j_offset, 0);
            eigen2 += strain_rates(i, j + j_offset, 1);
            M += 1;
          }
        }

        if (M > 0) {
          eigen1 /= M;
          eigen2 /= M;
        }

        PetscScalar
          calving_rate_horizontal = 0.0,
          eigenCalvOffset         = 0.0;
        // eigenCalvOffset allows adjusting the transition from
        // compressive to extensive flow regime

        // Calving law
        //
        // eigen1 * eigen2 has units [s^-2] and calving_rate_horizontal
        // [m*s^1] hence, eigen_calving_K has units [m*s]
        if (eigen2 > eigenCalvOffset &&
            eigen1 > 0.0) { // if spreading in all directions
          calving_rate_horizontal = eigen_calving_K * eigen1 * (eigen2 - eigenCalvOffset);
        }

        // calculate mass loss with respect to the associated ice thickness and the grid size:
        PetscScalar calving_rate = calving_rate_horizontal * H_average / grid.dx; // in m/s

        // apply calving rate at partially filled or empty grid cells
        if (calving_rate > 0.0) {
          vHref(i, j) -= calving_rate * dt; // in m

          if(vHref(i, j) < 0.0) { // i.e. partially filled grid cell has completely calved off

            thk_loss(i, j) = -vHref(i, j); // in m, corresponds to additional ice loss
            vHref(i, j)         = 0.0;

            // additional mass loss will be distributed among
            // N_floating_neighbors:
            if(N_floating_neighbors > 0)
              thk_loss(i, j) /= N_floating_neighbors;
          }
        }

      } // end of "if (ice_free_ocean && next_to_floating)"
    } // j-loop
  } // i-loop
  ierr = thk_loss.end_access(); CHKERRQ(ierr);

  ierr = thk_loss.update_ghosts(); CHKERRQ(ierr);

  ierr = thk_loss.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      PetscScalar thk_loss_ij = 0.0;

      if (mask.floating_ice(i, j) &&
          (thk_loss(i + 1, j) > 0.0 || thk_loss(i - 1, j) > 0.0 ||
           thk_loss(i, j + 1) > 0.0 || thk_loss(i, j - 1) > 0.0)) {

        thk_loss_ij = (thk_loss(i + 1, j) + thk_loss(i - 1, j) +
                       thk_loss(i, j + 1) + thk_loss(i, j - 1));     // in m/s

        // Note PetscMax: we do not account for further calving
        // ice-inwards! Alternatively CFL criterion for time stepping
        // could be adjusted to maximum of calving rate
        vHref(i, j) = PetscMax(vH(i, j) - thk_loss_ij, 0.0); // in m

        vH(i, j) = 0.0;
      }

    } // j-loop
  } // i-loop
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vHref.end_access(); CHKERRQ(ierr);
  ierr = strain_rates.end_access(); CHKERRQ(ierr);
  ierr = thk_loss.end_access(); CHKERRQ(ierr);

  ierr = vH.update_ghosts(vH); CHKERRQ(ierr);

  return 0;
}

//! \brief This calculates the CFL timestep according to the calving rate based
//! on principal strain rates ("eigencalving" with constant K).
PetscErrorCode IceModel::dt_from_eigenCalving() {
  PetscErrorCode ierr;

  double eigenCalvFactor = config.get("eigen_calving_K");

  if (PetscAbs(grid.dx - grid.dy) / PetscMin(grid.dx, grid.dy) > 1e-2) {
    PetscPrintf(grid.com,
                "PISMPIK_ERROR: -eigen_calving using a non-square grid cell does not work (yet);\n"
                "               since it has no direction!!!\n, dx = %f, dy = %f, rel. diff = %f",
                grid.dx, grid.dy, PetscAbs(grid.dx - grid.dy) / PetscMax(grid.dx, grid.dy));
    PISMEnd();
  }

  // About 9 hours which corresponds to 10000 km/year on a 10 km grid
  PetscScalar dt_from_eigencalving_min = 0.001;

  // Distance (grid cells) from calving front where strain rate is evaluated
  PetscInt offset = 2, i0 = 0, j0 = 0;
  PetscScalar
    my_calving_rate_max     = 0.0,
    my_calving_rate_mean    = 0.0,
    my_calving_rate_counter = 0.0;

  MaskQuery mask(vMask);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = strain_rates.begin_access(); CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      // Average of strain-rate eigenvalues in adjacent floating grid cells to
      // be used for eigencalving
      PetscScalar eigen1 = 0.0, eigen2 = 0.0;
      // Number of cells used in computing eigen1 and eigen2:
      PetscInt M = 0;

      // find partially filled or empty grid boxes on the icefree ocean, which
      // have floating ice neighbors after massContExplicitStep
      if ((mask.ice_free_ocean(i, j) &&
           mask.next_to_floating_ice(i, j)) == false)
        continue;

      PetscScalar
        calving_rate_horizontal = 0.0,
        eigenCalvOffset = 0.0;

      for (int p = -1; p < 2; p += 2) {
        int i_offset = p * offset;
        if (mask.floating_ice(i + i_offset, j) &&
            mask.ice_margin(i + i_offset, j) == false) {
          eigen1 += strain_rates(i + i_offset, j, 0);
          eigen2 += strain_rates(i + i_offset, j, 1);
          M += 1;
        }
      }

      for (int q = -1; q < 2; q += 2) {
        int j_offset = q * offset;
        if (mask.floating_ice(i, j + j_offset) &&
            mask.ice_margin(i,   j + j_offset) == false) {
          eigen1 += strain_rates(i, j + j_offset, 0);
          eigen2 += strain_rates(i, j + j_offset, 1);
          M += 1;
        }
      }

      if (M > 0) {
        eigen1 /= M;
        eigen2 /= M;
      }

      // calving law
      if (eigen2 > eigenCalvOffset && eigen1 > 0.0) { // if spreading in all directions
        calving_rate_horizontal = eigenCalvFactor * eigen1 * (eigen2 - eigenCalvOffset);
        my_calving_rate_counter += 1.0;
        my_calving_rate_mean += calving_rate_horizontal;
        if (my_calving_rate_max < calving_rate_horizontal) {
          i0 = i;
          j0 = j;
        }
        my_calving_rate_max = PetscMax(my_calving_rate_max, calving_rate_horizontal);
      } else calving_rate_horizontal = 0.0;

    } // i-loop
  } // j-loop

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = strain_rates.end_access(); CHKERRQ(ierr);

  PetscScalar calving_rate_max = 0.0, calving_rate_mean = 0.0, calving_rate_counter = 0.0;
  ierr = PISMGlobalSum(&my_calving_rate_mean, &calving_rate_mean, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&my_calving_rate_counter, &calving_rate_counter, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&my_calving_rate_max, &calving_rate_max, grid.com); CHKERRQ(ierr);

  calving_rate_mean /= calving_rate_counter;

  PetscScalar denom = calving_rate_max / grid.dx;
  const double epsilon = grid.convert(0.001 / (grid.dx + grid.dy), "seconds", "years");

  dt_from_eigencalving = 1.0 / (denom + epsilon);

  ierr = verbPrintf(2, grid.com,
                    "!!!!! c_rate = %.0f m/year ( dt=%.5f a ) at point %d, %d with mean_c=%.0f m/year over %.0f cells \n",
                    grid.convert(calving_rate_max, "m/s", "m/year"),
                    grid.convert(dt_from_eigencalving, "seconds", "years"),
                    i0, j0,
                    grid.convert(calving_rate_mean, "m/s", "m/year"),
                    calving_rate_counter); CHKERRQ(ierr);

  dt_from_eigencalving = PetscMax(dt_from_eigencalving, dt_from_eigencalving_min);

  return 0;
}

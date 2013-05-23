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
#include "PISMEigenCalving.hh"
#include "PISMStressBalance.hh"
#include "PISMIcebergRemover.hh"

PetscErrorCode IceModel::do_calving() {
  PetscErrorCode ierr;
  bool compute_cumulative_discharge = discharge_flux_2D_cumulative.was_created();
  IceModelVec2S old_H = vWork2d[0],
    old_Href = vWork2d[1];

  if (compute_cumulative_discharge) {
    ierr = vH.copy_to(old_H); CHKERRQ(ierr);
    ierr = vHref.copy_to(old_Href); CHKERRQ(ierr);
  }

  // eigen-calving should go first: it uses the ice velocity field,
  // which is defined at grid points that were icy at the *beginning*
  // of a time-step.
  if (eigen_calving != NULL) {
    ierr = eigen_calving->update(dt, vMask, vHref, vH); CHKERRQ(ierr);
  }

  if (ocean_kill_calving != NULL) {
    ierr = ocean_kill_calving->update(vMask, vH); CHKERRQ(ierr);
  }

  if (thickness_threshold_calving != NULL) {
    ierr = thickness_threshold_calving->update(vMask, vH); CHKERRQ(ierr);
  }

  if (iceberg_remover != NULL) {
    ierr = iceberg_remover->update(vMask, vH); CHKERRQ(ierr);
  }

  // The calls above modify ice thickness and update the mask
  // accordingly.
  ierr = update_surface_elevation(vbed, vH, vh); CHKERRQ(ierr);

  ierr = Href_cleanup(); CHKERRQ(ierr);

  if (compute_cumulative_discharge) {
    ierr = update_cumulative_discharge(vH, old_H, vHref, old_Href); CHKERRQ(ierr);
  }

  return 0;
}

/**
 * Clean up the Href field.
 *
 * Href(i,j) > 0 is allowed only if vH(i,j) == 0 and (i,j) has a
 * floating ice neighbor.
 */
PetscErrorCode IceModel::Href_cleanup() {
  PetscErrorCode ierr;

  if (vHref.was_created() == false)
    return 0;

  MaskQuery mask(vMask);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHref.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {

      if (vH(i, j) > 0 &&
          vHref(i, j) > 0) {
        vH(i, j) += vHref(i, j);
        vHref(i, j) = 0.0;
      }

      if (vHref(i, j) > 0.0 &&
          mask.next_to_floating_ice(i, j) == false) {
        vHref(i, j) = 0.0;
      }

    }
  }

  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vHref.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  return 0;
}

/**
 * Updates the cumulative ice discharge into the ocean.
 *
 * @param thickness current ice thickness
 * @param thickness_old old ice thickness
 * @param Href current "reference ice thickness"
 * @param Href_old old "reference ice thickness"
 *
 * @return 0 on success
 */
PetscErrorCode IceModel::update_cumulative_discharge(IceModelVec2S &thickness,
                                                     IceModelVec2S &thickness_old,
                                                     IceModelVec2S &Href,
                                                     IceModelVec2S &Href_old) {
  PetscErrorCode ierr;

  MaskQuery mask(vMask);

  const double ice_density = config.get("ice_density");

  ierr = thickness.begin_access(); CHKERRQ(ierr);
  ierr = thickness_old.begin_access(); CHKERRQ(ierr);
  ierr = Href.begin_access(); CHKERRQ(ierr);
  ierr = Href_old.begin_access(); CHKERRQ(ierr);
  ierr = discharge_flux_2D_cumulative.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (mask.ice_free_ocean(i,j)) {
        double
          delta_H    = thickness(i,j) - thickness_old(i,j),
          delta_Href = Href(i,j) - Href_old(i,j),
          discharge  = (delta_H + delta_Href) * cell_area(i,j) * ice_density;

        discharge_flux_2D_cumulative(i,j) += discharge;
      }
    }
  }

  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = discharge_flux_2D_cumulative.end_access(); CHKERRQ(ierr);
  ierr = Href_old.end_access(); CHKERRQ(ierr);
  ierr = Href.end_access(); CHKERRQ(ierr);
  ierr = thickness_old.end_access(); CHKERRQ(ierr);
  ierr = thickness.end_access(); CHKERRQ(ierr);

  return 0;
}

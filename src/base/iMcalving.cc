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

  if (config.get_flag("do_eigen_calving") && config.get_flag("use_ssa_velocity")) {
    ierr = eigen_calving->update(dt, vMask, vHref, vH); CHKERRQ(ierr);
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

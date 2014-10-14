// Copyright (C) 2004--2014 Torsten Albrecht and Constantine Khroulev
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
#include <assert.h>

#include "iceModel.hh"
#include "pism_signal.h"
#include "Mask.hh"
#include "PISMOcean.hh"
#include "PISMOceanKill.hh"
#include "PISMFloatKill.hh"
#include "PISMCalvingAtThickness.hh"
#include "PISMEigenCalving.hh"
#include "PISMStressBalance.hh"
#include "PISMIcebergRemover.hh"

namespace pism {

PetscErrorCode IceModel::do_calving() {
  PetscErrorCode ierr;
  bool compute_cumulative_discharge = discharge_flux_2D_cumulative.was_created();
  IceModelVec2S
    &old_H    = vWork2d[0],
    &old_Href = vWork2d[1];

  if (compute_cumulative_discharge) {
    ierr = ice_thickness.copy_to(old_H); CHKERRQ(ierr);
    if (vHref.was_created()) {
      ierr = vHref.copy_to(old_Href); CHKERRQ(ierr);
    }
  }

  // eigen-calving should go first: it uses the ice velocity field,
  // which is defined at grid points that were icy at the *beginning*
  // of a time-step.
  if (eigen_calving != NULL) {
    ierr = eigen_calving->update(dt, vMask, vHref, ice_thickness); CHKERRQ(ierr);
  }

  if (ocean_kill_calving != NULL) {
    ierr = ocean_kill_calving->update(vMask, ice_thickness); CHKERRQ(ierr);
  }

  if (float_kill_calving != NULL) {
    ierr = float_kill_calving->update(vMask, ice_thickness); CHKERRQ(ierr);
  }

  if (thickness_threshold_calving != NULL) {
    ierr = thickness_threshold_calving->update(vMask, ice_thickness); CHKERRQ(ierr);
  }

  // This call removes icebergs, too.
  ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);

  ierr = Href_cleanup(); CHKERRQ(ierr);

  ierr = update_cumulative_discharge(ice_thickness, old_H, vHref, old_Href); CHKERRQ(ierr);

  return 0;
}

/**
 * Clean up the Href field.
 *
 * Href(i,j) > 0 is allowed only if ice_thickness(i,j) == 0 and (i,j) has a
 * floating ice neighbor.
 */
PetscErrorCode IceModel::Href_cleanup() {
  if (vHref.was_created() == false)
    return 0;

  MaskQuery mask(vMask);

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(vHref);
  list.add(vMask);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ice_thickness(i, j) > 0 && vHref(i, j) > 0) {
      ice_thickness(i, j) += vHref(i, j);
      vHref(i, j) = 0.0;
    }

    if (vHref(i, j) > 0.0 && mask.next_to_ice(i, j) == false) {
      vHref(i, j) = 0.0;
    }
  }

  return 0;
}

/**
 * Updates the cumulative ice discharge into the ocean.
 *
 * Units: kg, computed as thickness [m] * cell_area [m2] * density [kg m-3].
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
  const bool
    update_2d_discharge = discharge_flux_2D_cumulative.was_created(),
    use_Href = Href.was_created() && Href_old.was_created();
  double my_total_discharge = 0.0, total_discharge;

  IceModelVec::AccessList list;
  list.add(thickness);
  list.add(thickness_old);
  list.add(vMask);
  list.add(cell_area);

  if (update_2d_discharge) {
    list.add(discharge_flux_2D_cumulative);
  }

  if (use_Href) {
    list.add(Href);
    list.add(Href_old);
  }

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free_ocean(i,j)) {
      double
        delta_H    = thickness(i,j) - thickness_old(i,j),
        delta_Href = 0.0,
        discharge  = 0.0;

      if (use_Href == true) {
        delta_Href = Href(i,j) - Href_old(i,j);
        // Only count mass loss. (A cell may stay "partially-filled"
        // for several time-steps as the calving front advances. In
        // this case delta_Href is real, but does not correspond to
        // either loss or gain of mass.)
        if (delta_Href > 0.0)
          delta_Href = 0.0;
      } else {
        delta_Href = 0.0;
      }

      discharge = (delta_H + delta_Href) * cell_area(i,j) * ice_density;

      if (update_2d_discharge)
        discharge_flux_2D_cumulative(i,j) += discharge;

      my_total_discharge += discharge;
    }
  }

  ierr = GlobalSum(grid.com, &my_total_discharge,  &total_discharge); CHKERRQ(ierr);

  this->discharge_flux_cumulative += total_discharge;

  return 0;
}

} // end of namespace pism

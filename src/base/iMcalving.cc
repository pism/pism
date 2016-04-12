// Copyright (C) 2004--2016 Torsten Albrecht and Constantine Khroulev
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
#include <cassert>

#include "iceModel.hh"
#include "base/calving/PISMCalvingAtThickness.hh"
#include "base/calving/PISMEigenCalving.hh"
#include "base/calving/PISMFloatKill.hh"
#include "base/calving/PISMIcebergRemover.hh"
#include "base/calving/PISMOceanKill.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/pism_const.hh"
#include "coupler/PISMOcean.hh"
#include "base/util/pism_utilities.hh"

namespace pism {

void IceModel::do_calving() {

  IceModelVec2S
    &old_H    = vWork2d[0],
    &old_Href = vWork2d[1];

  {
    old_H.copy_from(m_ice_thickness);
    if (vHref.was_created()) {
      old_Href.copy_from(vHref);
    }
  }

  // eigen-calving should go first: it uses the ice velocity field,
  // which is defined at grid points that were icy at the *beginning*
  // of a time-step.
  if (eigen_calving != NULL) {
    eigen_calving->update(m_dt, m_cell_type, vHref, m_ice_thickness);
  }

  if (ocean_kill_calving != NULL) {
    ocean_kill_calving->update(m_cell_type, m_ice_thickness);
  }

  if (float_kill_calving != NULL) {
    float_kill_calving->update(m_cell_type, m_ice_thickness);
  }

  if (thickness_threshold_calving != NULL) {
    thickness_threshold_calving->update(m_cell_type, m_ice_thickness);
  }

  // This call removes icebergs, too.
  updateSurfaceElevationAndMask();

  Href_cleanup();

  update_cumulative_discharge(m_ice_thickness, old_H, vHref, old_Href);
}

/**
 * Clean up the Href field.
 *
 * Href(i,j) > 0 is allowed only if ice_thickness(i,j) == 0 and (i,j) has a
 * floating ice neighbor.
 */
void IceModel::Href_cleanup() {
  if (not vHref.was_created()) {
    return;
  }

  IceModelVec::AccessList list;
  list.add(m_ice_thickness);
  list.add(vHref);
  list.add(m_cell_type);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_ice_thickness(i, j) > 0 && vHref(i, j) > 0) {
      m_ice_thickness(i, j) += vHref(i, j);
      vHref(i, j) = 0.0;
    }

    if (vHref(i, j) > 0.0 && not m_cell_type.next_to_ice(i, j)) {
      vHref(i, j) = 0.0;
    }
  }
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
void IceModel::update_cumulative_discharge(const IceModelVec2S &thickness,
                                           const IceModelVec2S &thickness_old,
                                           const IceModelVec2S &Href,
                                           const IceModelVec2S &Href_old) {

  const double ice_density = m_config->get_double("ice_density");
  const bool
    use_Href = Href.was_created() && Href_old.was_created();
  double my_total_discharge = 0.0, total_discharge;

  IceModelVec::AccessList list;
  list.add(thickness);
  list.add(thickness_old);
  list.add(m_cell_type);
  list.add(m_cell_area);


  list.add(m_discharge_flux_2D_cumulative);

  if (use_Href) {
    list.add(Href);
    list.add(Href_old);
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_cell_type.ice_free_ocean(i,j)) {
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
        if (delta_Href > 0.0) {
          delta_Href = 0.0;
        }
      } else {
        delta_Href = 0.0;
      }

      discharge = (delta_H + delta_Href) * m_cell_area(i,j) * ice_density;

      m_discharge_flux_2D_cumulative(i,j) += discharge;

      my_total_discharge += discharge;
    }
  }

  total_discharge = GlobalSum(m_grid->com, my_total_discharge);

  this->discharge_flux_cumulative += total_discharge;
}

} // end of namespace pism

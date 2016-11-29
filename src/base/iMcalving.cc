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
#include "base/calving/CalvingAtThickness.hh"
#include "base/calving/EigenCalving.hh"
#include "base/calving/vonMisesCalving.hh"
#include "base/calving/FrontalMelt.hh"
#include "base/calving/FloatKill.hh"
#include "base/calving/IcebergRemover.hh"
#include "base/calving/OceanKill.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/pism_const.hh"
#include "coupler/PISMOcean.hh"
#include "base/util/pism_utilities.hh"
#include "earth/PISMBedDef.hh"


namespace pism {

void IceModel::do_calving() {

  IceModelVec2S
    &old_H    = m_work2d[0],
    &old_Href = m_work2d[1];

  {
    old_H.copy_from(m_ice_thickness);
    if (m_Href.was_created()) {
      old_Href.copy_from(m_Href);
    } else {
      old_Href.set(0.0);
    }
  }

  // eigen-calving should go first: it uses the ice velocity field,
  // which is defined at grid points that were icy at the *beginning*
  // of a time-step.
  if (m_eigen_calving != NULL) {
    m_eigen_calving->update(m_dt,
                            m_ocean->sea_level_elevation(),
                            m_beddef->bed_elevation(),
                            m_cell_type,
                            m_Href,
                            m_ice_thickness);
  }

  if (m_vonmises_calving != NULL) {
    m_vonmises_calving->update(m_dt,
                               m_ocean->sea_level_elevation(),
                               m_beddef->bed_elevation(),
                               m_cell_type,
                               m_Href,
                               m_ice_thickness);
  }

  if (m_frontal_melt != NULL) {
    m_frontal_melt->update(m_dt,
                           m_ocean->sea_level_elevation(),
                           m_beddef->bed_elevation(),
                           m_cell_type,
                           m_Href,
                           m_ice_thickness);
  }

  if (m_ocean_kill_calving != NULL) {
    m_ocean_kill_calving->update(m_cell_type, m_ice_thickness);
  }

  if (m_float_kill_calving != NULL) {
    m_float_kill_calving->update(m_cell_type, m_ice_thickness);
  }

  if (m_thickness_threshold_calving != NULL) {
    m_thickness_threshold_calving->update(m_cell_type, m_ice_thickness);
  }

  // This call removes icebergs, too.
  enforce_consistency_of_geometry();

  Href_cleanup();

  // note that Href_cleanup() changes ice thickness, so we have to
  // update the mask and surface elevation.
  enforce_consistency_of_geometry();

  update_cumulative_discharge(m_ice_thickness, old_H, m_Href, old_Href);
}

/**
 * Clean up the Href field.
 *
 * Href(i,j) > 0 is allowed only if ice_thickness(i,j) == 0 and (i,j) has a
 * floating ice neighbor.
 */
void IceModel::Href_cleanup() {
  if (not m_Href.was_created()) {
    return;
  }

  IceModelVec::AccessList list;
  list.add(m_ice_thickness);
  list.add(m_Href);
  list.add(m_cell_type);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_ice_thickness(i, j) > 0 && m_Href(i, j) > 0) {
      m_ice_thickness(i, j) += m_Href(i, j);
      m_Href(i, j) = 0.0;
    }

    if (m_Href(i, j) > 0.0 && not m_cell_type.next_to_ice(i, j)) {
      m_Href(i, j) = 0.0;
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

  const double ice_density = m_config->get_double("constants.ice.density");
  const bool use_Href = Href.was_created() and Href_old.was_created();
  double total_discharge = 0.0;

  IceModelVec::AccessList list;
  list.add(thickness);
  list.add(thickness_old);
  list.add(m_cell_type);
  list.add(m_cell_area);


  list.add(m_cumulative_flux_fields.discharge);

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

      if (use_Href) {
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

      m_cumulative_flux_fields.discharge(i, j) += discharge;

      total_discharge += discharge;
    }
  }

  m_cumulative_fluxes.discharge += GlobalSum(m_grid->com, total_discharge);
}

} // end of namespace pism

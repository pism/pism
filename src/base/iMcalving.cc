// Copyright (C) 2004--2017 Torsten Albrecht and Constantine Khroulev
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

  // if no calving method was selected, stop
  if (m_config->get_string("calving.methods").empty()) {
    return;
  }

  if (not m_config->get_boolean("geometry.part_grid.enabled")) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "calving requires geometry.part_grid.enabled");
  }

  IceModelVec2S
    &old_H    = m_work2d[0],
    &old_Href = m_work2d[1];

  {
    old_H.copy_from(m_geometry.ice_thickness);
    old_Href.copy_from(m_geometry.ice_area_specific_volume);
  }

  // eigen-calving should go first: it uses the ice velocity field,
  // which is defined at grid points that were icy at the *beginning*
  // of a time-step.
  if (m_eigen_calving != NULL) {
    m_eigen_calving->update(m_dt,
                            m_ocean->sea_level_elevation(),
                            m_ssa_dirichlet_bc_mask,
                            m_geometry.bed_elevation,
                            m_geometry.cell_type,
                            m_geometry.ice_area_specific_volume,
                            m_geometry.ice_thickness);
  }

  if (m_vonmises_calving != NULL) {
    m_vonmises_calving->update(m_dt,
                               m_ocean->sea_level_elevation(),
                               m_ssa_dirichlet_bc_mask,
                               m_geometry.bed_elevation,
                               m_geometry.cell_type,
                               m_geometry.ice_area_specific_volume,
                               m_geometry.ice_thickness);
  }

  if (m_frontal_melt != NULL) {
    m_frontal_melt->update(m_dt,
                           m_ocean->sea_level_elevation(),
                           m_ssa_dirichlet_bc_mask,
                           m_geometry.bed_elevation,
                           m_geometry.cell_type,
                           m_geometry.ice_area_specific_volume,
                           m_geometry.ice_thickness);
  }

  if (m_ocean_kill_calving != NULL) {
    m_ocean_kill_calving->update(m_geometry.cell_type, m_geometry.ice_thickness);
  }

  if (m_float_kill_calving != NULL) {
    m_float_kill_calving->update(m_geometry.cell_type, m_geometry.ice_thickness);
  }

  if (m_thickness_threshold_calving != NULL) {
    m_thickness_threshold_calving->update(m_geometry.cell_type, m_geometry.ice_thickness);
  }

  // This call removes icebergs, too.
  enforce_consistency_of_geometry();

  Href_cleanup();

  // note that Href_cleanup() changes ice thickness, so we have to
  // update the mask and surface elevation.
  enforce_consistency_of_geometry();

  compute_discharge(m_geometry.ice_thickness,
                    m_geometry.ice_area_specific_volume,
                    old_H, old_Href,
                    m_dischange);
}

/**
 * Clean up the Href field.
 *
 * Href(i,j) > 0 is allowed only if ice_thickness(i,j) == 0 and (i,j) has a
 * floating ice neighbor.
 */
void IceModel::Href_cleanup() {

  IceModelVec2S
    &V = m_geometry.ice_area_specific_volume,
    &H = m_geometry.ice_thickness;

  IceModelVec::AccessList list{&H, &V, &m_geometry.cell_type};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (H(i, j) > 0 and V(i, j) > 0) {
      H(i, j) += V(i, j);
      V(i, j) = 0.0;
    }

    if (V(i, j) > 0.0 and not m_geometry.cell_type.next_to_ice(i, j)) {
      V(i, j) = 0.0;
    }
  }
}

/**
 * Compute the ice discharge into the ocean during the current time step.
 *
 * Units: kg, computed as thickness [m] * cell_area [m2] * density [kg m-3].
 *
 * @param thickness current ice thickness
 * @param Href current "reference ice thickness"
 * @param thickness_old old ice thickness
 * @param Href_old old "reference ice thickness"
 * @param[out] output computed discharge during the current time step
 */
void IceModel::compute_discharge(const IceModelVec2S &thickness,
                                 const IceModelVec2S &Href,
                                 const IceModelVec2S &thickness_old,
                                 const IceModelVec2S &Href_old,
                                 IceModelVec2S &output) {

  const double ice_density = m_config->get_double("constants.ice.density");

  IceModelVec::AccessList list{&thickness, &thickness_old,
      &Href, &Href_old, &m_geometry.cell_area, &output};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      H_old     = thickness_old(i, j) + Href_old(i, j),
      H_new     = thickness(i, j) + Href(i, j),
      discharge = (H_new - H_old) * m_geometry.cell_area(i, j) * ice_density;

    output(i, j) += discharge;
  }
}

} // end of namespace pism

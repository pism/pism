/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "CalvingAtThickness.hh"

#include "pism/util/Mask.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"


namespace pism {

//! @brief Calving and iceberg removal code.
namespace calving {

CalvingAtThickness::CalvingAtThickness(IceGrid::ConstPtr g)
  : Component(g) {

  m_old_mask.create(m_grid, "old_mask", WITH_GHOSTS, 1);

  m_calving_threshold.create(m_grid, "calving_threshold", WITHOUT_GHOSTS);

  m_calving_threshold.set_attrs("diagnostic",
                                "threshold used by the 'calving at threshold' calving method",
                                "m",
                                ""); // no standard name
  m_calving_threshold.set_time_independent(true);
}

CalvingAtThickness::~CalvingAtThickness() {
  // empty
}


void CalvingAtThickness::init() {

  m_log->message(2, "* Initializing the 'calving at a threshold thickness' mechanism...\n");

  std::string threshold_file = m_config->get_string("calving.thickness_calving.threshold_file");

  double calving_threshold = m_config->get_double("calving.thickness_calving.threshold");

  m_calving_threshold.set(calving_threshold);

  if (not threshold_file.empty()) {
    m_log->message(2,
                   "  Reading thickness calving threshold from file '%s'...\n",
                   threshold_file.c_str());

    m_calving_threshold.regrid(threshold_file, CRITICAL);
  } else {
    m_log->message(2,
                   "  Thickness threshold: %3.3f meters.\n", calving_threshold);
  }
}

/**
 * Updates ice cover mask and the ice thickness according to the
 * calving rule removing ice at the shelf front that is thinner than a
 * given threshold.
 *
 * @param[in,out] pism_mask ice cover mask
 * @param[in,out] ice_thickness ice thickness
 *
 * @return 0 on success
 */
void CalvingAtThickness::update(IceModelVec2CellType &pism_mask,
                                IceModelVec2S &ice_thickness) {

  // this call fills ghosts of m_old_mask
  m_old_mask.copy_from(pism_mask);

  IceModelVec::AccessList list{&pism_mask, &ice_thickness, &m_old_mask, &m_calving_threshold};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_old_mask.floating_ice(i, j)           &&
        m_old_mask.next_to_ice_free_ocean(i, j) &&
        ice_thickness(i, j) < m_calving_threshold(i, j)) {
      ice_thickness(i, j) = 0.0;
      pism_mask(i, j)     = MASK_ICE_FREE_OCEAN;
    }
  }

  pism_mask.update_ghosts();
  ice_thickness.update_ghosts();
}

const IceModelVec2S& CalvingAtThickness::threshold() const {
  return m_calving_threshold;
}

DiagnosticList CalvingAtThickness::diagnostics_impl() const {
  return {{"calving_threshold", Diagnostic::wrap(m_calving_threshold)}};
}

} // end of namespace calving
} // end of namespace pism

/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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
#include "base/util/pism_options.hh"
#include "base/util/Mask.hh"
#include "base/util/error_handling.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_const.hh"
#include "base/part_grid_threshold_thickness.hh"


namespace pism {

//! @brief Calving and iceberg removal code.
namespace calving {

CalvingAtThickness::CalvingAtThickness(IceGrid::ConstPtr g)
  : Component(g) {
  m_calving_threshold_scalar = m_config->get_double("calving.thickness_calving.threshold");

  m_old_mask.create(m_grid, "old_mask", WITH_GHOSTS, 1);

  m_calving_threshold.create(m_grid, "calving_threshold", WITH_GHOSTS,
                           m_config->get_double("grid.max_stencil_width"));

  m_calving_threshold.set_attrs("internal",
                              "calving thickness threshold",
                              "m", "");
  m_calving_threshold.set_time_independent(true);


}

CalvingAtThickness::~CalvingAtThickness() {
  // empty
}


void CalvingAtThickness::init() {
  
  m_log->message(2,
             "* Initializing the 'calving at a threshold thickness' mechanism...\n"
             "  thickness threshold: %3.3f meters\n", m_calving_threshold_scalar);

  options::String thickness_calving_threshold_file("-thickness_calving_threshold_file",
                                  "Specifies a file to get thickness calving threshold from");

  if (thickness_calving_threshold_file.is_set()) {

      m_log->message(2,
                     "* Option '-thickness_calving_threshold_file' found...\n"
                     "  reading thickness calving threshold from file %s\n"
                     "  this will override thickness threshold: %3.3f meters\n", thickness_calving_threshold_file->c_str(), m_calving_threshold_scalar);

      m_calving_threshold.regrid(thickness_calving_threshold_file, CRITICAL);

  } else {

    // set m_calving_threshold to scalar m_calving_threshold_scalar everywhere
    m_calving_threshold.set(m_calving_threshold_scalar);
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

std::map<std::string, Diagnostic::Ptr> CalvingAtThickness::diagnostics_impl() const {
  return {{"calving_threshold", Diagnostic::Ptr(new CalvingAtThickness_threshold(this))}};
}

CalvingAtThickness_threshold::CalvingAtThickness_threshold(const CalvingAtThickness *m)
  : Diag<CalvingAtThickness>(m) {

  /* set metadata: */
  m_vars.push_back(SpatialVariableMetadata(m_sys, "calving_threshold"));

  set_attrs("threshold used by the 'calving at threshold' calving method ", "",
            "", "", 0);
}

IceModelVec::Ptr CalvingAtThickness_threshold::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "calving_threshold", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];
  result->metadata(0).set_output_type(PISM_INT);

  result->copy_from(model->threshold());

  return result;
}

  
} // end of namespace calving
} // end of namespace pism

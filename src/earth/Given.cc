/* Copyright (C) 2020 PISM Authors
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

#include "Given.hh"

namespace pism {
namespace bed {

Given::Given(IceGrid::ConstPtr grid)
  : BedDef(grid),
    m_topg_reference(grid, "topg", WITHOUT_GHOSTS) {

  m_topg_reference.set_attrs("bed_deformation", "reference bed elevation",
                             "meters",
                             "meters",
                             "bedrock_altitude", 0);

  auto filename = m_config->get_string("bed_deformation.given.file");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    // point-wise time-series are not used
    unsigned int evaluations_per_year = 1;

    // periodic inputs are not supported
    bool periodic = false;

    File file(m_grid->com, filename, PISM_NETCDF3, PISM_READONLY);

    m_topg_delta = IceModelVec2T::ForcingField(m_grid,
                                               file,
                                               "topg_delta",
                                               "", // no standard name
                                               buffer_size,
                                               evaluations_per_year,
                                               periodic,
                                               LINEAR);
    m_topg_delta->set_attrs("bed_deformation",
                            "two-dimensional bed elevation changes",
                            "meters", "meters", "", 0);
  }
}

Given::~Given() {
  // empty
}

void Given::init_impl(const InputOptions &opts, const IceModelVec2S &ice_thickness,
                      const IceModelVec2S &sea_level_elevation) {
  (void) ice_thickness;
  (void) sea_level_elevation;

  m_log->message(2,
                 "* Initializing the bed model using a given topography change history...\n");

  BedDef::init_impl(opts, ice_thickness, sea_level_elevation);

  {
    auto reference_filename = m_config->get_string("bed_deformation.given.reference_file");
    m_topg_reference.regrid(reference_filename, CRITICAL); // fails if not found!
  }

  {
    auto filename = m_config->get_string("bed_deformation.given.file");
    m_topg_delta->init(filename, 0.0, 0.0);
  }
}

void Given::update_impl(const IceModelVec2S &ice_thickness,
                        const IceModelVec2S &sea_level_elevation,
                        double t, double dt) {
  (void) ice_thickness;
  (void) sea_level_elevation;

  m_topg_delta->update(t, dt);
  m_topg_delta->average(t, dt);

  m_topg_last.copy_from(m_topg);

  m_topg.copy_from(m_topg_reference);
  m_topg.add(1.0, *m_topg_delta);

  m_uplift.copy_from(m_topg_last);
  m_uplift.add(-1.0, m_topg);
}

} // end of namespace bed
} // end of namespace pism

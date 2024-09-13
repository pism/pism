/* Copyright (C) 2020, 2021, 2022, 2023, 2024 PISM Authors
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

#include "pism/earth/Given.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {
namespace bed {

Given::Given(std::shared_ptr<const Grid> grid)
  : BedDef(grid, "'prescribed topography change history'"),
    m_topg_reference(grid, "topg") {

  m_topg_reference.metadata(0)
      .long_name("reference bed elevation")
      .units("meters")
      .standard_name("bedrock_altitude");

  auto filename = m_config->get_string("bed_deformation.given.file");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    // periodic inputs are not supported
    bool periodic = false;

    File file(m_grid->com, filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_topg_delta = std::make_shared<array::Forcing>(m_grid, file, "topg_delta",
                                                    "", // no standard name
                                                    buffer_size, periodic, LINEAR);
    m_topg_delta->metadata(0)
        .long_name("two-dimensional bed elevation changes")
        .units("meters");
  }
}

void Given::init_impl(const InputOptions & /*opts*/, const array::Scalar & /*ice_thickness*/,
                      const array::Scalar & /*sea_level_elevation*/) {
  {
    auto reference_filename = m_config->get_string("bed_deformation.given.reference_file");
    m_topg_reference.regrid(reference_filename, io::Default::Nil()); // fails if not found!
  }

  {
    auto filename   = m_config->get_string("bed_deformation.given.file");
    bool periodic_p = false;
    m_topg_delta->init(filename, periodic_p);
  }
}

void Given::bootstrap_impl(const array::Scalar & /*bed_elevation*/,
                           const array::Scalar & /*bed_uplift*/,
                           const array::Scalar & /*ice_thickness*/,
                           const array::Scalar & /*sea_level_elevation*/) {
  // empty
}

void Given::update_impl(const array::Scalar & /*load*/, double t, double dt) {
  m_topg_delta->update(t, dt);
  m_topg_delta->average(t, dt);

  m_topg_reference.add(1.0, *m_topg_delta, m_topg);

  // mark m_topg as "modified"
  m_topg.inc_state_counter();
}

} // end of namespace bed
} // end of namespace pism

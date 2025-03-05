// Copyright (C) 2010, 2011, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2022, 2023, 2024 Constantine Khroulev
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

#include "pism/earth/BedDef.hh"
#include "pism/util/Grid.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace bed {

PointwiseIsostasy::PointwiseIsostasy(std::shared_ptr<const Grid> grid)
    : BedDef(grid, "pointwise isostasy"), m_load_last(m_grid, "load_last") {
  // empty
}

void PointwiseIsostasy::init_impl(const InputOptions &/*opts*/, const array::Scalar &ice_thickness,
                                  const array::Scalar &sea_level_elevation) {
  // store the initial load
  m_load_last.set(0.0);
  accumulate_load(m_topg, ice_thickness, sea_level_elevation, 1.0, m_load_last);
}


void PointwiseIsostasy::bootstrap_impl(const array::Scalar &bed_elevation,
                                       const array::Scalar &/*bed_uplift*/,
                                       const array::Scalar &ice_thickness,
                                       const array::Scalar &sea_level_elevation) {
  // store initial load and bed elevation
  m_load_last.set(0.0);
  accumulate_load(bed_elevation, ice_thickness, sea_level_elevation, 1.0, m_load_last);
  m_topg_last.copy_from(bed_elevation);
}

//! Updates the pointwise isostasy model.
/*!
 * Inputs:
 *
 * - ice thickness
 * - sea level
 * - old bed elevation
 * - old load
 *
 * Outputs:
 *
 * - new bed elevation
 * - updated "old" load
 */
void PointwiseIsostasy::update_impl(const array::Scalar &load,
                                    double /*t*/, double /*dt*/) {
  const double
    mantle_density = m_config->get_number("bed_deformation.mantle_density"),
    load_density   = m_config->get_number("constants.ice.density"),
    f              = load_density / mantle_density;

  // Our goal: topg_{n+1} = topg_{n} - f*(load(topg_{n}) - load_{n-1})

  array::AccessScope list{ &m_topg, &m_topg_last, &load, &m_load_last };

  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      m_topg(i, j) = m_topg_last(i, j) - f * (load(i, j) - m_load_last(i, j));
      m_load_last(i, j) = load(i, j);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  // mark m_topg as "modified"
  m_topg.inc_state_counter();
}

} // end of namespace bed
} // end of namespace pism

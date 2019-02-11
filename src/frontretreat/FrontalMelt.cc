/* Copyright (C) 2016, 2017, 2018, 2019 PISM Authors
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

#include "FrontalMelt.hh"

#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {

FrontalMelt::FrontalMelt(IceGrid::ConstPtr grid)
  : FrontRetreat(grid, 1) {
  // empty
}

FrontalMelt::~FrontalMelt() {
  // empty
}

void FrontalMelt::init() {
  m_log->message(2,
                 "* Initializing the frontal melt parameterization\n"
                 "  using sub-shelf mass flux from an ocean model...\n");
}

const IceModelVec2S &FrontalMelt::retreat_rate() const {
  return m_retreat_rate;
}

DiagnosticList FrontalMelt::diagnostics_impl() const {
  return {{"frontal_melt_rate", Diagnostic::wrap(m_retreat_rate)}};
}

/*!
 * Convert provided melt rate into the corresponding rate of retreat, considering which
 * part of the front is submerged.
 */
void FrontalMelt::update(const Geometry &geometry, const IceModelVec2S &frontal_melt_rate) {

  GeometryCalculator gc(*m_config);

  const IceModelVec2S
    &bed_elevation       = geometry.bed_elevation,
    &surface_elevation   = geometry.ice_surface_elevation,
    &ice_thickness       = geometry.ice_thickness,
    &sea_level_elevation = geometry.sea_level_elevation;
  const IceModelVec2CellType &cell_type = geometry.cell_type;

  const double
    ice_density = m_config->get_double("constants.ice.density"),
    alpha       = ice_density / m_config->get_double("constants.sea_water.density");

  IceModelVec::AccessList list{&cell_type, &frontal_melt_rate, &sea_level_elevation,
      &bed_elevation, &surface_elevation, &ice_thickness, &m_retreat_rate};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.ice_free_ocean(i, j) and cell_type.next_to_ice(i, j)) {
        const double
          bed       = bed_elevation(i, j),
          sea_level = sea_level_elevation(i, j);

        auto H = ice_thickness.star(i, j);
        auto h = surface_elevation.star(i, j);
        auto M = cell_type.int_star(i, j);

        double H_threshold = part_grid_threshold_thickness(M, H, h, bed);

        int m = gc.mask(sea_level, bed, H_threshold);

        double H_submerged = (mask::grounded(m) ?
                              std::max(sea_level - bed, 0.0) :
                              alpha * H_threshold);

        m_retreat_rate(i, j) = (H_submerged / H_threshold) * frontal_melt_rate(i, j);
      } else {
        m_retreat_rate(i, j) = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

} // end of namespace pism

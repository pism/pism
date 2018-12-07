/* Copyright (C) 2016, 2017, 2018 PISM Authors
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
  : CalvingFrontRetreat(grid, 1) {
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

DiagnosticList FrontalMelt::diagnostics_impl() const {
  return {{"frontal_melt_rate",
        Diagnostic::Ptr(new CalvingRate(this, "frontal_melt_rate",
                                        "horizontal front retreat rate due to melt"))}};
}

void FrontalMelt::compute_calving_rate(const CalvingInputs &inputs,
                                       IceModelVec2S &result) const {

  prepare_mask(inputs.geometry->cell_type, m_mask);

  GeometryCalculator gc(*m_config);

  const IceModelVec2S &shelf_base_mass_flux = *inputs.shelf_base_mass_flux;

  const IceModelVec2S
    &bed_elevation       = inputs.geometry->bed_elevation,
    &surface_elevation   = inputs.geometry->ice_surface_elevation,
    &ice_thickness       = inputs.geometry->ice_thickness,
    &sea_level_elevation = inputs.geometry->sea_level_elevation;

  const double
    ice_density = m_config->get_double("constants.ice.density"),
    alpha       = ice_density / m_config->get_double("constants.sea_water.density");

  IceModelVec::AccessList list{&m_mask, &shelf_base_mass_flux, &sea_level_elevation,
      &bed_elevation, &surface_elevation, &ice_thickness, &result};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_mask.ice_free_ocean(i, j) and m_mask.next_to_ice(i, j)) {
        const double
          bed       = bed_elevation(i, j),
          sea_level = sea_level_elevation(i, j);

        auto H = ice_thickness.star(i, j);
        auto h = surface_elevation.star(i, j);
        auto M = m_mask.int_star(i, j);

        double H_threshold = part_grid_threshold_thickness(M, H, h, bed);

        int m = gc.mask(sea_level, bed, H_threshold);

        double H_submerged = (mask::grounded(m) ?
                              std::max(sea_level - bed, 0.0) :
                              alpha * H_threshold);

        result(i, j) = (H_submerged / H_threshold) * shelf_base_mass_flux(i, j) / ice_density;
      } else {
        result(i, j) = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

} // end of namespace pism

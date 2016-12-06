/* Copyright (C) 2016 PISM Authors
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
#include "base/util/PISMVars.hh"
#include "base/part_grid_threshold_thickness.hh"

namespace pism {

FrontalMelt::FrontalMelt(IceGrid::ConstPtr g, const ocean::OceanModel *ocean_model)
  : CalvingFrontRetreat(g, 1), m_ocean(ocean_model) {
  m_shelf_base_mass_flux.create(m_grid, "shelf_base_mass_flux", WITHOUT_GHOSTS);
  m_shelf_base_mass_flux.set_attrs("internal", "sub-shelf mass flux", "kg m-2 s-1", "");
}

FrontalMelt::~FrontalMelt() {
  // empty
}

void FrontalMelt::init() {
  m_log->message(2,
                 "* Initializing the frontal melt parameterization\n"
                 "  using sub-shelf mass flux from an ocean model...\n");
}

std::map<std::string, Diagnostic::Ptr> FrontalMelt::diagnostics_impl() const {
  return {{"frontal_melt_rate",
        Diagnostic::Ptr(new CalvingRate(this, "frontal_melt_rate",
                                        "horizontal front retreat rate due to melt"))}};
}

void FrontalMelt::compute_calving_rate(const IceModelVec2CellType &mask,
                                       IceModelVec2S &result) const {
  (void) mask;
  (void) result;

  GeometryCalculator gc(*m_config);

  IceModelVec2S &shelf_base_mass_flux = m_tmp;
  m_ocean->shelf_base_mass_flux(shelf_base_mass_flux);

  const IceModelVec2S
    &bed_elevation     = *m_grid->variables().get_2d_scalar("bedrock_altitude"),
    &surface_elevation = *m_grid->variables().get_2d_scalar("surface_altitude"),
    &ice_thickness     = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  const double
    sea_level   = m_ocean->sea_level_elevation(),
    dx          = m_grid->dx(),
    ice_density = m_config->get_double("constants.ice.density"),
    alpha       = ice_density / m_config->get_double("constants.sea_water.density");

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(shelf_base_mass_flux);
  list.add(bed_elevation);
  list.add(surface_elevation);
  list.add(ice_thickness);
  list.add(result);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask.ice_free_ocean(i, j) and mask.next_to_ice(i, j)) {
        const double bed = bed_elevation(i, j);

        StarStencil<double> H = ice_thickness.star(i, j);
        StarStencil<double> h = surface_elevation.star(i, j);
        StarStencil<int> M = mask.int_star(i, j);

        const double H_threshold = part_grid_threshold_thickness(M, H, h, bed, dx, false);

        const int m = gc.mask(sea_level, bed, H_threshold);

        const double H_submerged = (mask::grounded(m) ? std::max(sea_level - bed, 0.0) :
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

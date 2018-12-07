/* Copyright (C) 2017, 2018 PISM Authors
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

#include "Geometry.hh"

#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/Mask.hh"
#include "pism/geometry/grounded_cell_fraction.hh"

namespace pism {

Geometry::Geometry(IceGrid::ConstPtr grid) {
  // FIXME: ideally these fields should be "global", i.e. without ghosts.
  // (However this may increase communication costs...)
  const unsigned int WIDE_STENCIL = grid->ctx()->config()->get_double("grid.max_stencil_width");

  latitude.create(grid, "lat", WITHOUT_GHOSTS);
  latitude.set_attrs("mapping", "latitude", "degree_north", "latitude");
  latitude.set_time_independent(true);
  latitude.metadata().set_string("coordinates", "");
  latitude.metadata().set_string("grid_mapping", "");
  latitude.metadata().set_doubles("valid_range", {-90.0, 90.0});

  longitude.create(grid, "lon", WITHOUT_GHOSTS);
  longitude.set_attrs("mapping", "longitude", "degree_east", "longitude");
  longitude.set_time_independent(true);
  longitude.metadata().set_string("coordinates", "");
  longitude.metadata().set_string("grid_mapping", "");
  longitude.metadata().set_doubles("valid_range", {-180.0, 180.0});

  bed_elevation.create(grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
  bed_elevation.set_attrs("model_state", "bedrock surface elevation",
                          "m", "bedrock_altitude");

  sea_level_elevation.create(grid, "sea_level", WITH_GHOSTS);
  sea_level_elevation.set_attrs("model_state",
                                "sea level elevation above reference ellipsoid", "meters",
                                "sea_surface_height_above_reference_ellipsoid");

  ice_thickness.create(grid, "thk", WITH_GHOSTS, WIDE_STENCIL);
  ice_thickness.set_attrs("model_state", "land ice thickness",
                          "m", "land_ice_thickness");
  ice_thickness.metadata().set_double("valid_min", 0.0);

  ice_area_specific_volume.create(grid, "ice_area_specific_volume", WITH_GHOSTS);
  ice_area_specific_volume.set_attrs("model_state",
                                     "ice-volume-per-area in partially-filled grid cells",
                                     "m3/m2", "");
  ice_area_specific_volume.metadata().set_string("comment",
                                                 "this variable represents the amount of ice "
                                                 "in a partially-filled cell and not "
                                                 "the corresponding geometry, so thinking "
                                                 "about it as 'thickness' is not helpful");

  cell_type.create(grid, "mask", WITH_GHOSTS, WIDE_STENCIL);
  cell_type.set_attrs("diagnostic", "ice-type (ice-free/grounded/floating/ocean) integer mask",
                      "", "");
  std::vector<double> mask_values = {
    MASK_ICE_FREE_BEDROCK,
    MASK_GROUNDED,
    MASK_FLOATING,
    MASK_ICE_FREE_OCEAN};

  cell_type.metadata().set_doubles("flag_values", mask_values);
  cell_type.metadata().set_string("flag_meanings",
                                  "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
  cell_type.metadata().set_output_type(PISM_BYTE);

  cell_grounded_fraction.create(grid, "cell_grounded_fraction", WITHOUT_GHOSTS);
  cell_grounded_fraction.set_attrs("internal",
                                   "fractional grounded/floating mask (floating=0, grounded=1)",
                                   "", "");

  ice_surface_elevation.create(grid, "usurf", WITH_GHOSTS, WIDE_STENCIL);
  ice_surface_elevation.set_attrs("diagnostic", "ice upper surface elevation",
                                  "m", "surface_altitude");
}

void check_minimum_ice_thickness(const IceModelVec2S &ice_thickness) {
  IceGrid::ConstPtr grid = ice_thickness.grid();

  IceModelVec::AccessList list(ice_thickness);

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (ice_thickness(i, j) < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "H = %e (negative) at point i=%d, j=%d",
                                      ice_thickness(i, j), i, j);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

void Geometry::ensure_consistency(double ice_free_thickness_threshold) {
  IceGrid::ConstPtr grid = ice_thickness.grid();
  Config::ConstPtr config = grid->ctx()->config();

  check_minimum_ice_thickness(ice_thickness);

  IceModelVec::AccessList list{&sea_level_elevation, &bed_elevation,
      &ice_thickness, &ice_area_specific_volume,
      &cell_type, &ice_surface_elevation};

  // first ensure that ice_area_specific_volume is 0 if ice_thickness > 0.
  {
    ParallelSection loop(grid->com);
    try {
      for (Points p(*grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        if (ice_thickness(i, j) > 0.0 and ice_area_specific_volume(i, j) > 0.0) {
          ice_thickness(i, j) += ice_area_specific_volume(i, j);
          ice_area_specific_volume(i, j) = 0.0;
        }
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();
  }

  // compute cell type and surface elevation
  {
    GeometryCalculator gc(*config);
    gc.set_icefree_thickness(ice_free_thickness_threshold);

    ParallelSection loop(grid->com);
    try {
      for (Points p(*grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        int mask = 0;
        gc.compute(sea_level_elevation(i, j), bed_elevation(i, j), ice_thickness(i, j),
                   &mask, &ice_surface_elevation(i, j));
        cell_type(i, j) = mask;
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();
  }

  ice_thickness.update_ghosts();
  ice_area_specific_volume.update_ghosts();
  cell_type.update_ghosts();
  ice_surface_elevation.update_ghosts();

  const double
    ice_density = config->get_double("constants.ice.density"),
    ocean_density = config->get_double("constants.sea_water.density");

  compute_grounded_cell_fraction(ice_density,
                                 ocean_density,
                                 sea_level_elevation,
                                 ice_thickness,
                                 bed_elevation,
                                 cell_grounded_fraction);
}

} // end of namespace pism

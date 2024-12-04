/* Copyright (C) 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024 PISM Authors
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

#include <functional>

#include "pism/geometry/Geometry.hh"

#include "pism/util/array/CellType.hh"
#include "pism/util/Mask.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/geometry/grounded_cell_fraction.hh"
#include "pism/util/Context.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {

Geometry::Geometry(const std::shared_ptr<const Grid> &grid)
  // FIXME: ideally these fields should be "global", i.e. without ghosts.
  // (However this may increase communication costs...)
  : latitude(grid, "lat"),
    longitude(grid, "lon"),
    bed_elevation(grid, "topg"),
    sea_level_elevation(grid, "sea_level"),
    ice_thickness(grid, "thk"),
    ice_area_specific_volume(grid, "ice_area_specific_volume"),
    cell_type(grid, "mask"),
    cell_grounded_fraction(grid, "cell_grounded_fraction"),
    ice_surface_elevation(grid, "usurf") {

  latitude.metadata(0)
      .long_name("latitude")
      .units("degree_north")
      .standard_name("latitude")
      .set_time_independent(true);
  latitude.metadata()["grid_mapping"] = "";
  latitude.metadata()["valid_range"]  = { -90.0, 90.0 };

  longitude.metadata(0)
      .long_name("longitude")
      .units("degree_east")
      .standard_name("longitude")
      .set_time_independent(true);
  longitude.metadata()["grid_mapping"] = "";
  longitude.metadata()["valid_range"]  = { -180.0, 180.0 };

  bed_elevation.metadata(0)
      .long_name("bedrock surface elevation")
      .units("m")
      .standard_name("bedrock_altitude");

  sea_level_elevation.metadata(0)
      .long_name("sea level elevation above reference ellipsoid")
      .units("meters")
      .standard_name("sea_surface_height_above_reference_ellipsoid");

  ice_thickness.metadata(0)
      .long_name("land ice thickness")
      .units("m")
      .standard_name("land_ice_thickness");
  ice_thickness.metadata()["valid_min"] = { 0.0 };

  ice_area_specific_volume.metadata(0)
      .long_name("ice-volume-per-area in partially-filled grid cells")
      .units("m^3/m^2");
  ice_area_specific_volume.metadata()["comment"] =
      "this variable represents the amount of ice in a partially-filled cell and not "
      "the corresponding geometry, so thinking about it as 'thickness' is not helpful";

  cell_type.metadata(0)
      .long_name("ice-type (ice-free/grounded/floating/ocean) integer mask")
      .set_output_type(io::PISM_INT);
  cell_type.metadata()["flag_values"] = { MASK_ICE_FREE_BEDROCK, MASK_GROUNDED, MASK_FLOATING,
                                          MASK_ICE_FREE_OCEAN };
  cell_type.metadata()["flag_meanings"] =
      "ice_free_bedrock grounded_ice floating_ice ice_free_ocean";

  cell_grounded_fraction.metadata(0).long_name(
      "fractional grounded/floating mask (floating=0, grounded=1)");

  ice_surface_elevation.metadata(0)
      .long_name("ice upper surface elevation")
      .units("m")
      .standard_name("surface_altitude");

  // make sure all the fields are initialized
  latitude.set(0.0);
  longitude.set(0.0);
  bed_elevation.set(0.0);
  sea_level_elevation.set(0.0);
  ice_thickness.set(0.0);
  ice_area_specific_volume.set(0.0);
  ensure_consistency(0.0);
}

void Geometry::ensure_consistency(double ice_free_thickness_threshold) {
  auto grid = ice_thickness.grid();
  Config::ConstPtr config = grid->ctx()->config();

  array::AccessScope list{&sea_level_elevation, &bed_elevation,
      &ice_thickness, &ice_area_specific_volume,
      &cell_type, &ice_surface_elevation};

  // first ensure that ice_area_specific_volume is 0 if ice_thickness > 0.
  {
    ParallelSection loop(grid->com);
    try {
      for (auto p = grid->points(); p; p.next()) {
        const int i = p.i(), j = p.j();

        if (ice_thickness(i, j) < 0.0) {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "H = %e (negative) at point i=%d, j=%d",
                                        ice_thickness(i, j), i, j);
        }

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
      for (auto p = grid->points(); p; p.next()) {
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
    ice_density = config->get_number("constants.ice.density"),
    ocean_density = config->get_number("constants.sea_water.density");

  try {
    compute_grounded_cell_fraction(ice_density,
                                   ocean_density,
                                   sea_level_elevation,
                                   ice_thickness,
                                   bed_elevation,
                                   cell_grounded_fraction);
  } catch (RuntimeError &e) {
    e.add_context("computing the grounded cell fraction");

    std::string output_file = config->get_string("output.file");
    std::string o_file = filename_add_suffix(output_file,
                                             "_grounded_cell_fraction_failed", "");
    // save geometry to a file for debugging
    dump(o_file.c_str());
    throw;
  }
}

void Geometry::dump(const char *filename) const {
  auto grid = ice_thickness.grid();

  File file(grid->com, filename,
            string_to_backend(grid->ctx()->config()->get_string("output.format")),
            io::PISM_READWRITE_CLOBBER);

  io::define_time(file, *grid->ctx());
  io::append_time(file, *grid->ctx()->config(), 0.0);

  latitude.write(file);
  longitude.write(file);
  bed_elevation.write(file);
  sea_level_elevation.write(file);
  ice_thickness.write(file);
  ice_area_specific_volume.write(file);
  cell_type.write(file);
  cell_grounded_fraction.write(file);
  ice_surface_elevation.write(file);
}

/*! Compute the elevation of the bottom surface of the ice.
 */
void ice_bottom_surface(const Geometry &geometry, array::Scalar &result) {

  auto grid = result.grid();
  auto config = grid->ctx()->config();

  double
    ice_density   = config->get_number("constants.ice.density"),
    water_density = config->get_number("constants.sea_water.density"),
    alpha         = ice_density / water_density;

  const array::Scalar &ice_thickness = geometry.ice_thickness;
  const array::Scalar &bed_elevation = geometry.bed_elevation;
  const array::Scalar &sea_level     = geometry.sea_level_elevation;

  array::AccessScope list{&ice_thickness, &bed_elevation, &sea_level, &result};

  ParallelSection loop(grid->com);
  try {
    for (auto p = grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      double
        b_grounded = bed_elevation(i, j),
        b_floating = sea_level(i, j) - alpha * ice_thickness(i, j);

      result(i, j) = std::max(b_grounded, b_floating);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result.update_ghosts();
}

//! Computes the ice volume, in m^3.
double ice_volume(const Geometry &geometry, double thickness_threshold) {
  auto grid = geometry.ice_thickness.grid();
  auto config = grid->ctx()->config();

  array::AccessScope list{&geometry.ice_thickness};

  double volume = 0.0;

  auto cell_area = grid->cell_area();

  {
    for (auto p = grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (geometry.ice_thickness(i,j) >= thickness_threshold) {
        volume += geometry.ice_thickness(i,j) * cell_area;
      }
    }
  }

  // Add the volume of the ice in Href:
  if (config->get_flag("geometry.part_grid.enabled")) {
    list.add(geometry.ice_area_specific_volume);
    for (auto p = grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      volume += geometry.ice_area_specific_volume(i,j) * cell_area;
    }
  }

  return GlobalSum(grid->com, volume);
}

double ice_volume_not_displacing_seawater(const Geometry &geometry,
                                          double thickness_threshold) {
  auto grid = geometry.ice_thickness.grid();
  auto config = grid->ctx()->config();

  const double
    sea_water_density = config->get_number("constants.sea_water.density"),
    ice_density       = config->get_number("constants.ice.density"),
    cell_area         = grid->cell_area();

  array::AccessScope list{&geometry.cell_type, &geometry.ice_thickness,
      &geometry.bed_elevation, &geometry.sea_level_elevation};

  double volume = 0.0;

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      bed       = geometry.bed_elevation(i, j),
      thickness = geometry.ice_thickness(i, j),
      sea_level = geometry.sea_level_elevation(i, j);

    if (geometry.cell_type.grounded(i, j) and thickness > thickness_threshold) {
      double max_floating_thickness =
          std::max(sea_level - bed, 0.0) * (sea_water_density / ice_density);
      volume += cell_area * (thickness - max_floating_thickness);
    }
  } // end of the loop over grid points

  return GlobalSum(grid->com, volume);
}

static double compute_area(const Grid &grid, std::function<bool(int, int)> condition) {
  double cell_area = grid.cell_area();
  double area = 0.0;

  for (auto p = grid.points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (condition(i, j)) {
      area += cell_area;
    }
  }

  return GlobalSum(grid.com, area);
}

//! Computes ice area, in m^2.
double ice_area(const Geometry &geometry, double thickness_threshold) {
  array::AccessScope list{ &geometry.ice_thickness };
  return compute_area(*geometry.ice_thickness.grid(), [&](int i, int j) {
    return geometry.ice_thickness(i, j) >= thickness_threshold;
  });
}

//! Computes grounded ice area, in m^2.
double ice_area_grounded(const Geometry &geometry, double thickness_threshold) {
  array::AccessScope list{ &geometry.cell_type, &geometry.ice_thickness };
  return compute_area(*geometry.ice_thickness.grid(), [&](int i, int j) {
    return (geometry.cell_type.grounded(i, j) and
            geometry.ice_thickness(i, j) >= thickness_threshold);
  });
}

//! Computes floating ice area, in m^2.
double ice_area_floating(const Geometry &geometry, double thickness_threshold) {
  array::AccessScope list{ &geometry.cell_type, &geometry.ice_thickness };
  return compute_area(*geometry.ice_thickness.grid(), [&](int i, int j) {
    return (geometry.cell_type.ocean(i, j) and geometry.ice_thickness(i, j) >= thickness_threshold);
  });
}


//! Computes the sea level rise that would result if all the ice were melted.
double sea_level_rise_potential(const Geometry &geometry, double thickness_threshold) {
  auto config = geometry.ice_thickness.grid()->ctx()->config();

  const double
    water_density = config->get_number("constants.fresh_water.density"),
    ice_density   = config->get_number("constants.ice.density"),
    ocean_area    = config->get_number("constants.global_ocean_area");

  const double
    volume                  = ice_volume_not_displacing_seawater(geometry,
                                                                 thickness_threshold),
    additional_water_volume = (ice_density / water_density) * volume,
    sea_level_change        = additional_water_volume / ocean_area;

  return sea_level_change;
}


/*!
 * @brief Set no_model_mask variable to have value 1 in strip of width 'strip' m around
 * edge of computational domain, and value 0 otherwise.
 */
void set_no_model_strip(const Grid &grid, double width, array::Scalar &result) {

  if (width <= 0.0) {
    return;
  }

  array::AccessScope list(result);

  for (auto p = grid.points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (grid::in_null_strip(grid, i, j, width)) {
      result(i, j) = 1;
    } else {
      result(i, j) = 0;
    }
  }

  result.update_ghosts();
}

} // end of namespace pism

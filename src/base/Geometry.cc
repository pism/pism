/* Copyright (C) 2017 PISM Authors
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

#include "base/util/iceModelVec.hh"
#include "base/util/IceModelVec2CellType.hh"
#include "base/util/Mask.hh"
#include "base/grounded_cell_fraction.hh"

namespace pism {

Geometry::Geometry(IceGrid::ConstPtr grid) {
  // FIXME: these fields should be "global", i.e. without ghosts.
  const unsigned int WIDE_STENCIL = grid->ctx()->config()->get_double("grid.max_stencil_width");

  m_cell_area.create(grid, "cell_area", WITHOUT_GHOSTS);
  m_cell_area.set_attrs("diagnostic", "cell areas", "m2", "");
  m_cell_area.metadata().set_string("comment",
                                    "values are equal to dx*dy "
                                    "if projection parameters are not available; "
                                    "otherwise WGS84 ellipsoid is used");
  m_cell_area.set_time_independent(true);
  m_cell_area.metadata().set_string("glaciological_units", "km2");
  m_cell_area.write_in_glaciological_units = true;

  m_latitude.create(grid, "lat", WITH_GHOSTS); // has ghosts so that we can compute cell areas
  m_latitude.set_attrs("mapping", "latitude", "degree_north", "latitude");
  m_latitude.set_time_independent(true);
  m_latitude.metadata().set_string("coordinates", "");
  m_latitude.metadata().set_string("grid_mapping", "");
  m_latitude.metadata().set_double("valid_min", -90.0);
  m_latitude.metadata().set_double("valid_max",  90.0);

  m_longitude.create(grid, "lon", WITH_GHOSTS);
  m_longitude.set_attrs("mapping", "longitude", "degree_east", "longitude");
  m_longitude.set_time_independent(true);
  m_longitude.metadata().set_string("coordinates", "");
  m_longitude.metadata().set_string("grid_mapping", "");
  m_longitude.metadata().set_double("valid_min", -180.0);
  m_longitude.metadata().set_double("valid_max",  180.0);

  m_bed_elevation.create(grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
  m_bed_elevation.set_attrs("model_state", "bedrock surface elevation",
                            "m", "bedrock_altitude");

  m_sea_level_elevation.create(grid, "sea_level_elevation", WITHOUT_GHOSTS);
  m_sea_level_elevation.set_attrs("model_state",
                                  "sea level elevation above reference ellipsoid", "meters",
                                  "sea_surface_height_above_reference_ellipsoid");

  m_ice_thickness.create(grid, "thk", WITH_GHOSTS, WIDE_STENCIL);
  m_ice_thickness.set_attrs("model_state", "land ice thickness",
                            "m", "land_ice_thickness");
  m_ice_thickness.metadata().set_double("valid_min", 0.0);

  m_ice_area_specific_volume.create(grid, "ice_area_specific_volume", WITH_GHOSTS);
  m_ice_area_specific_volume.set_attrs("model_state",
                                       "ice-volume-per-area in partially-filled grid cells",
                                       "m3/m2", "");

  m_cell_type.create(grid, "mask", WITH_GHOSTS, WIDE_STENCIL);
  m_cell_type.set_attrs("diagnostic", "ice-type (ice-free/grounded/floating/ocean) integer mask",
                        "", "");
  std::vector<double> mask_values = {
    MASK_ICE_FREE_BEDROCK,
    MASK_GROUNDED,
    MASK_FLOATING,
    MASK_ICE_FREE_OCEAN};

  m_cell_type.metadata().set_doubles("flag_values", mask_values);
  m_cell_type.metadata().set_string("flag_meanings",
                                    "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");

  m_cell_grounded_fraction.create(grid, "cell_grounded_fraction", WITHOUT_GHOSTS);
  m_cell_grounded_fraction.set_attrs("internal",
                                     "fractional grounded/floating mask (floating=0, grounded=1)",
                                     "", "");

  m_ice_surface_elevation.create(grid, "usurf", WITH_GHOSTS, WIDE_STENCIL);
  m_ice_surface_elevation.set_attrs("diagnostic", "ice upper surface elevation",
                                    "m", "surface_altitude");
}

void Geometry::ensure_consistency(double ice_free_thickness_threshold) {
  IceGrid::ConstPtr grid = m_ice_thickness.get_grid();
  Config::ConstPtr config = grid->ctx()->config();

  IceModelVec::AccessList list{&m_sea_level_elevation, &m_bed_elevation,
      &m_ice_thickness, &m_ice_area_specific_volume,
      &m_cell_type, &m_ice_surface_elevation};

  // first ensure that ice_area_specific_volume is 0 if ice_thickness > 0.
  {
    ParallelSection loop(grid->com);
    try {
      for (Points p(*grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        if (m_ice_thickness(i, j) > 0.0 and m_ice_area_specific_volume(i, j) > 0.0) {
          m_ice_thickness(i, j) += m_ice_area_specific_volume(i, j);
          m_ice_area_specific_volume(i, j) = 0.0;
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
        gc.compute(m_sea_level_elevation(i, j), m_bed_elevation(i, j), m_ice_thickness(i, j),
                   &mask, &m_ice_surface_elevation(i, j));
        m_cell_type(i, j) = mask;
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();
  }

  m_ice_thickness.update_ghosts();
  m_ice_area_specific_volume.update_ghosts();
  m_cell_type.update_ghosts();
  m_ice_surface_elevation.update_ghosts();

  const double
    ice_density = config->get_double("constants.ice.density"),
    ocean_density = config->get_double("constants.sea_water.density"),
    sea_level = 0;              // FIXME: use the 2D field

  compute_grounded_cell_fraction(ice_density,
                                 ocean_density,
                                 sea_level,
                                 m_ice_thickness,
                                 m_bed_elevation,
                                 m_cell_type,
                                 m_cell_grounded_fraction,
                                 NULL, NULL);
}

IceModelVec2S& Geometry::cell_area() {
  return m_cell_area;
}

IceModelVec2S& Geometry::latitude() {
  return m_latitude;
}

IceModelVec2S& Geometry::longitude() {
  return m_longitude;
}

IceModelVec2S& Geometry::bed_elevation() {
  return m_bed_elevation;
}

IceModelVec2S& Geometry::sea_level_elevation() {
  return m_sea_level_elevation;
}

IceModelVec2S& Geometry::ice_thickness() {
  return m_ice_thickness;
}

IceModelVec2S& Geometry::ice_area_specific_volume() {
  return m_ice_area_specific_volume;
}

IceModelVec2S& Geometry::cell_grounded_fraction() {
  return m_cell_grounded_fraction;
}

IceModelVec2S& Geometry::ice_surface_elevation() {
  return m_ice_surface_elevation;
}

IceModelVec2CellType& Geometry::cell_type() {
  return m_cell_type;
}

} // end of namespace pism

/* Copyright (C) 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023 PISM Authors
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

#ifndef PISM_GEOMETRY_H
#define PISM_GEOMETRY_H

#include "pism/util/array/CellType.hh"

namespace pism {

class Grid;

class Geometry {
public:
  Geometry(const std::shared_ptr<const Grid> &grid);

  /*!
   * Ensures consistency of ice geometry by re-computing cell type, cell grounded fraction, and ice
   * surface elevation.
   */
  void ensure_consistency(double ice_free_thickness_threshold);

  // This is grid information, which is not (strictly speaking) ice geometry, but it should be
  // available everywhere we use ice geometry.
  array::Scalar latitude;
  array::Scalar longitude;

  // Part of ice geometry, but managed by the bed model and the ocean model. From the point of view
  // of the code updating ice geometry, these are inputs. These fields should be filled in before
  // passing a Geometry instance to the code that uses it.
  array::Scalar2 bed_elevation;
  array::Scalar1 sea_level_elevation;

  // the minimal "state"
  array::Scalar2 ice_thickness;
  array::Scalar1 ice_area_specific_volume;

  // redundant fields (can be computed using the ones above)
  array::CellType2 cell_type;
  array::Scalar cell_grounded_fraction;
  array::Scalar2 ice_surface_elevation;

  void dump(const char *filename) const;
};

void ice_bottom_surface(const Geometry &geometry, array::Scalar &result);

double ice_volume(const Geometry &geometry, double thickness_threshold);
double ice_area_floating(const Geometry &geometry, double thickness_threshold);
double ice_area_grounded(const Geometry &geometry, double thickness_threshold);
double ice_area(const Geometry &geometry, double thickness_threshold);
double ice_volume_not_displacing_seawater(const Geometry &geometry,
                                          double thickness_threshold);
double sea_level_rise_potential(const Geometry &geometry, double thickness_threshold);

void set_no_model_strip(const Grid &grid, double width, array::Scalar &result);

} // end of namespace pism

#endif /* PISM_GEOMETRY_H */

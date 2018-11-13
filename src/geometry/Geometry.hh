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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {

class Geometry {
public:
  Geometry(IceGrid::ConstPtr grid);

  /*!
   * Ensures consistency of ice geometry by re-computing cell type, cell grounded fraction, and ice
   * surface elevation.
   */
  void ensure_consistency(double ice_free_thickness_threshold);

  // This is grid information, which is not (strictly speaking) ice geometry, but it should be
  // available everywhere we use ice geometry.
  IceModelVec2S latitude;
  IceModelVec2S longitude;

  // Part of ice geometry, but managed by the bed model and the ocean model. From the point of view
  // of the code updating ice geometry, these are inputs. These fields should be filled in before
  // passing a Geometry instance to the code that uses it.
  IceModelVec2S bed_elevation;
  IceModelVec2S sea_level_elevation;

  // the minimal "state"
  IceModelVec2S ice_thickness;
  IceModelVec2S ice_area_specific_volume; // previously known as Href

  // redundant fields (can be computed using the ones above)
  IceModelVec2CellType cell_type;
  IceModelVec2S cell_grounded_fraction;
  IceModelVec2S ice_surface_elevation;
};

} // end of namespace pism

#endif /* GEOMETRY_H */

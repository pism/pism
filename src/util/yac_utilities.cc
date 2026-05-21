/* Copyright (C) 2026 PISM Authors
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

#include "pism/util/yac_utilities.hh"
#include "pism/util/LonLatGrid.hh"
#include "pism/util/error_handling.hh"

#if (Pism_USE_YAC == 0)
#error "This code requires YAC"
#endif

extern "C" {
#include "yac.h"
}

namespace pism {

/*!
 * Define the PISM grid. Each PE defines its own subdomain.
 *
 * Returns the point ID that can be used to define a "field".
 */
int define_yac_grid(const std::vector<double> &x_cell, const std::vector<double> &y_cell,
                    const std::string &grid_name, const std::string &projection) {

  if (projection.empty()) {
    throw pism::RuntimeError::formatted(
        PISM_ERROR_LOCATION, "grid '%s' has no projection information", grid_name.c_str());
  }

  // Shift x and y by half a grid spacing and add one more row and column to get
  // coordinates of corners of cells in the local sub-domain:
  std::vector<double> x_node(x_cell.size() + 1), y_node(y_cell.size() + 1);
  {
    // note: dx and dy may be negative here
    double dx = x_cell[1] - x_cell[0];
    double dy = y_cell[1] - y_cell[0];

    for (size_t k = 0; k < x_cell.size(); ++k) {
      x_node[k] = x_cell[k] - 0.5 * dx;
    }
    x_node.back() = x_cell.back() + 0.5 * dx;

    for (size_t k = 0; k < y_cell.size(); ++k) {
      y_node[k] = y_cell[k] - 0.5 * dy;
    }
    y_node.back() = y_cell.back() + 0.5 * dy;
  }

  // Compute lon,lat coordinates of cell centers:
  LonLatGrid cells(x_cell, y_cell, projection);
  // Compute lon,lat coordinates of cell corners:
  LonLatGrid nodes(x_node, y_node, projection);

  int point_id = 0;
  {
    int cyclic[] = { 0, 0 };

    int grid_id = 0;

    int n_nodes[2] = { (int)x_node.size(), (int)y_node.size() };
    yac_cdef_grid_curve2d(grid_name.c_str(), n_nodes, cyclic, nodes.lon.data(), nodes.lat.data(),
                          &grid_id);

    int n_cells[2] = { (int)x_cell.size(), (int)y_cell.size() };
    yac_cdef_points_curve2d(grid_id, n_cells, YAC_LOCATION_CELL, cells.lon.data(), cells.lat.data(),
                            &point_id);
  }
  return point_id;
}
} // namespace pism

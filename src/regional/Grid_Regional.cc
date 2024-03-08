/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2023, 2024 PISM Authors
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

#include <vector>
#include <cassert>

#include <gsl/gsl_interp.h>

#include "pism/util/pism_options.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Grid.hh"
#include "pism/util/io/File.hh"
#include "pism/util/Component.hh" // process_input_options
#include "pism/util/Context.hh"

namespace pism {

static void validate_range(const std::string &axis,
                           const std::vector<double> &x,
                           double x_min, double x_max) {
  if (x_min >= x_max) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid -%s_range: %s_min >= %s_max (%f >= %f).",
                                  axis.c_str(), axis.c_str(), axis.c_str(),
                                  x_min, x_max);
  }

  if (x_min >= x.back()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid -%s_range: %s_min >= max(%s) (%f >= %f).",
                                  axis.c_str(), axis.c_str(), axis.c_str(),
                                  x_min, x.back());
  }

  if (x_max <= x.front()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid -%s-range: %s_max <= min(%s) (%f <= %f).",
                                  axis.c_str(), axis.c_str(), axis.c_str(),
                                  x_max, x.front());
  }
}

static void subset_extent(const std::string& axis,
                          const std::vector<double> &x,
                          double x_min, double x_max,
                          double &x0, double &Lx, unsigned int &Mx) {

  validate_range(axis, x, x_min, x_max);

  size_t x_start = gsl_interp_bsearch(x.data(), x_min, 0, x.size() - 1);
  // include one more point if we can
  if (x_start > 0) {
    x_start -= 1;
  }

  size_t x_end = gsl_interp_bsearch(x.data(), x_max, 0, x.size() - 1);
  // include one more point if we can
  x_end = std::min(x.size() - 1, x_end + 1);

  // NOTE: this assumes the CELL_CORNER grid registration
  Lx = (x[x_end] - x[x_start]) / 2.0;

  x0 = (x[x_start] + x[x_end]) / 2.0;

  assert(x_end + 1 >= x_start);
  Mx = x_end + 1 - x_start;
}

//! Create a grid using command-line options and (possibly) an input file.
/** Processes options -i, -bootstrap, -Mx, -My, -Mz, -Lx, -Ly, -Lz, -x_range, -y_range.
 */
std::shared_ptr<Grid> regional_grid_from_options(std::shared_ptr<Context> ctx) {

  auto options = process_input_options(ctx->com(), ctx->config());

  const options::RealList x_range("-x_range",
                                  "range of X coordinates in the selected subset", {});
  const options::RealList y_range("-y_range",
                                  "range of Y coordinates in the selected subset", {});

  const options::Integer refinement_factor("-refinement_factor",
                                           "Grid refinement factor (applies to the horizontal grid)", 1);

  if (options.type == INIT_BOOTSTRAP and x_range.is_set() and y_range.is_set()) {
    // bootstrapping; get domain size defaults from an input file, allow overriding all grid
    // parameters using command-line options

    if (x_range->size() != 2) {
      throw RuntimeError(PISM_ERROR_LOCATION, "invalid -x_range argument: need 2 numbers.");
    }

    if (y_range->size() != 2) {
      throw RuntimeError(PISM_ERROR_LOCATION, "invalid -y_range argument: need 2 numbers.");
    }

    grid::Parameters input_grid(*ctx->config());

    std::vector<std::string> names = {"enthalpy", "temp", "land_ice_thickness",
                                      "bedrock_altitude", "thk", "topg"};
    bool grid_info_found = false;

    File file(ctx->com(), options.filename, io::PISM_NETCDF3, io::PISM_READONLY);
    for (const auto& name : names) {

      grid_info_found = file.variable_exists(name);
      if (not grid_info_found) {
        // Failed to find using a short name. Try using name as a
        // standard name...
        grid_info_found = file.find_variable("unlikely_name", name).exists;
      }

      if (grid_info_found) {
        input_grid = grid::Parameters(*ctx, file, name, grid::CELL_CORNER);

        auto full_grid = grid::InputGridInfo(file, name, ctx->unit_system(), grid::CELL_CORNER);

        // x direction
        subset_extent("x", full_grid.x, x_range[0], x_range[1],
                      input_grid.x0, input_grid.Lx, input_grid.Mx);
        // y direction
        subset_extent("y", full_grid.y, y_range[0], y_range[1],
                      input_grid.y0, input_grid.Ly, input_grid.My);

        // Set registration to "CELL_CORNER" so that Grid computes
        // coordinates correctly.
        input_grid.registration = grid::CELL_CORNER;

        break;
      }
    }

    if (not grid_info_found) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "no geometry information found in '%s'",
                                    options.filename.c_str());
    }

    if (refinement_factor > 1) {
      input_grid.Mx = (input_grid.Mx - 1) * refinement_factor + 1;
      input_grid.My = (input_grid.My - 1) * refinement_factor + 1;
    }

    // process options controlling vertical grid parameters, overriding values read from a file
    input_grid.vertical_grid_from_options(ctx->config());

    // process options controlling ownership ranges
    input_grid.ownership_ranges_from_options(ctx->size());

    return std::make_shared<Grid>(ctx, input_grid);
  }

  if (x_range.is_set() ^ y_range.is_set()) {
    throw RuntimeError(PISM_ERROR_LOCATION, "Please set both -x_range and -y_range.");
  }

  return Grid::FromOptions(ctx);
}


} // end of namespace pism

/* Copyright (C) 2015 PISM Authors
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

#include "base/util/error_handling.hh"
#include "base/util/IceGrid.hh"
#include "base/util/io/PIO.hh"

namespace pism {

static void subset_parameters(const std::vector<double> &x,
                              double x_min, double x_max,
                              double &x0, double &Lx, unsigned int &Mx) {

  assert(x_max > x_min);

  size_t x_start = gsl_interp_bsearch(&x[0], x_min, 0, x.size() - 1);
  // include one more point if we can
  if (x_start > 0) {
    x_start -= 1;
  }

  size_t x_end = gsl_interp_bsearch(&x[0], x_max, 0, x.size() - 1);
  // include one more point if we can
  x_end = std::min(x.size() - 1, x_end + 1);

  Lx = (x[x_end] - x[x_start]) / 2.0;

  x0 = (x[x_start] + x[x_end]) / 2.0;

  assert(x_end + 1 >= x_start);
  Mx = x_end + 1 - x_start;
}

static void validate_range(const std::string &axis, std::vector<double> &x,
                           double x_min, double x_max) {
  if (x_min >= x_max) {
    throw RuntimeError::formatted("invalid -%s_range: %s_min >= %s_max (%f >= %f).",
                                  axis.c_str(), axis.c_str(), axis.c_str(),
                                  x_min, x_max);
  }

  if (x_min >= x.back()) {
    throw RuntimeError::formatted("invalid -%s_range: %s_min >= max(%s) (%f >= %f).",
                                  axis.c_str(), axis.c_str(), axis.c_str(),
                                  x_min, x.back());
  }

  if (x_max <= x.front()) {
    throw RuntimeError::formatted("invalid -%s-range: %s_max <= min(%s) (%f <= %f).",
                                  axis.c_str(), axis.c_str(), axis.c_str(),
                                  x_max, x.front());
  }
}

IceGrid::Ptr regional_grid(Context::Ptr ctx,
                           const std::string &filename,
                           double x_min, double x_max,
                           double y_min, double y_max) {

  // FIXME: we should add periodicity to x and y coordinate variables in PISM output files.
  Periodicity p = string_to_periodicity(ctx->config()->get_string("grid_periodicity"));

  GridParameters params;

  PIO file(ctx->com(), "netcdf3");
  file.open(filename, PISM_READONLY); // will be closed automatically

  grid_info full;
  if (file.inq_var("enthalpy")) {
    full = grid_info(file, "enthalpy", ctx->unit_system(), p);
  } else if (file.inq_var("temp")) {
    full = grid_info(file, "temp", ctx->unit_system(), p);
  } else {
    throw RuntimeError::formatted("neither 'enthalpy' nor 'temp' was found in %s.",
                                  filename.c_str());
  }

  // x direction
  validate_range("x", full.x, x_min, x_max);
  subset_parameters(full.x, x_min, x_max, params.x0, params.Lx, params.Mx);
  // y direction
  validate_range("y", full.y, y_min, y_max);
  subset_parameters(full.y, y_min, y_max, params.y0, params.Ly, params.My);

  // vertical grid parameters
  params.z = full.z;
  params.periodicity = NONE;

  // set ownership ranges
  params.ownership_ranges_from_options(ctx->size());

  return IceGrid::Ptr(new IceGrid(ctx, params));
}

} // end of namespace pism

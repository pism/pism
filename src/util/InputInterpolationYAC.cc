/* Copyright (C) 2024 PISM Authors
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

#include <cstddef>
#include <memory>
#include <vector>
#include <cmath>

#include "InputInterpolation.hh"
#include "Interpolation1D.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/projection.hh"
#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/InputInterpolationYAC.hh"
#include "pism/util/pism_utilities.hh" // GlobalMin()

#if (Pism_USE_PROJ == 0)
#error "This code requires PROJ"
#endif

#include <proj.h>

#include "pism/util/Proj.hh"

#if (Pism_USE_YAC_INTERPOLATION == 0)
#error "This code requires YAC"
#endif

extern "C" {
#include "yac.h"
}

namespace pism {

/*!
 * Utility class converting `x,y` coordinates in a projection to a `lon,lat` pair.
 *
 * Requires the `PROJ` library.
 */
class LonLatCalculator {
public:
  LonLatCalculator(const std::string &proj_string)
      : m_coordinate_mapping(proj_string, "EPSG:4326") {
  }

  std::array<double, 2> lonlat(double x, double y) {
    PJ_COORD in, out;

    in.xy = { x, y };
    out   = proj_trans(m_coordinate_mapping, PJ_FWD, in);

    return { out.lp.phi, out.lp.lam };
  }

private:
  Proj m_coordinate_mapping;
};

/*!
 * Grid definition using coordinates in radians.
 */
struct LonLatGrid {
  std::vector<double> lon;
  std::vector<double> lat;

  /*!
   *
   * Converts a Cartesian grid in a `projection` that uses coordinates
   * `x` and `y` in meters into the form that can be used to define a
   * curvilinear grid in YAC.
   *
   * The `projection` string has to use the format compatible with PROJ.
   */
  LonLatGrid(const std::vector<double> &x, const std::vector<double> &y,
             const std::string &projection) {

    int nrow = y.size();
    int ncol = x.size();
    int N    = nrow * ncol;

    lon.resize(N);
    lat.resize(N);

    // convert from (row, col) to the linear index in "cell" arrays
    auto C = [ncol](int row, int col) { return row * ncol + col; };

    // convert from degrees to radians
    auto deg2rad = [](double degree) { return degree * M_PI / 180; };

    pism::LonLatCalculator mapping(projection);

    for (int row = 0; row < nrow; ++row) {
      for (int col = 0; col < ncol; ++col) {
        auto coords = mapping.lonlat(x[col], y[row]);

        lon[C(row, col)] = deg2rad(coords[0]);
        lat[C(row, col)] = deg2rad(coords[1]);
      }
    }
  }
};


/*!
 * Define the PISM grid. Each PE defines its own subdomain.
 *
 * Returns the point ID that can be used to define a "field".
 */
int InputInterpolationYAC::define_grid(const std::vector<double> &x_cell,
                                       const std::vector<double> &y_cell,
                                       const std::string &grid_name,
                                       const std::string &projection) {

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

/*!
 * Return the YAC field_id corresponding to a given PISM grid and its
 * projection info.
 *
 * @param[in] component_id YAC component ID
 * @param[in] pism_grid PISM's grid
 * @param[in] name string describing this grid and field
 */
int InputInterpolationYAC::define_field(int component_id, const std::vector<double> &x,
                                        const std::vector<double> &y,
                                        const std::string &proj_string, const std::string &name) {

  int point_id = define_grid(x, y, name, proj_string);

  const char *time_step_length = "1";
  const int point_set_size     = 1;
  const int collection_size    = 1;

  int field_id = 0;
  yac_cdef_field(name.c_str(), component_id, &point_id, point_set_size, collection_size,
                 time_step_length, YAC_TIME_UNIT_SECOND, &field_id);
  return field_id;
}

static void pism_yac_error_handler(MPI_Comm /* unused */, const char *msg, const char *source,
                                   int line) {
  throw pism::RuntimeError::formatted(pism::ErrorLocation(source, line), "YAC error: %s", msg);
}

/*!
 * Extract the "local" (corresponding to the current sub-domain) grid subset.
 */
static std::vector<double> grid_subset(int xs, int xm, const std::vector<double> &coords) {
  std::vector<double> result(xm);
  for (int k = 0; k < xm; ++k) {
    result[k] = coords[xs + k];
  }

  return result;
}

static double dx_estimate(Proj &mapping, double x1, double x2, double y) {
  PJ_COORD p1, p2;
  p1.lp = {proj_torad(x1), proj_torad(y)};
  p2.lp = {proj_torad(x2), proj_torad(y)};

  return proj_lp_dist(mapping, p1, p2);
}

static double dx_min(const std::string &proj_string,
                     const std::vector<double> &x,
                     const std::vector<double> &y) {
  size_t Nx = x.size();
  size_t Ny = y.size();

  Proj mapping(proj_string, "EPSG:4326");

  double dx = dx_estimate(mapping, x[0], x[1], y[0]);
  for (size_t j = 0; j < Ny; ++j) {   // y
    for (size_t i = 1; i < Nx; ++i) { // x; note: starts from 1
      dx = std::min(dx, dx_estimate(mapping, x[i-1], x[i], y[j]));
    }
  }
  return dx;
}

static double dy_estimate(Proj &mapping, double x, double y1, double y2) {
  PJ_COORD p1, p2;
  p1.lp = {proj_torad(x), proj_torad(y1)};
  p2.lp = {proj_torad(x), proj_torad(y2)};

  return proj_lp_dist(mapping, p1, p2);
}

static double dy_min(const std::string &proj_string,
                     const std::vector<double> &x,
                     const std::vector<double> &y) {
  size_t Nx = x.size();
  size_t Ny = y.size();

  Proj mapping(proj_string, "EPSG:4326");

  double dy = dy_estimate(mapping, x[0], x[1], y[0]);
  for (size_t i = 0; i < Nx; ++i) { // x
    for (size_t j = 1; j < Ny; ++j) { // y; note: starts from 1
      dy = std::min(dy, dy_estimate(mapping, x[i], y[j-1], y[j]));
    }
  }
  return dy;
}

InputInterpolationYAC::InputInterpolationYAC(const pism::Grid &target_grid,
                                             const pism::File &input_file,
                                             const std::string &variable_name,
                                             InterpolationType type)
    : m_instance_id(0), m_source_field_id(0), m_target_field_id(0) {
  auto ctx = target_grid.ctx();

  if (target_grid.get_mapping_info().proj_string.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "internal grid projection is not known");
  }

  yac_set_abort_handler((yac_abort_func)pism_yac_error_handler);

  try {
    auto log = ctx->log();

    auto source_grid_name = grid_name(input_file, variable_name, ctx->unit_system());

    log->message(
        2, "* Initializing 2D interpolation on the sphere from '%s' to the internal grid...\n",
        source_grid_name.c_str());

    grid::InputGridInfo source_grid_info(input_file, variable_name, ctx->unit_system(),
                                         pism::grid::CELL_CENTER);

    auto source_grid_mapping = MappingInfo::FromFile(input_file, variable_name, ctx->unit_system());

    std::string grid_mapping_name = source_grid_mapping.cf_mapping["grid_mapping_name"];

    log->message(2, "Input grid:\n");
    if (not grid_mapping_name.empty()) {
      log->message(2, " Grid mapping: %s\n", grid_mapping_name.c_str());
    }
    if (not source_grid_mapping.proj_string.empty()) {
      log->message(2, " PROJ string: '%s'\n", source_grid_mapping.proj_string.c_str());
    }
    source_grid_info.report(*log, 2, ctx->unit_system());

    grid::Parameters P(*ctx->config(), source_grid_info.x.size(), source_grid_info.y.size(),
                       source_grid_info.Lx, source_grid_info.Ly);

    P.x0 = source_grid_info.x0;
    P.y0 = source_grid_info.y0;
    P.registration = grid::CELL_CENTER;
    P.variable_name = variable_name;
    P.vertical_grid_from_options(*ctx->config());
    P.ownership_ranges_from_options(*ctx->config(), ctx->size());

    auto source_grid = std::make_shared<Grid>(ctx, P);

    source_grid->set_mapping_info(source_grid_mapping);

    if (source_grid->get_mapping_info().proj_string.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "unsupported or missing projection info for the grid '%s'",
                                    source_grid_name.c_str());
    }

    m_buffer = std::make_shared<pism::array::Scalar>(source_grid, variable_name);

    std::string target_grid_name = "internal for " + source_grid_name;
    double fill_value            = NAN;
    {
      // Initialize YAC:
      {
        yac_cinit_instance(&m_instance_id);
        yac_cdef_calendar(YAC_YEAR_OF_365_DAYS);
        // Note: zero-padding of months and days *is* required.
        yac_cdef_datetime_instance(m_instance_id, "-1-01-01", "+1-01-01");
      }

      // Define components: this has to be done using *one* call
      // (cannot call yac_cdef_comp?_instance() more than once)
      const int n_comps               = 2;
      const char *comp_names[n_comps] = { "source_component", "target_component" };
      int comp_ids[n_comps]           = { 0, 0 };
      yac_cdef_comps_instance(m_instance_id, comp_names, n_comps, comp_ids);

      double source_resolution = 0.0;
      log->message(2, "Defining the source grid (%s)...\n", source_grid_name.c_str());
      {
        auto x = grid_subset(source_grid->xs(), source_grid->xm(), source_grid_info.x);
        auto y = grid_subset(source_grid->ys(), source_grid->ym(), source_grid_info.y);

        m_source_field_id = define_field(
            comp_ids[0], x, y, source_grid->get_mapping_info().proj_string, source_grid_name);

        double dx = 0.0;
        double dy = 0.0;
        if (source_grid_info.longitude_latitude) {
          dx = dx_min(source_grid_mapping.proj_string, x, y);
          dy = dy_min(source_grid_mapping.proj_string, x, y);
        } else {
          dx = std::abs(x[1] - x[0]);
          dy = std::abs(y[1] - y[0]);
        }
        source_resolution = GlobalMin(ctx->com(), std::min(dx, dy));
        log->message(2, "Source grid resolution: ~%3.3f m\n", source_resolution);
      }

      log->message(2, "Defining the target grid (%s)...\n", target_grid_name.c_str());
      {
        auto x = grid_subset(target_grid.xs(), target_grid.xm(), target_grid.x());
        auto y = grid_subset(target_grid.ys(), target_grid.ym(), target_grid.y());

        m_target_field_id = define_field(
            comp_ids[1], x, y, target_grid.get_mapping_info().proj_string, target_grid_name);
      }

      // Define the interpolation stack:
      {
        std::string method = "nearest neighbor";
        int interp_stack_id = 0;
        yac_cget_interp_stack_config(&interp_stack_id);

        if (type != PIECEWISE_CONSTANT) {
          method = "2nd order conservative";

          int order                = 2;
          int enforce_conservation = 0;
          int partial_coverage     = 0;

          yac_cadd_interp_stack_config_conservative(interp_stack_id, order, enforce_conservation,
                                                    partial_coverage, YAC_CONSERV_DESTAREA);

          // use average over source grid nodes containing a target point as a backup:
          yac_cadd_interp_stack_config_average(interp_stack_id, YAC_AVG_BARY, partial_coverage);
        }

        // use nearest neighbor interpolation as a backup and to interpolate integer
        // fields:
        {
          // nearest neighbor
          int n_neighbors            = 1;
          double scaling             = 1.0;
          double max_search_distance = 0.0; // unlimited
          yac_cadd_interp_stack_config_nnn(interp_stack_id, YAC_NNN_DIST, n_neighbors,
                                           max_search_distance, scaling);
        }

        log->message(2, "Interpolation method: %s\n", method.c_str());

        // last resort: fill with `fill_value`
        yac_cadd_interp_stack_config_fixed(interp_stack_id, fill_value);

        // Define the coupling between fields:
        const int src_lag = 0;
        const int tgt_lag = 0;
        yac_cdef_couple_instance(m_instance_id,
                                 "source_component",       // source component name
                                 source_grid_name.c_str(), // source grid name
                                 source_grid_name.c_str(), // source field name
                                 "target_component",       // target component name
                                 target_grid_name.c_str(), // target grid name
                                 target_grid_name.c_str(), // target field name
                                 "1",                      // time step length in units below
                                 YAC_TIME_UNIT_SECOND,     // time step length units
                                 YAC_REDUCTION_TIME_NONE,  // reduction in time (for
                                                           // asynchronous coupling)
                                 interp_stack_id, src_lag, tgt_lag);

        // free the interpolation stack config now that we defined the coupling
        yac_cfree_interp_stack_config(interp_stack_id);
      }

      double start = MPI_Wtime();
      yac_cenddef_instance(m_instance_id);
      double end = MPI_Wtime();
      log->message(2, "Initialized interpolation from %s in %f seconds.\n",
                   source_grid_name.c_str(), end - start);
    }
  } catch (pism::RuntimeError &e) {
    e.add_context("initializing interpolation from %s to the internal grid",
                  input_file.name().c_str());
    throw;
  }
}

InputInterpolationYAC::~InputInterpolationYAC() {
  yac_ccleanup_instance(m_instance_id);
}

double InputInterpolationYAC::interpolate(const pism::array::Scalar &source,
                                          pism::petsc::Vec &target) const {

  pism::petsc::VecArray input_array(source.vec());
  pism::petsc::VecArray output_array(target);

  double *send_field_    = input_array.get();
  double **send_field[1] = { &send_field_ };

  double *recv_field[1] = { output_array.get() };

  int ierror          = 0;
  int send_info       = 0;
  int recv_info       = 0;
  int collection_size = 1;
  double start        = MPI_Wtime();
  yac_cexchange(m_source_field_id, m_target_field_id, collection_size, send_field, recv_field,
                &send_info, &recv_info, &ierror);
  double end = MPI_Wtime();

  return end - start;
}


void InputInterpolationYAC::regrid(const pism::File &file, pism::array::Scalar &output) const {

  double time_spent =
      InputInterpolation::regrid(output.metadata(0), file, -1, *output.grid(), output.vec());

  auto log = output.grid()->ctx()->log();

  log->message(2, "Interpolation took %f seconds.\n", time_spent);
}


double InputInterpolationYAC::regrid_impl(const SpatialVariableMetadata &metadata,
                                          const pism::File &file, int record_index,
                                          const Grid & /* target_grid (unused) */,
                                          petsc::Vec &output) const {

  // set metadata to help the following call find the variable, convert units, etc
  m_buffer->metadata(0) = metadata;
  m_buffer->read(file, record_index);

  double time_spent = interpolate(*m_buffer, output);

  return time_spent;
}

} // namespace pism

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
#include "pism/util/VariableMetadata.hh"
#include "pism/util/projection.hh"
#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/InputInterpolationYAC.hh"

#if (Pism_USE_PROJ == 0)
#error "This code requires PROJ"
#endif

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
    out   = proj_trans(*m_coordinate_mapping, PJ_FWD, in);

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
int InputInterpolationYAC::define_grid(const pism::Grid &grid, const std::string &grid_name,
                                       const std::string &projection) {

  if (projection.empty()) {
    throw pism::RuntimeError::formatted(
        PISM_ERROR_LOCATION, "grid '%s' has no projection information", grid_name.c_str());
  }

  int xs = grid.xs();
  int ys = grid.ys();
  int xm = grid.xm();
  int ym = grid.ym();

  std::vector<double> x(xm);
  std::vector<double> y(ym);

  // Set x and y to coordinates of centers of cells in the local sub-domain:
  {
    for (int k = 0; k < xm; ++k) {
      x[k] = grid.x(xs + k);
    }
    for (int k = 0; k < ym; ++k) {
      y[k] = grid.y(ys + k);
    }
  }

  // Compute lon,lat coordinates of cell centers:
  LonLatGrid cells(x, y, projection);

  // Shift x and y by half a grid spacing and add one more row and column to get
  // coordinates of corners of cells in the local sub-domain:
  {
    double dx = x[1] - x[0];
    double dy = y[1] - y[0];

    double x_last = x.back() + 0.5 * dx;
    for (size_t k = 0; k < x.size(); ++k) {
      x[k] -= 0.5 * dx;
    }
    x.push_back(x_last);

    double y_last = y.back() + 0.5 * dy;
    for (size_t k = 0; k < y.size(); ++k) {
      y[k] -= 0.5 * dy;
    }
    y.push_back(y_last);
  }

  // Compute lon,lat coordinates of cell corners:
  LonLatGrid nodes(x, y, projection);
  int n_nodes[2] = { (int)x.size(), (int)y.size() };

  std::vector<int> cell_global_index(xm * ym);
  {
    int Mx = grid.Mx();
    int k  = 0;
    for (int j = ys; j < ys + ym; ++j) {
      for (int i = xs; i < xs + xm; ++i) {
        cell_global_index[k] = j * Mx + i;
        ++k;
      }
    }
  }

  int point_id = 0;
  {
    int cyclic[] = { 0, 0 };

    int grid_id = 0;

    yac_cdef_grid_curve2d(grid_name.c_str(), n_nodes, cyclic, nodes.lon.data(), nodes.lat.data(),
                          &grid_id);

    yac_cset_global_index(cell_global_index.data(), YAC_LOCATION_CELL, grid_id);

    int n_cells[2] = { xm, ym };
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
int InputInterpolationYAC::define_field(int component_id, const pism::Grid &grid,
                                        const std::string &name) {

  int point_id = define_grid(grid, name, grid.get_mapping_info().proj_string);

  const char *time_step_length = "1";
  const int point_set_size     = 1;
  const int collection_size    = 1;

  int field_id = 0;
  yac_cdef_field(name.c_str(), component_id, &point_id, point_set_size, collection_size,
                 time_step_length, YAC_TIME_UNIT_SECOND, &field_id);
  return field_id;
}

int InputInterpolationYAC::interpolation_coarse_to_fine(double missing_value) {
  int id = 0;
  yac_cget_interp_stack_config(&id);

  int partial_coverage = 0;

  // average over source grid nodes containing a target point, weighted using
  // barycentric local coordinates
  yac_cadd_interp_stack_config_average(id, YAC_AVG_BARY, partial_coverage);

  // nearest neighbor
  int n_neighbors = 1;
  double scaling  = 1.0;
  double max_search_distance = 0.0; // unlimited
  yac_cadd_interp_stack_config_nnn(id, YAC_NNN_DIST, n_neighbors, max_search_distance, scaling);

  // constant if all of the above failed
  yac_cadd_interp_stack_config_fixed(id, missing_value);

  return id;
}

int InputInterpolationYAC::interpolation_fine_to_coarse(double missing_value) {
  int id = 0;
  yac_cget_interp_stack_config(&id);

  int order                = 1;
  int enforce_conservation = 1;
  int partial_coverage     = 0;

  // conservative
  yac_cadd_interp_stack_config_conservative(id, order, enforce_conservation, partial_coverage,
                                            YAC_CONSERV_DESTAREA);

  // average over source grid nodes containing a target point, weighted using
  // barycentric local coordinates
  yac_cadd_interp_stack_config_average(id, YAC_AVG_BARY, partial_coverage);

  // nearest neighbor
  int n_neighbors = 1;
  double scaling  = 1.0;
  double max_search_distance = 0.0; // unlimited
  yac_cadd_interp_stack_config_nnn(id, YAC_NNN_DIST, n_neighbors, max_search_distance, scaling);

  // constant if all of the above failed
  yac_cadd_interp_stack_config_fixed(id, missing_value);

  return id;
}

static void pism_yac_error_handler(MPI_Comm /* unused */, const char *msg, const char *source,
                                   int line) {
  throw pism::RuntimeError::formatted(pism::ErrorLocation(source, line), "YAC error: %s", msg);
}

InputInterpolationYAC::InputInterpolationYAC(const pism::Grid &target_grid,
                                             const pism::File &input_file,
                                             const std::string &variable_name) {
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

    grid::InputGridInfo info(input_file, variable_name, ctx->unit_system(),
                             pism::grid::CELL_CENTER);

    auto mapping = MappingInfo::FromFile(input_file, variable_name, ctx->unit_system());

    std::string grid_mapping_name = mapping.cf_mapping["grid_mapping_name"];

    log->message(2, "Input grid:\n");
    if (not grid_mapping_name.empty()) {
      log->message(2, " Grid mapping: %s\n", grid_mapping_name.c_str());
    }
    if (not mapping.proj_string.empty()) {
      log->message(2, " PROJ string: '%s'\n", mapping.proj_string.c_str());
    }
    info.report(*log, 2, ctx->unit_system());

    grid::Parameters P(*ctx->config(), info.x.size(), info.y.size(), info.Lx, info.Ly);
    P.x0 = info.x0;
    P.y0 = info.y0;
    P.registration = grid::CELL_CENTER;
    P.variable_name = variable_name;
    P.vertical_grid_from_options(*ctx->config());
    P.ownership_ranges_from_options(*ctx->config(), ctx->size());

    auto source_grid = std::make_shared<Grid>(ctx, P);

    source_grid->set_mapping_info(mapping);

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

      log->message(2, "Defining the source grid (%s)...\n", source_grid_name.c_str());
      m_source_field_id = define_field(comp_ids[0], *source_grid, source_grid_name);

      log->message(2, "Defining the target grid (%s)...\n", target_grid_name.c_str());
      m_target_field_id = define_field(comp_ids[1], target_grid, target_grid_name);

      // Define the interpolation stack:
      {
        std::string direction;
        int interp_stack_id = 0;
        // FIXME: this will almost always choose conservative interpolation when going
        // from lon,lat to a projected grid because if (...) below compares quantities
        // that have different units.
        if (source_grid->dx() < target_grid.dx() or source_grid->dy() < target_grid.dy()) {
          interp_stack_id = interpolation_fine_to_coarse(fill_value);
          direction       = "fine to coarse (conservative)";
        } else {
          interp_stack_id = interpolation_coarse_to_fine(fill_value);
          direction       = "coarse to fine (distance-weighted sum of neighbors)";
        }
        log->message(2, "Interpolation direction: %s\n", direction.c_str());

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
                                          const Grid &/* target_grid (unused) */,
                                          petsc::Vec &output) const {

  // set metadata to help the following call find the variable, convert units, etc
  m_buffer->metadata(0) = metadata;
  m_buffer->read(file, record_index);

  double time_spent = interpolate(*m_buffer, output);

  return time_spent;
}

} // namespace pism

/* Copyright (C) 2024, 2025, 2026 PISM Authors
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
#include <mpi.h>
#include <vector>
#include <cmath>

#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/InputInterpolation.hh"
#include "pism/util/InputInterpolationYAC.hh"
#include "pism/util/Logger.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/pism_utilities.hh" // GlobalMin(), clip()
#include "pism/util/projection.hh"
#include "pism/util/yac_utilities.hh"

#if (Pism_USE_PROJ == 0)
#error "This code requires PROJ"
#endif

#include <proj.h>

#include "pism/util/Proj.hh"

#if (Pism_USE_YAC == 0)
#error "This code requires YAC"
#endif

extern "C" {
#include "yac.h"
}

namespace pism {


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

  int point_id = pism::define_yac_grid(x, y, name, proj_string);

  const char *time_step_length = "1";
  const int point_set_size     = 1;
  const int collection_size    = 1;

  int field_id = 0;
  yac_cdef_field(name.c_str(), component_id, &point_id, point_set_size, collection_size,
                 time_step_length, YAC_TIME_UNIT_SECOND, &field_id);
  return field_id;
}

static double dx_estimate(Proj &mapping, double x1, double x2, double y) {
  PJ_COORD p1 = proj_coord(0, 0, 0, 0), p2 = proj_coord(0, 0, 0, 0);
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
  PJ_COORD p1 = proj_coord(0, 0, 0, 0), p2 = proj_coord(0, 0, 0, 0);
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
    : m_instance_id(0), m_source_field_id(0), m_target_field_id(0), m_split_comm(MPI_COMM_NULL) {
  auto ctx = target_grid.ctx();

  std::string target_proj_params = target_grid.get_mapping_info()["proj_params"];

  if (target_proj_params.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "internal grid projection is not known");
  }

  try {
    // Initialize a YAC instance:
    {
      yac_cinit_comm_instance(target_grid.com, &m_instance_id);
      yac_cdef_calendar(YAC_YEAR_OF_365_DAYS);
      // Note: zero-padding of months and days *is* required.
      yac_cdef_datetime_instance(m_instance_id, "-1-01-01", "+1-01-01");
    }

    auto log = ctx->log();

    // Note: `input_file` is created on the communicator corresponding to target_grid, so
    // all ranks of target_grid.com have to call functions that use `input_file`:
    auto source_grid_name = grid_name(input_file, variable_name, ctx->unit_system(),
                                      type == PIECEWISE_CONSTANT);
    auto target_grid_name = "internal for " + source_grid_name;
    double target_grid_spacing = std::min(target_grid.dx(), target_grid.dy());

    log->message(
        2, "* Initializing 2D interpolation on the sphere from '%s' to the internal grid...\n",
        source_grid_name.c_str());

    log->message(2, " Internal grid spacing: %3.3f m\n", target_grid_spacing);

    grid::InputGridInfo source_grid_info(input_file, variable_name, ctx->unit_system(),
                                         pism::grid::CELL_CENTER);

    auto source_grid_mapping = mapping_info_from_file(input_file, variable_name, ctx->unit_system());

    std::string grid_mapping_name = source_grid_mapping["grid_mapping_name"];

    std::string source_proj_params = source_grid_mapping["proj_params"];

    if (source_proj_params.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "unsupported or missing projection info for the grid '%s'",
                                    source_grid_name.c_str());
    }

    // Number of source grid rows per MPI process.
    //
    // FIXME: come up with a better way of choosing this number.
    const int rows_per_proc = 100;

    // Many data-sets are stored in the (y,x) order, i.e. the number of grid points in the
    // Y direction is the number of rows. This does not really matter, though -- we just
    // need a way to split the source grid into reasonably-sized chunks. We could split it
    // in *both* X and Y directions, but I don't think that would be any better.
    int nrows = (int)source_grid_info.y.size();

    // Number of MPI processes for the source grid (make sure we use at least one, but no
    // more than ctx->size()):
    int Ny = pism::clip(nrows / rows_per_proc, 1, ctx->size());

    bool io_subcomm = ctx->rank() < Ny;
    MPI_Comm_split(ctx->com(), io_subcomm ? 1 : 0, 0, &m_split_comm);

    // define source and target communicators:
    int target_comp_id = 0;
    int source_comp_id = 0;
    if (io_subcomm) {
      const int n_comps = 2;
      const char *comp_names[n_comps] = {"source_component", "target_component"};
      int comp_ids[n_comps] = {0, 0};
      yac_cdef_comps_instance(m_instance_id, comp_names, n_comps, comp_ids);
      source_comp_id = comp_ids[0];
      target_comp_id = comp_ids[1];
    } else {
      yac_cdef_comp_instance(m_instance_id, "target_component", &target_comp_id);
    }

    // define the target field (performed by all ranks in target_grid.com):
    {
      auto x = grid::subset(target_grid.xs(), target_grid.xm(), target_grid.x());
      auto y = grid::subset(target_grid.ys(), target_grid.ym(), target_grid.y());

      m_target_field_id = define_field(target_comp_id, x, y, target_proj_params, target_grid_name);
    }

    // define the source field on the io_subcomm:
    if (io_subcomm) {
      // create a restriction of this context to a sub-communicator:
      auto io_ctx = ctx->restrict_to_subcomm(m_split_comm, "pism_input_reader");

      auto io_log = io_ctx->log();

      io_log->message(2, "Input grid:\n");
      if (not grid_mapping_name.empty()) {
        io_log->message(2, " Grid mapping: %s\n", grid_mapping_name.c_str());
      }

      io_log->message(2, " PROJ string: '%s'\n", source_proj_params.c_str());

      source_grid_info.report(*io_log, 2, io_ctx->unit_system());

      int Mx = (int)source_grid_info.x.size();
      int My = (int)source_grid_info.y.size();
      grid::Parameters P(*io_ctx->config(), Mx, My, source_grid_info.Lx, source_grid_info.Ly);

      {
        P.x0                = source_grid_info.x0;
        P.y0                = source_grid_info.y0;
        P.registration      = grid::CELL_CENTER;
        P.variable_name     = variable_name;
        P.z                 = { 0.0, 1.0 }; // dummy vertical grid (unused)
        P.max_stencil_width = 0;            // the source grid does not need to support stencils

        // All sub-domains of the source grid are Mx grid points wide in the X direction,
        // i.e. all sub-domains are "strips" of width Mx:
        P.procs_x = { P.Mx };

        // compute widths of sub-domains in the Y direction
        P.procs_y = grid::ownership_ranges(P.My, Ny);
      }

      auto source_grid = std::make_shared<Grid>(io_ctx, P);

      source_grid->set_mapping_info(source_grid_mapping);

      auto x = grid::subset(source_grid->xs(), source_grid->xm(), source_grid_info.x);
      auto y = grid::subset(source_grid->ys(), source_grid->ym(), source_grid_info.y);

      double source_grid_spacing = 0;
      {
        double dx = 0.0;
        double dy = 0.0;
        if (source_grid_info.longitude_latitude) {
          dx = dx_min(source_proj_params, x, y);
          dy = dy_min(source_proj_params, x, y);
        } else {
          dx = std::abs(x[1] - x[0]);
          dy = std::abs(y[1] - y[0]);
        }
        source_grid_spacing = GlobalMin(source_grid->com, std::min(dx, dy));
      }

      m_buffer = std::make_shared<pism::array::Scalar>(source_grid, variable_name);

      io_log->message(2, "Defining the input grid (%s)...\n", source_grid_name.c_str());
      io_log->message(2, " Using %d MPI process%s\n", Ny, Ny > 1 ? "es" : "");
      io_log->message(2, " Input grid spacing: ~%3.3f m\n", source_grid_spacing);
      {
        m_source_field_id =
            define_field(source_comp_id, x, y, source_proj_params, source_grid_name);
      }

      // Define the interpolation stack and the "couple":
      {
        std::string method;
        int interp_stack_id = 0;
        yac_cget_interp_stack_config(&interp_stack_id);

        if (type == PIECEWISE_CONSTANT) {
          method = "nearest neighbor";

          // use nearest neighbor interpolation to interpolate integer fields:
          {
            // nearest neighbor
            int n_neighbors            = 1; // only one neighbor
            double scaling             = 1.0;
            double max_search_distance = 0.0; // unlimited
            yac_cadd_interp_stack_config_nnn(interp_stack_id, YAC_NNN_DIST, n_neighbors,
                                             max_search_distance, scaling);
          }
        } else {
          int partial_coverage = 0;
          if (source_grid_spacing < target_grid_spacing) {
            method = "1st order conservative";

            int order                = 1;
            int enforce_conservation = 1;

            yac_cadd_interp_stack_config_conservative(interp_stack_id, order, enforce_conservation,
                                                      partial_coverage, YAC_CONSERV_DESTAREA);
          } else {
            method = "weighted average of source cell nodes";

            // use average over source grid nodes containing a target point as a backup:
            yac_cadd_interp_stack_config_average(interp_stack_id, YAC_AVG_BARY, partial_coverage);
          }

          {
            // nearest neighbor as a fallback
            int n_neighbors            = 1;
            double scaling             = 1.0;
            double max_search_distance = 0.0; // unlimited
            yac_cadd_interp_stack_config_nnn(interp_stack_id, YAC_NNN_DIST, n_neighbors,
                                             max_search_distance, scaling);
          }
        }

        io_log->message(2, "Interpolation method: %s\n", method.c_str());

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
      } // end of the block defining the interpolation stack and the "couple"

    } // end of "if (io_subcomm) {...}"

    // end definitions:
    {
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
  MPI_Comm_free(&m_split_comm);
  yac_cfinalize_instance(m_instance_id);
}

double InputInterpolationYAC::interpolate(const pism::array::Scalar *source,
                                          pism::petsc::Vec &target) const {
  double start  = MPI_Wtime();
  {
    int collection_size = 1;
    if (source != nullptr) {
      petsc::VecArray input_array(source->vec());
      double *send_field_    = input_array.get();
      double **send_field[1] = { &send_field_ };

      int ierror    = 0;
      int send_info = 0;

      yac_cput(m_source_field_id, collection_size, send_field, &send_info, &ierror);
    }

    pism::petsc::VecArray output_array(target);
    double *recv_field[1] = { output_array.get() };

    int recv_info = 0;
    int ierror    = 0;
    yac_cget(m_target_field_id, collection_size, recv_field, &recv_info, &ierror);
  }
  double end = MPI_Wtime();

  return end - start;
}


void InputInterpolationYAC::regrid(const pism::File &file, pism::array::Scalar &output) const {

  double time_spent =
      InputInterpolation::regrid(output.metadata(0), file, -1, *output.grid(), output.vec());

  auto log = output.grid()->ctx()->log();

  log->message(2, "Interpolation took %f seconds.\n", time_spent);
}


double InputInterpolationYAC::regrid_impl(const VariableMetadata &metadata,
                                          const pism::File &file, int record_index,
                                          const Grid & /* target_grid (unused) */,
                                          petsc::Vec &output) const {

  if (m_buffer != nullptr) {
    // set metadata to help the following call find the variable, convert units, etc
    m_buffer->metadata(0) = metadata;

    // open the file using the sub-communicator used for reading, instead of the
    // communicator corresponding to the target grid (in the `file` argument):
    pism::File input_file(m_split_comm, file.name(), io::PISM_GUESS, io::PISM_READONLY);
    m_buffer->read(input_file, record_index);
  }

  double time_spent = interpolate(m_buffer.get(), output);

  return time_spent;
}

} // namespace pism

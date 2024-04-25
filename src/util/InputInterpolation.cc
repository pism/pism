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

#include "pism/util/InputInterpolation.hh"
#include "pism/pism_config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/projection.hh"
#include "pism/util/Logger.hh"

#if (Pism_USE_YAC_INTERPOLATION == 1)
#include "InputInterpolationYAC.hh"
#endif

#include <memory>

namespace pism {

double InputInterpolation::regrid(const SpatialVariableMetadata &metadata, const pism::File &file,
                                  int record_index, const Grid &grid, petsc::Vec &output) const {
  if (record_index == -1) {
    auto nrecords =
      file.nrecords(metadata.get_name(), metadata["standard_name"], metadata.unit_system());

    record_index = (int)nrecords - 1;
  }

  return regrid_impl(metadata, file, record_index, grid, output);
}

InputInterpolation::InputInterpolation() {
  // empty
}

InputInterpolation3D::InputInterpolation3D(const Grid &target_grid,
                                           const std::vector<double> &levels,
                                           const File &input_file, const std::string &variable_name,
                                           InterpolationType type) {

  auto log         = target_grid.ctx()->log();
  auto unit_system = target_grid.ctx()->unit_system();

  auto name = grid_name(input_file, variable_name, unit_system);

  log->message(2, "* Initializing bi- or tri-linear interpolation from '%s' to the internal grid...\n",
               name.c_str());

  grid::InputGridInfo input_grid(input_file, variable_name, unit_system,
                                 target_grid.registration());

  input_grid.report(*log, 4, unit_system);

  io::check_input_grid(input_grid, target_grid, levels);

  m_interp_context = std::make_shared<LocalInterpCtx>(input_grid, target_grid, levels, type);
}


double InputInterpolation3D::regrid_impl(const SpatialVariableMetadata &metadata,
                                         const pism::File &file,
                                         int record_index,
                                         const Grid &target_grid,
                                         petsc::Vec &output) const {

  petsc::VecArray output_array(output);

  double start = get_time(target_grid.com);
  {
    int old_t_start                 = m_interp_context->start[T_AXIS];
    m_interp_context->start[T_AXIS] = record_index;

    io::regrid_spatial_variable(metadata, target_grid, *m_interp_context, file,
                                output_array.get());

    m_interp_context->start[T_AXIS] = old_t_start;
  }
  double end = get_time(target_grid.com);

  return end - start;
}


std::shared_ptr<InputInterpolation>
InputInterpolation::create(const Grid &target_grid,
                           const std::vector<double> &levels, const File &input_file,
                           const std::string &variable_name, InterpolationType type) {

  auto source_projection =
      get_projection_info(input_file, variable_name, target_grid.ctx()->unit_system()).proj;

  auto target_projection = target_grid.get_mapping_info().proj;

#if (Pism_USE_YAC_INTERPOLATION == 1)
  if (levels.size() < 2 and (not source_projection.empty()) and (not target_projection.empty())) {
    return std::make_shared<InputInterpolationYAC>(target_grid, input_file, variable_name);
  }
#endif

  return std::make_shared<InputInterpolation3D>(target_grid, levels, input_file, variable_name,
                                                type);
}


} // namespace pism

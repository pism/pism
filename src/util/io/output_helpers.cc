/* Copyright (C) 2025 PISM Authors
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
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_utilities.hh"
#include <vector>

namespace pism {
namespace io {

/*!
 * Define a 2D or 3D NetCDF variable and set attributes. Do nothing if a variable is
 * already present.
 */
void define_spatial_variable(const SpatialVariableMetadata &metadata,
                             const VariableMetadata &cf_mapping, const Config &config,
                             const File &file) {

  // Make a copy of `metadata` so we can modify it:
  auto var  = metadata;
  auto name = var.get_name();

  if (file.variable_exists(name)) {
    return;
  }

  if (config.get_flag("output.use_MKS")) {
    var.output_units(var["units"]);
  }

  // add the "grid_mapping" attribute if the grid has an associated mapping. Variables lat, lon,
  // lat_bnds, and lon_bnds should not have the grid_mapping attribute to support CDO (see issue
  // #384).
  if (cf_mapping.has_attributes() and not member(name, { "lat_bnds", "lon_bnds", "lat", "lon" })) {
    var["grid_mapping"] = cf_mapping.get_name();
  }

  std::vector<std::string> dims;

  if (not var.get_time_independent()) {
    dims = { config.get_string("time.dimension_name") };
  }

  // define dimensions and coordinate variables; assemble the list of dimension names:
  //
  // Note the order: y,x,z
  for (const auto &dimension : { var.y(), var.x(), var.z() }) {

    auto dimension_name = dimension.get_name();

    if (dimension_name.empty()) {
      continue;
    }

    dims.push_back(dimension_name);
    file.define_dimension(dimension_name, dimension.length());
    file.define_variable(dimension_name, PISM_DOUBLE, { dimension_name });
    for (const auto& arg : dimension.all_strings()) {
      file.write_attribute(dimension_name, arg.first, arg.second);
    }
  }

  assert(dims.size() > 1);

  // define the variable itself:
  file.define_variable(var.get_name(), var.get_output_type(), dims);
  for (const auto &arg : var.all_strings()) {
    file.write_attribute(var.get_name(), arg.first, arg.second);
  }
  for (const auto &arg : var.all_doubles()) {
    file.write_attribute(var.get_name(), arg.first, var.get_output_type(), arg.second);
  }
}

/*!
 * Write a 1D or 2D array redundantly stored on all MPI ranks.
 *
 * @param[in] file file to write to
 * @param[in] name variable name
 * @param[in] start starting index of the first dimension variable `name` depends on
 * @param[in] M number of elements along the first dimension
 * @param[in] N number of elements along the second dimension (1 for 1D arrays)
 * @param[in] data array to write
 */
void write_array(const File &file, const std::string &variable_name,
                 const std::vector<unsigned int> &start, const std::vector<unsigned int> &count,
                 const std::vector<double> &input) {
  file.write_variable(variable_name, start, count, input.data());
}

/*!
 * Write a 1D or 2D array redundantly stored on all MPI ranks, converting to output units
 * first.
 *
 * @param[in] file file to write to
 * @param[in] metadata variable metadata
 * @param[in] start starting index of the first dimension variable `name` depends on
 * @param[in] M number of elements along the first dimension
 * @param[in] N number of elements along the second dimension (1 for 1D arrays)
 * @param[in] input array to write
 *
 */
void write_array(const File &file, const VariableMetadata &metadata,
                 const std::vector<unsigned int> &start, const std::vector<unsigned int> &count,
                 const std::vector<double> &input) {
  // create a copy of "data" to change units
  std::vector<double> data = input;

  units::Converter(metadata.unit_system(), metadata["units"], metadata["output_units"])
      .convert_doubles(data.data(), data.size());

  write_array(file, metadata.get_name(), start, count, data);
}

//! Append to the time dimension.
void append_time(const File &file, const std::string &name, double value) {
  write_array(file, name, {file.dimension_length(name)}, {1}, { value });
}

//! 
/*! Write a distributed 2D or 3D to a file.
 *
 * Ensures that all coordinate variables this 2D or 3D variable depends on (except time)
 * are written to the file as well.
 *
 * Converts units if internal and "output" units are different.
 *
 * Avoids writing time-independent variables more than once.
 */
void write_spatial_variable(const SpatialVariableMetadata &metadata,
                            const grid::DistributedGridInfo &grid, const Config &config,
                            const File &file, const double *input) {
  auto var = metadata;
  const auto &name = var.get_name();

  if (not file.variable_exists(name)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' in '%s'.", name.c_str(),
                                  file.name().c_str());
  }

  // check if we need to write this variable
  bool time_independent = var.get_time_independent();
  // avoid writing time-independent variables more than once (saves time when writing to
  // extra_files)
  if (time_independent and file.get_variable_was_written(name)) {
    return;
  }

  // write dimensions:
  {
    std::map<std::string, const std::vector<double>&> data = {{var.x().get_name(), grid.x},
                                                              {var.y().get_name(), grid.y},
                                                              {var.z().get_name(), var.levels()}};
    for (const auto &p : data) {
      const auto &dimension = p.first;
      const auto &coordinates = p.second;
      bool exists = file.dimension_exists(dimension);
      bool written = file.get_variable_was_written(dimension);
      if (exists and not written) {
        write_array(file, dimension, { 0 }, { (unsigned int)coordinates.size() }, coordinates);
        file.set_variable_was_written(dimension);
      }
    }
  }

  // override `output_units` if "output.use_MKS" is set
  if (config.get_flag("output.use_MKS")) {
    var.output_units(var["units"]);
  }

  // make sure we have at least one level
  unsigned int nlevels = std::max(var.levels().size(), (size_t)1);

  std::string units = var["units"], output_units = var["output_units"];

  if (units != output_units) {
    size_t data_size = grid.xm * grid.ym * nlevels;

    // create a temporary array, convert to output units, and
    // save
    std::vector<double> tmp(data_size);
    for (size_t k = 0; k < data_size; ++k) {
      tmp[k] = input[k];
    }

    units::Converter(var.unit_system(), units, output_units)
        .convert_doubles(tmp.data(), tmp.size());

    file.write_distributed_array(name, grid, nlevels, not time_independent, tmp.data());
  } else {
    file.write_distributed_array(name, grid, nlevels, not time_independent, input);
  }
  file.set_variable_was_written(name);
}

//! \brief Moves the file aside (file.nc -> file.nc~).
/*!
 * Note: only one processor does the renaming.
 */
void move_if_exists(MPI_Comm com, const std::string &file_to_move, int rank_to_use) {
  int stat = 0, rank = 0;
  MPI_Comm_rank(com, &rank);
  std::string backup_filename = file_to_move + "~";

  if (rank == rank_to_use) {
    bool exists = false;

    // Check if the file exists:
    if (FILE *f = fopen(file_to_move.c_str(), "r")) {
      fclose(f);
      exists = true;
    } else {
      exists = false;
    }

    if (exists) {
      stat = rename(file_to_move.c_str(), backup_filename.c_str());
    }
  } // end of "if (rank == rank_to_use)"

  int global_stat = 0;
  MPI_Allreduce(&stat, &global_stat, 1, MPI_INT, MPI_SUM, com);

  if (global_stat != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "PISM ERROR: can't move '%s' to '%s'",
                                  file_to_move.c_str(), backup_filename.c_str());
  }
}

//! \brief Check if a file is present are remove it.
/*!
 * Note: only processor 0 does the job.
 */
void remove_if_exists(MPI_Comm com, const std::string &file_to_remove, int rank_to_use) {
  int stat = 0, rank = 0;
  MPI_Comm_rank(com, &rank);

  if (rank == rank_to_use) {
    bool exists = false;

    // Check if the file exists:
    if (FILE *f = fopen(file_to_remove.c_str(), "r")) {
      fclose(f);
      exists = true;
    } else {
      exists = false;
    }

    if (exists) {
      stat = remove(file_to_remove.c_str());
    }
  } // end of "if (rank == rank_to_use)"

  int global_stat = 0;
  MPI_Allreduce(&stat, &global_stat, 1, MPI_INT, MPI_SUM, com);

  if (global_stat != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "PISM ERROR: can't remove '%s'",
                                  file_to_remove.c_str());
  }
}

} // namespace io
} // namespace pism

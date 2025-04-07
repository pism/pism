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
#include "pism/util/io/File.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Context.hh"
#include "pism/util/Config.hh"
#include "pism/util/Time.hh"
#include "pism/util/Grid.hh"
#include "pism/util/projection.hh"

namespace pism {
namespace io {

static void define_variable(const File &file, const std::vector<std::string> &dims,
                            io::Type type,
                            const VariableMetadata &metadata) {

  if (file.variable_exists(metadata.get_name())) {
    return;
  }

  file.define_variable(metadata.get_name(), type, dims);

  write_attributes(file, metadata, type);
}

//! \brief Define a dimension \b and the associated coordinate variable. Set attributes.
void define_dimension(const File &file, unsigned long int length,
                      const VariableMetadata &metadata) {
  std::string name = metadata.get_name();
  try {
    if (file.dimension_exists(name)) {
      return;
    }

    file.define_dimension(name, length);

    define_variable(file, { name }, PISM_DOUBLE, metadata);

  } catch (RuntimeError &e) {
    e.add_context("defining dimension '%s' in '%s'", name.c_str(), file.name().c_str());
    throw;
  }
}

//! Prepare a file for output.
void define_time(const File &file, const Context &ctx) {
  const Time &time     = *ctx.time();
  const Config &config = *ctx.config();

  define_time(file, config.get_string("time.dimension_name"), time.calendar(), time.units_string(),
              ctx.unit_system());
}

/*!
 * Define a time dimension and the corresponding coordinate variable. Does nothing if the time
 * variable is already present.
 */
void define_time(const File &file, const std::string &name, const std::string &calendar,
                 const std::string &units, units::System::Ptr unit_system) {
  // time
  VariableMetadata time(name, unit_system);
  time["long_name"] = "time";
  time["calendar"]  = calendar;
  time["units"]     = units;
  time["axis"]      = "T";

  define_dimension(file, PISM_UNLIMITED, time);
}

//! Prepare a file for output.
void append_time(const File &file, const Config &config, double time_seconds) {
  append_time(file, config.get_string("time.dimension_name"), time_seconds);
}

//! \brief Append to the time dimension.
void append_time(const File &file, const std::string &name, double value) {
  try {
    unsigned int start = file.dimension_length(name);

    file.write_variable(name, { start }, { 1 }, &value);
  } catch (RuntimeError &e) {
    e.add_context("appending to the time dimension in \"" + file.name() + "\"");
    throw;
  }
}

//! \brief Define dimensions a variable depends on.
static void define_dimensions(const SpatialVariableMetadata &var,
                              const grid::GridInfo &grid, const File &file) {

  // x
  {
    auto x = var.x();
    x["spacing_meters"] = { grid.x[1] - grid.x[0] };

    define_dimension(file, grid.x.size(), x);
  }

  // y
  {
    auto y = var.y();
    y["spacing_meters"] = { grid.y[1] - grid.y[0] };

    define_dimension(file, grid.y.size(), var.y());
  }

  // z
  {
    auto z = var.z();

    if (not z.get_name().empty()) {
      const auto &levels = var.levels();

      // make sure we have at least one level
      unsigned int nlevels = std::max(levels.size(), (size_t)1);

      bool spatial_dim = not var.z().get_string("axis").empty();

      if (nlevels > 1 and spatial_dim) {
        double dz_max = levels[1] - levels[0];
        double dz_min = levels.back() - levels.front();

        for (unsigned int k = 0; k < nlevels - 1; ++k) {
          double dz = levels[k + 1] - levels[k];
          dz_max    = std::max(dz_max, dz);
          dz_min    = std::min(dz_min, dz);
        }

        z["spacing_min_meters"] = { dz_min };
        z["spacing_max_meters"] = { dz_max };
      }

      define_dimension(file, nlevels, z);
    }
  }
}

static void write_dimension_data(const File &file, const std::string &name,
                                 const std::vector<double> &data) {
  bool exists = file.dimension_exists(name);
  bool written = file.get_variable_was_written(name);
  if (exists and not written) {
    file.write_variable(name, { 0 }, { (unsigned int)data.size() }, data.data());
    file.set_variable_was_written(name);
  }
}

//! Define a NetCDF variable corresponding to a VariableMetadata object.
void define_spatial_variable(const SpatialVariableMetadata &metadata, const Grid &grid,
                             const File &file, io::Type default_type) {
  const auto &config = *grid.ctx()->config();

  // Make a copy of `metadata` so we can modify it:
  auto var = metadata;
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
  const auto &mapping = grid.get_mapping_info().cf_mapping;
  if (mapping.has_attributes() and not member(name, { "lat_bnds", "lon_bnds", "lat", "lon" })) {
    var["grid_mapping"] = mapping.get_name();
  }

  define_dimensions(var, grid.info(), file);

  std::string x = var.x().get_name(), y = var.y().get_name(), z = var.z().get_name();

  std::vector<std::string> dims;

  if (not var.get_time_independent()) {
    auto time_name = config.get_string("time.dimension_name");
    dims = {time_name, y, x};
  } else {
    dims = {y, x};
  }

  if (not z.empty()) {
    dims.push_back(z);
  }

  assert(dims.size() > 1);

  auto type = var.get_output_type();
  if (type == PISM_NAT) {
    type = default_type;
  }

  define_variable(file, dims, type, var);
}

//! \brief Write a double array to a file.
/*!
  Converts units if internal and "output" units are different.
 */
void write_spatial_variable(const SpatialVariableMetadata &metadata,
                            const grid::DistributedGridInfo &grid,
                            const Config &config,
                            const File &file, const double *input) {
  // make a copy of `metadata` so we can override `output_units` if "output.use_MKS" is
  // set.
  SpatialVariableMetadata var = metadata;
  if (config.get_flag("output.use_MKS")) {
    var.output_units(var["units"]);
  }

  // write dimensions:
  {
    // x
    write_dimension_data(file, var.x().get_name(), grid.x);

    // y
    write_dimension_data(file, var.y().get_name(), grid.y);

    // z
    write_dimension_data(file, var.z().get_name(), var.levels());
  }

  bool time_independent = var.get_time_independent();
  bool written = file.get_variable_was_written(var.get_name());

  // avoid writing time-independent variables more than once (saves time when writing to
  // extra_files)
  if (written and time_independent) {
    return;
  }

  const auto &name = var.get_name();

  if (not file.variable_exists(name)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' in '%s'.", name.c_str(),
                                  file.name().c_str());
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
  file.set_variable_was_written(var.get_name());
}

//! Define a NetCDF variable corresponding to a time-series.
void define_timeseries(const VariableMetadata &var, const std::string &dimension_name,
                       const File &file, io::Type output_type) {

  define_dimension(file, PISM_UNLIMITED, { dimension_name, var.unit_system() });

  define_variable(file, { dimension_name }, output_type, var);
}

/** @brief Write a time-series `data` to a file.
 *
 * Always use output units when saving time-series.
 */
void write_timeseries(const File &file, const VariableMetadata &metadata, size_t t_start,
                      const std::vector<double> &data) {

  std::string name = metadata.get_name();
  try {
    // create a copy of "data":
    std::vector<double> tmp = data;

    // convert to output units:
    units::Converter(metadata.unit_system(), metadata["units"], metadata["output_units"])
        .convert_doubles(tmp.data(), tmp.size());

    file.write_variable(name, {(unsigned int)t_start}, {(unsigned int)tmp.size()}, tmp.data());

  } catch (RuntimeError &e) {
    e.add_context("writing time-series variable '%s' to '%s'", name.c_str(),
                  file.name().c_str());
    throw;
  }
}

void define_time_bounds(const VariableMetadata& var,
                        const std::string &time_name,
                        const std::string &bounds_name,
                        const File &file, io::Type output_type) {
  define_dimension(file, 2, { bounds_name, var.unit_system() });
  define_variable(file, { time_name, bounds_name }, output_type, var);
}

void write_time_bounds(const File &file, const VariableMetadata &metadata, size_t t_start,
                       const std::vector<double> &bounds) {

  const auto &name = metadata.get_name();
  try {
    // make a copy of "data"
    auto data = bounds;

    // convert to output units:
    units::Converter(metadata.unit_system(), metadata["units"], metadata["output_units"])
        .convert_doubles(data.data(), data.size());

    file.write_variable(name, { (unsigned int)t_start, 0 }, { (unsigned int)data.size() / 2, 2 },
                        data.data());

  } catch (RuntimeError &e) {
    e.add_context("writing time-bounds variable '%s' to '%s'", name.c_str(), file.name().c_str());
    throw;
  }
}

static void write_attributes(const File &file, const std::string &var_name,
                             const std::map<std::string, std::string> &strings,
                             const std::map<std::string, std::vector<double> > &numbers,
                             io::Type output_type) {
  // Write text attributes:
  for (const auto &s : strings) {
    const auto &name  = s.first;
    const auto &value = s.second;

    if (value.empty()) {
      continue;
    }

    file.write_attribute(var_name, name, value);
  }

  // Write double attributes:
  for (const auto &d : numbers) {
    const auto &name   = d.first;
    const auto &values = d.second;

    if (values.empty()) {
      continue;
    }

    file.write_attribute(var_name, name, output_type, values);
  }
}

//! Write variable attributes to a NetCDF file.
/*!
  - If both valid_min and valid_max are set, then valid_range is written
  instead of the valid_min, valid_max pair.

  - Skips empty text attributes.
*/
void write_attributes(const File &file, const VariableMetadata &metadata, io::Type output_type) {

  // make a copy so we can edit metadata
  auto variable = metadata;
  const auto &var_name = variable.get_name();

  try {
    std::string units = variable["units"], output_units = variable["output_units"];

    // output units should never be written to a file
    variable["output_units"] = "";

    bool use_output_units = (not units.empty() and
                             not output_units.empty() and
                             units != output_units);
    
    if (use_output_units) {
      // Replace "units" with "output_units" and clear "output_units"
      if (variable.has_attribute("units")) {
        variable["units"]        = output_units;
      }

      units::Converter c(variable.unit_system(), units, output_units);

      // We need to convert units of valid_min, valid_max and valid_range:
      {
        if (variable.has_attribute("valid_range")) {
          auto bounds = variable.get_numbers("valid_range");

          variable["valid_range"] = { c(bounds[0]), c(bounds[1]) };
        } else {
          if (variable.has_attribute("valid_min")) {
            auto min = variable.get_number("valid_min");
            variable["valid_min"] = { c(min) };
          }
          if (variable.has_attribute("valid_max")) {
            auto max = variable.get_number("valid_max");
            variable["valid_max"] = { c(max) };
          }
        }

        if (variable.has_attribute("_FillValue")) {
          auto fill = variable.get_number("_FillValue");
          variable["_FillValue"] = { c(fill) };
        }
      }
    }

    write_attributes(file, variable.get_name(), variable.all_strings(), variable.all_doubles(),
                     output_type);

  } catch (RuntimeError &e) {
    e.add_context("writing attributes of variable '%s' to '%s'",
                  var_name.c_str(), file.name().c_str());
    throw;
  }
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

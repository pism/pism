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

//! \brief Define a dimension \b and the associated coordinate variable. Set attributes.
void define_dimension(const File &file, unsigned long int length,
                      const VariableMetadata &metadata) {
  std::string name = metadata.get_name();
  try {
    file.define_dimension(name, length);

    file.define_variable(name, PISM_DOUBLE, { name });

    write_attributes(file, metadata, PISM_DOUBLE);

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
  try {
    if (file.variable_exists(name)) {
      return;
    }

    // time
    VariableMetadata time(name, unit_system);
    time["long_name"] = "time";
    time["calendar"]  = calendar;
    time["units"]     = units;
    time["axis"]      = "T";

    define_dimension(file, PISM_UNLIMITED, time);
  } catch (RuntimeError &e) {
    e.add_context("defining the time dimension in \"" + file.name() + "\"");
    throw;
  }
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
  std::string x_name = var.x().get_name();
  if (not file.dimension_exists(x_name)) {
    define_dimension(file, grid.x.size(), var.x());
    file.write_attribute(x_name, "spacing_meters", PISM_DOUBLE, { grid.x[1] - grid.x[0] });
  }

  // y
  std::string y_name = var.y().get_name();
  if (not file.dimension_exists(y_name)) {
    define_dimension(file, grid.y.size(), var.y());
    file.write_attribute(y_name, "spacing_meters", PISM_DOUBLE, { grid.y[1] - grid.y[0] });
  }

  // z
  std::string z_name = var.z().get_name();
  if (not z_name.empty()) {
    if (not file.dimension_exists(z_name)) {
      const std::vector<double> &levels = var.levels();
      // make sure we have at least one level
      unsigned int nlevels = std::max(levels.size(), (size_t)1);
      define_dimension(file, nlevels, var.z());

      bool spatial_dim = not var.z().get_string("axis").empty();

      if (nlevels > 1 and spatial_dim) {
        double dz_max = levels[1] - levels[0];
        double dz_min = levels.back() - levels.front();

        for (unsigned int k = 0; k < nlevels - 1; ++k) {
          double dz = levels[k + 1] - levels[k];
          dz_max    = std::max(dz_max, dz);
          dz_min    = std::min(dz_min, dz);
        }

        file.write_attribute(z_name, "spacing_min_meters", PISM_DOUBLE, { dz_min });
        file.write_attribute(z_name, "spacing_max_meters", PISM_DOUBLE, { dz_max });
      }
    }
  }
}

static void write_dimension_data(const File &file, const std::string &name,
                                 const std::vector<double> &data) {
  bool written = file.get_variable_was_written(name);
  if (not written) {
    file.write_variable(name, { 0 }, { (unsigned int)data.size() }, data.data());
    file.set_variable_was_written(name);
  }
}

void write_dimensions(const SpatialVariableMetadata &var, const grid::GridInfo &grid,
                      const File &file) {
  // x
  std::string x_name = var.x().get_name();
  if (file.dimension_exists(x_name)) {
    write_dimension_data(file, x_name, grid.x);
  }

  // y
  std::string y_name = var.y().get_name();
  if (file.dimension_exists(y_name)) {
    write_dimension_data(file, y_name, grid.y);
  }

  // z
  std::string z_name = var.z().get_name();
  if (file.dimension_exists(z_name)) {
    write_dimension_data(file, z_name, var.levels());
  }
}

//! Define a NetCDF variable corresponding to a VariableMetadata object.
void define_spatial_variable(const SpatialVariableMetadata &metadata, const Grid &grid,
                             const File &file, io::Type default_type) {
  auto config = grid.ctx()->config();

  // make a copy of `metadata` so we can override `output_units` if "output.use_MKS" is
  // set.
  SpatialVariableMetadata var = metadata;
  if (config->get_flag("output.use_MKS")) {
    var.output_units(var["units"]);
  }

  std::vector<std::string> dims;
  std::string name = var.get_name();

  if (file.variable_exists(name)) {
    return;
  }

  define_dimensions(var, grid.info(), file);

  std::string x = var.x().get_name(), y = var.y().get_name(), z = var.z().get_name();

  if (not var.get_time_independent()) {
    dims.push_back(config->get_string("time.dimension_name"));
  }

  dims.push_back(y);
  dims.push_back(x);

  if (not z.empty()) {
    dims.push_back(z);
  }

  assert(dims.size() > 1);

  io::Type type = var.get_output_type();
  if (type == PISM_NAT) {
    type = default_type;
  }
  file.define_variable(name, type, dims);

  write_attributes(file, var, type);

  // add the "grid_mapping" attribute if the grid has an associated mapping. Variables lat, lon,
  // lat_bnds, and lon_bnds should not have the grid_mapping attribute to support CDO (see issue
  // #384).
  const VariableMetadata &mapping = grid.get_mapping_info().cf_mapping;
  if (mapping.has_attributes() and not member(name, { "lat_bnds", "lon_bnds", "lat", "lon" })) {
    file.write_attribute(var.get_name(), "grid_mapping", mapping.get_name());
  }
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

  auto name = var.get_name();

  if (not file.variable_exists(name)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' in '%s'.", name.c_str(),
                                  file.name().c_str());
  }

  write_dimensions(var, grid, file);

  bool time_independent = var.get_time_independent();
  bool written = file.get_variable_was_written(var.get_name());

  // avoid writing time-independent variables more than once (saves time when writing to
  // extra_files)
  if (written and time_independent) {
    return;
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
                       const File &file, io::Type nctype) {

  std::string name = var.get_name();

  if (file.variable_exists(name)) {
    return;
  }

  if (not file.dimension_exists(dimension_name)) {
    define_dimension(file, PISM_UNLIMITED, VariableMetadata(dimension_name, var.unit_system()));
  }

  if (not file.variable_exists(name)) {
    file.define_variable(name, nctype, { dimension_name });
  }

  write_attributes(file, var, nctype);
}

/** @brief Write a time-series `data` to a file.
 *
 * Always use output units when saving time-series.
 */
void write_timeseries(const File &file, const VariableMetadata &metadata, size_t t_start,
                      const std::vector<double> &data) {

  std::string name = metadata.get_name();
  try {
    if (not file.variable_exists(name)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "variable '%s' not found", name.c_str());
    }

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
                        const std::string &dimension_name,
                        const std::string &bounds_name,
                        const File &file, io::Type nctype) {
  std::string name = var.get_name();

  if (file.variable_exists(name)) {
    return;
  }

  if (not file.dimension_exists(dimension_name)) {
    file.define_dimension(dimension_name, PISM_UNLIMITED);
  }

  if (not file.dimension_exists(bounds_name)) {
    file.define_dimension(bounds_name, 2);
  }

  file.define_variable(name, nctype, {dimension_name, bounds_name});

  write_attributes(file, var, nctype);
}

void write_time_bounds(const File &file, const VariableMetadata &metadata,
                       size_t t_start, const std::vector<double> &data) {

  VariableMetadata var = metadata;

  std::string name = var.get_name();
  try {
    bool variable_exists = file.variable_exists(name);
    if (not variable_exists) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "variable '%s' not found",
                                    name.c_str());
    }

    // make a copy of "data"
    std::vector<double> tmp = data;

    // convert to output units:
    units::Converter(var.unit_system(), var["units"], var["output_units"])
        .convert_doubles(tmp.data(), tmp.size());

    file.write_variable(name,
                        {(unsigned int)t_start, 0},
                        {(unsigned int)tmp.size() / 2, 2},
                        tmp.data());

  } catch (RuntimeError &e) {
    e.add_context("writing time-bounds variable '%s' to '%s'", name.c_str(),
                  file.name().c_str());
    throw;
  }
}

//! Write variable attributes to a NetCDF file.
/*!
  - If both valid_min and valid_max are set, then valid_range is written
  instead of the valid_min, valid_max pair.

  - Skips empty text attributes.
*/
void write_attributes(const File &file, const VariableMetadata &variable, io::Type nctype) {
  std::string var_name = variable.get_name();

  try {
    std::string
      units               = variable["units"],
      output_units = variable["output_units"];

    bool use_output_units = units != output_units;

    // units, valid_min, valid_max and valid_range need special treatment:
    if (variable.has_attribute("units")) {
      file.write_attribute(var_name, "units", use_output_units ? output_units : units);
    }

    std::vector<double> bounds(2);
    if (variable.has_attribute("valid_range")) {
      bounds = variable.get_numbers("valid_range");
    } else {
      if (variable.has_attribute("valid_min")) {
        bounds[0]  = variable.get_number("valid_min");
      }
      if (variable.has_attribute("valid_max")) {
        bounds[1]  = variable.get_number("valid_max");
      }
    }

    double fill_value = 0.0;
    if (variable.has_attribute("_FillValue")) {
      fill_value = variable.get_number("_FillValue");
    }

    // We need to save valid_min, valid_max and valid_range in the units
    // matching the ones in the output.
    if (use_output_units) {

      units::Converter c(variable.unit_system(), units, output_units);

      bounds[0]  = c(bounds[0]);
      bounds[1]  = c(bounds[1]);
      fill_value = c(fill_value);
    }

    if (variable.has_attribute("_FillValue")) {
      file.write_attribute(var_name, "_FillValue", nctype, {fill_value});
    }

    if (variable.has_attribute("valid_range")) {
      file.write_attribute(var_name, "valid_range", nctype, bounds);
    } else if (variable.has_attribute("valid_min") and
               variable.has_attribute("valid_max")) {
      file.write_attribute(var_name, "valid_range", nctype, bounds);
    } else if (variable.has_attribute("valid_min")) {
      file.write_attribute(var_name, "valid_min",   nctype, {bounds[0]});
    } else if (variable.has_attribute("valid_max")) {
      file.write_attribute(var_name, "valid_max",   nctype, {bounds[1]});
    }

    // Write text attributes:
    for (const auto& s : variable.all_strings()) {
      std::string
        name  = s.first,
        value = s.second;

      if (name == "units" or
          name == "output_units" or
          value.empty()) {
        continue;
      }

      file.write_attribute(var_name, name, value);
    }

    // Write double attributes:
    for (const auto& d : variable.all_doubles()) {
      std::string name  = d.first;
      std::vector<double> values = d.second;

      if (member(name, {"valid_min", "valid_max", "valid_range", "_FillValue"}) or
          values.empty()) {
        continue;
      }

      file.write_attribute(var_name, name, nctype, values);
    }

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

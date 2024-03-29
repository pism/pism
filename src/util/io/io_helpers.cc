/* Copyright (C) 2015--2023 PISM Authors
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

#include <cassert>
#include <cmath> // isfinite
#include <cstddef>
#include <memory>
#include <array>
#include <vector>

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/Time.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/interpolation.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/LocalInterpCtx.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/projection.hh"

namespace pism {
namespace io {

//! \brief Bi-(or tri-)linear interpolation.
/*!
 * This is the interpolation code itself.
 *
 * Note that its inputs are (essentially)
 * - the definition of the input grid
 * - the definition of the output grid
 * - input array (lic->buffer)
 * - output array (double *output_array)
 *
 * The `output_array` is expected to be big enough to contain
 * `grid.xm()*`grid.ym()*length(zlevels_out)` numbers.
 *
 * We should be able to switch to using an external interpolation library
 * fairly easily...
 */
static void regrid(const Grid &grid, const LocalInterpCtx &lic, double const *input_array,
                   double *output_array) {
  // We'll work with the raw storage here so that the array we are filling is
  // indexed the same way as the buffer we are pulling from (input_array)

  const int X = X_AXIS,
            Z = Z_AXIS; // indices, just for clarity

  unsigned int nlevels = lic.z->n_output();

  // array sizes for mapping from logical to "flat" indices
  int x_count = lic.count[X], z_count = lic.count[Z];

  for (auto p = grid.points(); p; p.next()) {
    const int i_global = p.i(), j_global = p.j();

    const int i = i_global - grid.xs(), j = j_global - grid.ys();

    // Indices of neighboring points.
    const int X_m = lic.x->left(i), X_p = lic.x->right(i), Y_m = lic.y->left(j),
              Y_p = lic.y->right(j);

    for (unsigned int k = 0; k < nlevels; k++) {

      double a_mm = 0.0, a_mp = 0.0, a_pm = 0.0, a_pp = 0.0;

      if (nlevels > 1) {
        const int Z_m = lic.z->left(k), Z_p = lic.z->right(k);

        const double alpha_z = lic.z->alpha(k);

        // We pretend that there are always 8 neighbors (4 in the map plane,
        // 2 vertical levels). And compute the indices into the input_array for
        // those neighbors.
        const int mmm = (Y_m * x_count + X_m) * z_count + Z_m,
                  mmp = (Y_m * x_count + X_m) * z_count + Z_p,
                  mpm = (Y_m * x_count + X_p) * z_count + Z_m,
                  mpp = (Y_m * x_count + X_p) * z_count + Z_p,
                  pmm = (Y_p * x_count + X_m) * z_count + Z_m,
                  pmp = (Y_p * x_count + X_m) * z_count + Z_p,
                  ppm = (Y_p * x_count + X_p) * z_count + Z_m,
                  ppp = (Y_p * x_count + X_p) * z_count + Z_p;

        // linear interpolation in the z-direction
        a_mm = input_array[mmm] * (1.0 - alpha_z) + input_array[mmp] * alpha_z;
        a_mp = input_array[mpm] * (1.0 - alpha_z) + input_array[mpp] * alpha_z;
        a_pm = input_array[pmm] * (1.0 - alpha_z) + input_array[pmp] * alpha_z;
        a_pp = input_array[ppm] * (1.0 - alpha_z) + input_array[ppp] * alpha_z;
      } else {
        // we don't need to interpolate vertically for the 2-D case
        a_mm = input_array[Y_m * x_count + X_m];
        a_mp = input_array[Y_m * x_count + X_p];
        a_pm = input_array[Y_p * x_count + X_m];
        a_pp = input_array[Y_p * x_count + X_p];
      }

      // interpolation coefficient in the x direction
      const double x_alpha = lic.x->alpha(i);
      // interpolation coefficient in the y direction
      const double y_alpha = lic.y->alpha(j);

      // interpolate in x direction
      const double a_m = a_mm * (1.0 - x_alpha) + a_mp * x_alpha;
      const double a_p = a_pm * (1.0 - x_alpha) + a_pp * x_alpha;

      int index = (j * grid.xm() + i) * nlevels + k;

      // index into the new array and interpolate in x direction
      output_array[index] = a_m * (1.0 - y_alpha) + a_p * y_alpha;
      // done with the point at (x,y,z)
    }
  }
}

struct StartCountInfo {
  std::vector<unsigned int> start;
  std::vector<unsigned int> count;
  std::vector<unsigned int> imap;
};

static StartCountInfo compute_start_and_count(std::vector<AxisType> &dim_types,
                                              std::array<int, 4> start_in,
                                              std::array<int, 4> count_in) {

  auto x_start = start_in[X_AXIS];
  auto x_count = count_in[X_AXIS];
  auto y_start = start_in[Y_AXIS];
  auto y_count = count_in[Y_AXIS];
  auto z_start = start_in[Z_AXIS];
  auto z_count = count_in[Z_AXIS];

  StartCountInfo result;

  auto ndims = dim_types.size();

  // Resize output vectors:
  result.start.resize(ndims);
  result.count.resize(ndims);
  result.imap.resize(ndims);

  // Assemble start, count and imap:
  for (unsigned int j = 0; j < ndims; j++) {
    AxisType dimtype = dim_types[j];

    switch (dimtype) {
    case T_AXIS:
      result.start[j] = start_in[T_AXIS];
      result.count[j] = count_in[T_AXIS];
      result.imap[j]  = x_count * y_count * z_count;
      break;
    case Y_AXIS:
      result.start[j] = y_start;
      result.count[j] = y_count;
      result.imap[j]  = x_count * z_count;
      break;
    case X_AXIS:
      result.start[j] = x_start;
      result.count[j] = x_count;
      result.imap[j]  = z_count;
      break;
    default:
    case Z_AXIS:
      result.start[j] = z_start;
      result.count[j] = z_count;
      result.imap[j]  = 1;
      break;
    }
  }

  return result;
}

//! \brief Define a dimension \b and the associated coordinate variable. Set attributes.
void define_dimension(const File &file, unsigned long int length,
                      const VariableMetadata &metadata) {
  std::string name = metadata.get_name();
  try {
    file.define_dimension(name, length);

    file.define_variable(name, PISM_DOUBLE, { name });

    write_attributes(file, metadata, PISM_DOUBLE);

  } catch (RuntimeError &e) {
    e.add_context("defining dimension '%s' in '%s'", name.c_str(), file.filename().c_str());
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
    if (file.find_variable(name)) {
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
    e.add_context("defining the time dimension in \"" + file.filename() + "\"");
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

    // PIO's I/O type PnetCDF requires this
    file.sync();
  } catch (RuntimeError &e) {
    e.add_context("appending to the time dimension in \"" + file.filename() + "\"");
    throw;
  }
}

//! \brief Define dimensions a variable depends on.
static void define_dimensions(const SpatialVariableMetadata &var, const Grid &grid,
                              const File &file) {

  // x
  std::string x_name = var.x().get_name();
  if (not file.find_dimension(x_name)) {
    define_dimension(file, grid.Mx(), var.x());
    file.write_attribute(x_name, "spacing_meters", PISM_DOUBLE, { grid.x(1) - grid.x(0) });
    file.write_attribute(x_name, "not_written", PISM_INT, { 1.0 });
  }

  // y
  std::string y_name = var.y().get_name();
  if (not file.find_dimension(y_name)) {
    define_dimension(file, grid.My(), var.y());
    file.write_attribute(y_name, "spacing_meters", PISM_DOUBLE, { grid.y(1) - grid.y(0) });
    file.write_attribute(y_name, "not_written", PISM_INT, { 1.0 });
  }

  // z
  std::string z_name = var.z().get_name();
  if (not z_name.empty()) {
    if (not file.find_dimension(z_name)) {
      const std::vector<double> &levels = var.levels();
      // make sure we have at least one level
      unsigned int nlevels = std::max(levels.size(), (size_t)1);
      define_dimension(file, nlevels, var.z());
      file.write_attribute(z_name, "not_written", PISM_INT, { 1.0 });

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
  bool written = file.attribute_type(name, "not_written") == PISM_NAT;
  if (not written) {
    file.write_variable(name, { 0 }, { (unsigned int)data.size() }, data.data());
    file.redef();
    file.remove_attribute(name, "not_written");
  }
}

void write_dimensions(const SpatialVariableMetadata &var, const Grid &grid, const File &file) {
  // x
  std::string x_name = var.x().get_name();
  if (file.find_dimension(x_name)) {
    write_dimension_data(file, x_name, grid.x());
  }

  // y
  std::string y_name = var.y().get_name();
  if (file.find_dimension(y_name)) {
    write_dimension_data(file, y_name, grid.y());
  }

  // z
  std::string z_name = var.z().get_name();
  if (file.find_dimension(z_name)) {
    write_dimension_data(file, z_name, var.levels());
  }
}

/**
 * Check if the storage order of a variable in the current file
 * matches the memory storage order used by PISM.
 *
 * @param[in] file input file
 * @param var_name name of the variable to check
 * @returns false if storage orders match, true otherwise
 */
static bool use_transposed_io(std::vector<AxisType> dimension_types) {

  std::vector<AxisType> storage, memory = { Y_AXIS, X_AXIS };

  bool first = true;
  for (auto dimtype : dimension_types) {

    if (first && dimtype == T_AXIS) {
      // ignore the time dimension, but only if it is the first
      // dimension in the list
      first = false;
      continue;
    }

    if (dimtype == X_AXIS || dimtype == Y_AXIS) {
      storage.push_back(dimtype);
    } else if (dimtype == Z_AXIS) {
      memory.push_back(dimtype); // now memory = {Y_AXIS, X_AXIS, Z_AXIS}
      // assume that this variable has only one Z_AXIS in the file
      storage.push_back(dimtype);
    } else {
      // an UNKNOWN_AXIS or T_AXIS at index != 0 was found, use transposed I/O
      return true;
    }
  }

  // we support 2D and 3D in-memory arrays, but not 4D
  assert(memory.size() <= 3);

  return storage != memory;
}

static std::vector<AxisType> dimension_types(const File &file, const std::string &var_name,
                                             std::shared_ptr<units::System> unit_system) {
  std::vector<AxisType> result;
  for (const auto &dimension : file.dimensions(var_name)) {
    result.push_back(file.dimension_type(dimension, unit_system));
  }
  return result;
}

//! \brief Read an array distributed according to the grid.
static void read_distributed_array(const File &file, const Grid &grid, const std::string &var_name,
                                   unsigned int z_count, unsigned int t_start, double *output) {
  try {
    auto dim_types = dimension_types(file, var_name, grid.ctx()->unit_system());

    auto sc = compute_start_and_count(dim_types,
                                      { (int)t_start, grid.xs(), grid.ys(), 0 },
                                      { 1, grid.xm(), grid.ym(), (int)z_count });

    if (use_transposed_io(dim_types)) {
      file.read_variable_transposed(var_name, sc.start, sc.count, sc.imap, output);
    } else {
      file.read_variable(var_name, sc.start, sc.count, output);
    }

  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' from '%s'", var_name.c_str(), file.filename().c_str());
    throw;
  }
}

/** Regrid `variable_name` from a file, possibly replacing missing values with `default_value`.
 *
 * @param file input file
 * @param variable_name variable to regrid
 * @param internal_grid computational grid; used to initialize interpolation
 */
static std::vector<double> read_for_interpolation(const File &file,
                                                  const std::string &variable_name,
                                                  const Grid &internal_grid,
                                                  const LocalInterpCtx &lic) {

  try {
    auto unit_system = internal_grid.ctx()->unit_system();

    std::vector<double> buffer(lic.buffer_size());

    auto dim_types = dimension_types(file, variable_name, internal_grid.ctx()->unit_system());

    auto sc = compute_start_and_count(dim_types, lic.start, lic.count);

    if (use_transposed_io(dim_types)) {
      file.read_variable_transposed(variable_name, sc.start, sc.count, sc.imap, buffer.data());
    } else {
      file.read_variable(variable_name, sc.start, sc.count, buffer.data());
    }

    // Stop with an error message if some values match the _FillValue attribute:
    {
      auto attribute = file.read_double_attribute(variable_name, "_FillValue");
      if (attribute.size() == 1) {
        double fill_value = attribute[0], epsilon = 1e-12;

        for (const auto &value : buffer) {
          if (fabs(value - fill_value) < epsilon) {
            throw RuntimeError::formatted(
                PISM_ERROR_LOCATION, "Some values of '%s' in '%s' match the _FillValue attribute.",
                variable_name.c_str(), file.filename().c_str());
          }
        }
      }
    }

    return buffer;
  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' from '%s'", variable_name.c_str(),
                  file.filename().c_str());
    throw;
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

  if (file.find_variable(name)) {
    return;
  }

  define_dimensions(var, grid, file);

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
  const VariableMetadata &mapping = grid.get_mapping_info().mapping;
  if (mapping.has_attributes() and not member(name, { "lat_bnds", "lon_bnds", "lat", "lon" })) {
    file.write_attribute(var.get_name(), "grid_mapping", mapping.get_name());
  }

  if (var.get_time_independent()) {
    // mark this variable as "not written" so that write_spatial_variable can avoid
    // writing it more than once.
    file.write_attribute(var.get_name(), "not_written", PISM_INT, { 1.0 });
  }
}

//! Read a variable from a file into an array `output`.
/*! This also converts data from input units to internal units if needed.
 */
void read_spatial_variable(const SpatialVariableMetadata &variable, const Grid &grid,
                           const File &file, unsigned int time, double *output) {

  const Logger &log = *grid.ctx()->log();

  // Find the variable:
  auto var = file.find_variable(variable.get_name(), variable["standard_name"]);

  if (not var.exists) {
    throw RuntimeError::formatted(
        PISM_ERROR_LOCATION, "Can't find '%s' (%s) in '%s'.", variable.get_name().c_str(),
        variable.get_string("standard_name").c_str(), file.filename().c_str());
  }

  // Sanity check: the variable in an input file should have the expected
  // number of spatial dimensions.
  {
    // Set of spatial dimensions for this field:
    std::set<int> axes{ X_AXIS, Y_AXIS };
    if (axis_type_from_string(variable.z()["axis"]) == Z_AXIS) {
      axes.insert(Z_AXIS);
    }

    int input_spatial_dim_count = 0; // number of spatial dimensions (input file)
    size_t matching_dim_count   = 0; // number of matching dimensions

    auto input_dims = file.dimensions(var.name);
    for (const auto &d : input_dims) {
      auto dim_type = file.dimension_type(d, variable.unit_system());

      if (dim_type != T_AXIS) {
        ++input_spatial_dim_count;
      }

      if (axes.find(dim_type) != axes.end()) {
        ++matching_dim_count;
      }
    }

    if (axes.size() != matching_dim_count) {

      // Print the error message and stop:
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "found the %dD variable %s (%s) in '%s' while trying to read\n"
          "'%s' ('%s'), which is %d-dimensional.",
          input_spatial_dim_count, var.name.c_str(), join(input_dims, ",").c_str(),
          file.filename().c_str(), variable.get_name().c_str(),
          variable.get_string("long_name").c_str(), static_cast<int>(axes.size()));
    }
  }

  // make sure we have at least one level
  size_t nlevels = std::max(variable.levels().size(), (size_t)1);

  read_distributed_array(file, grid, var.name, nlevels, time, output);

  std::string input_units           = file.read_text_attribute(var.name, "units");
  const std::string &internal_units = variable["units"];

  if (input_units.empty() and not internal_units.empty()) {
    const std::string &long_name = variable["long_name"];
    log.message(2,
                "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                "              Assuming that it is in '%s'.\n",
                variable.get_name().c_str(), long_name.c_str(), internal_units.c_str());
    input_units = internal_units;
  }

  // Convert data:
  size_t size = grid.xm() * grid.ym() * nlevels;

  // stop if some values match the _FillValue attribute
  {
    auto att = file.read_double_attribute(var.name, "_FillValue");
    if (att.size() == 1) {
      double fill_value = att[0];
      for (size_t k = 0; k < size; ++k) {
        if (output[k] == fill_value) {
          throw RuntimeError::formatted(
              PISM_ERROR_LOCATION,
              "Some values of the variable '%s' in '%s' match the _FillValue attribute.",
              var.name.c_str(), file.filename().c_str());
        }
      }
    }
  }

  units::Converter(variable.unit_system(), input_units, internal_units)
      .convert_doubles(output, size);
}

//! \brief Write a double array to a file.
/*!
  Converts units if internal and "output" units are different.
 */
void write_spatial_variable(const SpatialVariableMetadata &metadata, const Grid &grid,
                            const File &file, const double *input) {
  auto config = grid.ctx()->config();

  // make a copy of `metadata` so we can override `output_units` if "output.use_MKS" is
  // set.
  SpatialVariableMetadata var = metadata;
  if (config->get_flag("output.use_MKS")) {
    var.output_units(var["units"]);
  }

  auto name = var.get_name();

  if (not file.find_variable(name)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' in '%s'.", name.c_str(),
                                  file.filename().c_str());
  }

  write_dimensions(var, grid, file);

  bool time_dependent = not var.get_time_independent();

  // avoid writing time-independent variables more than once (saves time when writing to
  // extra_files)
  if (not time_dependent) {
    bool written = file.attribute_type(var.get_name(), "not_written") == PISM_NAT;
    if (written) {
      return;
    }

    file.redef();
    file.remove_attribute(var.get_name(), "not_written");
  }

  // make sure we have at least one level
  unsigned int nlevels = std::max(var.levels().size(), (size_t)1);

  std::string units = var["units"], output_units = var["output_units"];

  if (units != output_units) {
    size_t data_size = grid.xm() * grid.ym() * nlevels;

    // create a temporary array, convert to output units, and
    // save
    std::vector<double> tmp(data_size);
    for (size_t k = 0; k < data_size; ++k) {
      tmp[k] = input[k];
    }

    units::Converter(var.unit_system(), units, output_units)
        .convert_doubles(tmp.data(), tmp.size());

    file.write_distributed_array(name, grid, nlevels, time_dependent, tmp.data());
  } else {
    file.write_distributed_array(name, grid, nlevels, time_dependent, input);
  }
}

/*!
 * Check the overlap of the input grid and the internal grid.
 *
 * Set `allow_extrapolation` to `true` to "extend" the vertical grid during "bootstrapping".
 */
static void check_grid_overlap(const grid::InputGridInfo &input, const Grid &internal,
                               const std::vector<double> &z_internal) {

  // Grid spacing (assume that the grid is equally-spaced) and the
  // extent of the domain. To compute the extent of the domain, assume
  // that the grid is cell-centered, so edge of the domain is half of
  // the grid spacing away from grid points at the edge.
  //
  // Note that x_min is not the same as internal.x(0).
  const double x_min = internal.x0() - internal.Lx(), x_max = internal.x0() + internal.Lx(),
               y_min = internal.y0() - internal.Ly(), y_max = internal.y0() + internal.Ly(),
               input_x_min = input.x0 - input.Lx, input_x_max = input.x0 + input.Lx,
               input_y_min = input.y0 - input.Ly, input_y_max = input.y0 + input.Ly;

  // tolerance (one micron)
  double eps = 1e-6;

  // horizontal grid extent
  if (not(x_min >= input_x_min - eps and x_max <= input_x_max + eps and
          y_min >= input_y_min - eps and y_max <= input_y_max + eps)) {
    throw RuntimeError::formatted(
        PISM_ERROR_LOCATION,
        "PISM's computational domain is not a subset of the domain in '%s'\n"
        "PISM grid:       x: [%3.3f, %3.3f] y: [%3.3f, %3.3f] meters\n"
        "input file grid: x: [%3.3f, %3.3f] y: [%3.3f, %3.3f] meters",
        input.filename.c_str(), x_min, x_max, y_min, y_max, input_x_min, input_x_max, input_y_min,
        input_y_max);
  }

  if (z_internal.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Internal vertical grid has 0 levels. This should never happen.");
  }

  if (z_internal.size() == 1 and input.z.size() > 1) {
    // internal field is 2D or 3D with one level, input variable is 3D with more than one
    // vertical level
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "trying to read in a 2D field but the input file %s contains\n"
                                  "a 3D field with %d levels",
                                  input.filename.c_str(), static_cast<int>(input.z.size()));
  }

  if (z_internal.size() > 1 and input.z.size() <= 1) {
    // internal field is 3D with more than one vertical level, input variable is 2D or 3D
    // with 1 level
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "trying to read in a 3D field but the input file %s contains\n"
                                  "a 2D field",
                                  input.filename.c_str());
  }

  if (z_internal.size() > 1 and (not input.z.empty())) {
    // both internal field and input variable are 3D: check vertical grid extent
    // Note: in PISM 2D fields have one vertical level (z = 0).
    const double input_z_min = input.z.front(), input_z_max = input.z.back(),
                 z_min = z_internal.front(), z_max = z_internal.back();

    if (not(z_min >= input.z_min - eps and z_max <= input.z_max + eps)) {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "PISM's computational domain is not a subset of the domain in '%s'\n"
          "PISM grid:       z: [%3.3f, %3.3f] meters\n"
          "input file grid: z: [%3.3f, %3.3f] meters",
          input.filename.c_str(), z_min, z_max, input_z_min, input_z_max);
    }
  }
}

/*! @brief Check that x, y, and z coordinates of the input grid are strictly increasing. */
void check_input_grid(const grid::InputGridInfo &input_grid,
                      const Grid& internal_grid,
                      const std::vector<double> &internal_z_levels) {
  if (not is_increasing(input_grid.x)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "input x coordinate has to be strictly increasing");
  }

  if (not is_increasing(input_grid.y)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "input y coordinate has to be strictly increasing");
  }

  if (not is_increasing(input_grid.z)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "input vertical grid has to be strictly increasing");
  }

  bool allow_extrapolation = internal_grid.ctx()->config()->get_flag("grid.allow_extrapolation");

  if (not allow_extrapolation) {
    check_grid_overlap(input_grid, internal_grid, internal_z_levels);
  }
}

//! \brief Regrid from a NetCDF file into a distributed array `output`.
/*!
  - if `flag` is `CRITICAL` or `CRITICAL_FILL_MISSING`, stops if the
  variable was not found in the input file
  - if `flag` is one of `CRITICAL_FILL_MISSING` and
  `OPTIONAL_FILL_MISSING`, replace _FillValue with `default_value`.
  - sets `v` to `default_value` if `flag` is `OPTIONAL` and the
  variable was not found in the input file
  - uses the last record in the file
*/

void regrid_spatial_variable(SpatialVariableMetadata &variable,
                             const Grid &internal_grid,
                             const LocalInterpCtx &lic, const File &file,
                             double *output) {

  auto var_info = file.find_variable(variable.get_name(), variable["standard_name"]);
  auto variable_name = var_info.name;

  const Profiling &profiling = internal_grid.ctx()->profiling();

  profiling.begin("io.regridding.read");
  auto buffer = read_for_interpolation(file, variable_name, internal_grid, lic);
  profiling.end("io.regridding.read");

  // interpolate
  profiling.begin("io.regridding.interpolate");
  regrid(internal_grid, lic, buffer.data(), output);
  profiling.end("io.regridding.interpolate");

  // Get the units string from the file and convert the units:
  {
    std::string input_units    = file.read_text_attribute(variable_name, "units");
    std::string internal_units = variable["units"];

    if (input_units.empty() and not internal_units.empty()) {
      const Logger &log = *internal_grid.ctx()->log();
      log.message(2,
                  "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                  "              Assuming that it is in '%s'.\n",
                  variable.get_name().c_str(), variable.get_string("long_name").c_str(),
                  internal_units.c_str());
      input_units = internal_units;
    }

    const size_t data_size = internal_grid.xm() * internal_grid.ym() * lic.z->n_output();

    // Convert data:
    units::Converter(variable.unit_system(), input_units, internal_units)
        .convert_doubles(output, data_size);
  }

  read_valid_range(file, variable_name, variable);
}


//! Define a NetCDF variable corresponding to a time-series.
void define_timeseries(const VariableMetadata &var, const std::string &dimension_name,
                       const File &file, io::Type nctype) {

  std::string name = var.get_name();

  if (file.find_variable(name)) {
    return;
  }

  if (not file.find_dimension(dimension_name)) {
    define_dimension(file, PISM_UNLIMITED, VariableMetadata(dimension_name, var.unit_system()));
  }

  if (not file.find_variable(name)) {
    file.define_variable(name, nctype, { dimension_name });
  }

  write_attributes(file, var, nctype);
}

//! Read a time-series variable from a NetCDF file to a vector of doubles.
void read_timeseries(const File &file, const VariableMetadata &metadata, const Logger &log,
                     std::vector<double> &data) {

  std::string name = metadata.get_name();

  try {
    // Find the variable:
    std::string long_name = metadata["long_name"], standard_name = metadata["standard_name"];

    auto var = file.find_variable(name, standard_name);

    if (not var.exists) {
      throw RuntimeError(PISM_ERROR_LOCATION, "variable " + name + " is missing");
    }

    auto dims = file.dimensions(var.name);
    if (dims.size() != 1) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "variable '%s' in '%s' should to have 1 dimension (got %d)",
                                    name.c_str(), file.filename().c_str(), (int)dims.size());
    }

    auto dimension_name = dims[0];

    unsigned int length = file.dimension_length(dimension_name);
    if (length <= 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "dimension " + dimension_name + " has length zero");
    }

    data.resize(length); // memory allocation happens here

    file.read_variable(var.name, { 0 }, { length }, data.data());

    units::System::Ptr system = metadata.unit_system();
    units::Unit internal_units(system, metadata["units"]), input_units(system, "1");

    std::string input_units_string = file.read_text_attribute(var.name, "units");

    bool input_has_units = not input_units_string.empty();

    if (input_has_units) {
      input_units = units::Unit(system, input_units_string);
    }

    if (metadata.has_attribute("units") && not input_has_units) {
      std::string units_string = internal_units.format();
      log.message(2,
                  "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                  "              Assuming that it is in '%s'.\n",
                  name.c_str(), long_name.c_str(), units_string.c_str());
      input_units = internal_units;
    }

    units::Converter(input_units, internal_units).convert_doubles(data.data(), data.size());

  } catch (RuntimeError &e) {
    e.add_context("reading time-series variable '%s' from '%s'", name.c_str(),
                  file.filename().c_str());
    throw;
  }
}

/** @brief Write a time-series `data` to a file.
 *
 * Always use output units when saving time-series.
 */
void write_timeseries(const File &file, const VariableMetadata &metadata, size_t t_start,
                      const std::vector<double> &data) {

  std::string name = metadata.get_name();
  try {
    if (not file.find_variable(name)) {
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
                  file.filename().c_str());
    throw;
  }
}

void define_time_bounds(const VariableMetadata& var,
                        const std::string &dimension_name,
                        const std::string &bounds_name,
                        const File &file, io::Type nctype) {
  std::string name = var.get_name();

  if (file.find_variable(name)) {
    return;
  }

  if (not file.find_dimension(dimension_name)) {
    file.define_dimension(dimension_name, PISM_UNLIMITED);
  }

  if (not file.find_dimension(bounds_name)) {
    file.define_dimension(bounds_name, 2);
  }

  file.define_variable(name, nctype, {dimension_name, bounds_name});

  write_attributes(file, var, nctype);
}

void read_time_bounds(const File &file,
                      const VariableMetadata &metadata,
                      const Logger &log,
                      std::vector<double> &data) {

  std::string name = metadata.get_name();

  try {
    auto system = metadata.unit_system();
    units::Unit internal_units(system, metadata["units"]);

    // Find the variable:
    if (not file.find_variable(name)) {
      throw RuntimeError(PISM_ERROR_LOCATION, "variable " + name + " is missing");
    }

    std::vector<std::string> dims = file.dimensions(name);

    if (dims.size() != 2) {
      throw RuntimeError(PISM_ERROR_LOCATION, "variable " + name + " has to has two dimensions");
    }

    std::string
      &dimension_name = dims[0],
      &bounds_name    = dims[1];

    // Check that we have 2 vertices (interval end-points) per time record.
    unsigned int length = file.dimension_length(bounds_name);
    if (length != 2) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "time-bounds variable " + name + " has to have exactly 2 bounds per time record");
    }

    // Get the number of time records.
    length = file.dimension_length(dimension_name);
    if (length <= 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "dimension " + dimension_name + " has length zero");
    }

    data.resize(2*length);                // memory allocation happens here

    file.read_variable(name, {0, 0}, {length, 2}, data.data());

    // Find the corresponding 'time' variable. (We get units from the 'time'
    // variable, because according to CF-1.5 section 7.1 a "boundary variable"
    // may not have metadata set.)
    if (not file.find_variable(dimension_name)) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "time coordinate variable " + dimension_name + " is missing");
    }

    bool input_has_units = false;
    units::Unit input_units(internal_units.system(), "1");

    std::string input_units_string = file.read_text_attribute(dimension_name, "units");
    if (input_units_string.empty()) {
      input_has_units = false;
    } else {
      input_units = units::Unit(internal_units.system(), input_units_string);
      input_has_units = true;
    }

    if (metadata.has_attribute("units") && not input_has_units) {
      std::string units_string = internal_units.format();
      log.message(2,
                  "PISM WARNING: Variable '%s' does not have the units attribute.\n"
                  "              Assuming that it is in '%s'.\n",
                  dimension_name.c_str(),
                  units_string.c_str());
      input_units = internal_units;
    }

    units::Converter(input_units, internal_units).convert_doubles(data.data(), data.size());

    // FIXME: check that time intervals described by the time bounds
    // variable are contiguous (without gaps) and stop if they are not.
  } catch (RuntimeError &e) {
    e.add_context("reading time bounds variable '%s' from '%s'", name.c_str(),
                  file.filename().c_str());
    throw;
  }
}

void write_time_bounds(const File &file, const VariableMetadata &metadata,
                       size_t t_start, const std::vector<double> &data) {

  VariableMetadata var = metadata;

  std::string name = var.get_name();
  try {
    bool variable_exists = file.find_variable(name);
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
                  file.filename().c_str());
    throw;
  }
}

/*!
 * Reads and validates times and time bounds.
 */
void read_time_info(const Logger &log,
                    std::shared_ptr<units::System> unit_system,
                    const File &file,
                    const std::string &time_name,
                    const std::string &time_units,
                    std::vector<double> &times,
                    std::vector<double> &bounds) {

  size_t N = 0;
  {
    VariableMetadata time_variable(time_name, unit_system);
    time_variable["units"] = time_units;

    io::read_timeseries(file, time_variable, log, times);

    if (not is_increasing(times)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "times have to be strictly increasing");
    }
    N = times.size();
  }

  // Read time bounds
  {
    std::string time_bounds_name = file.read_text_attribute(time_name, "bounds");

    if (time_bounds_name.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "please provide time bounds for '%s'",
                                    time_name.c_str());
    }

    VariableMetadata bounds_variable(time_bounds_name, unit_system);
    bounds_variable["units"] = time_units;

    io::read_time_bounds(file, bounds_variable, log, bounds);

    if (2 * N != bounds.size()) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "each time record has to have 2 bounds");
    }

    for (size_t k = 0; k < N; ++k) {
      if (not (times[k] >= bounds[2 * k + 0] and
               times[k] <= bounds[2 * k + 1])) {
        throw RuntimeError(PISM_ERROR_LOCATION,
                           "each time has to be contained in its time interval");
      }
    }
  } // end of the block reading time bounds
}

bool file_exists(MPI_Comm com, const std::string &filename) {
  int file_exists_flag = 0, rank = 0;
  MPI_Comm_rank(com, &rank);

  if (rank == 0) {
    // Check if the file exists:
    if (FILE *f = fopen(filename.c_str(), "r")) {
      file_exists_flag = 1;
      fclose(f);
    } else {
      file_exists_flag = 0;
    }
  }
  MPI_Bcast(&file_exists_flag, 1, MPI_INT, 0, com);

  return file_exists_flag == 1;
}

void read_attributes(const File &file,
                     const std::string &variable_name,
                     VariableMetadata &variable) {
  try {
    if (not file.find_variable(variable_name)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "variable '%s' is missing", variable_name.c_str());
    }

    variable.clear();

    unsigned int nattrs = file.nattributes(variable_name);

    for (unsigned int j = 0; j < nattrs; ++j) {
      std::string attribute_name = file.attribute_name(variable_name, j);
      io::Type nctype = file.attribute_type(variable_name, attribute_name);

      if (nctype == PISM_CHAR) {
        variable[attribute_name] = file.read_text_attribute(variable_name,
                                                            attribute_name);
      } else {
        variable[attribute_name] = file.read_double_attribute(variable_name,
                                                              attribute_name);
      }
    } // end of for (int j = 0; j < nattrs; ++j)
  } catch (RuntimeError &e) {
    e.add_context("reading attributes of variable '%s' from '%s'",
                  variable_name.c_str(), file.filename().c_str());
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
                  var_name.c_str(), file.filename().c_str());
    throw;
  }
}

//! Read the valid range information from a file.
/*! Reads `valid_min`, `valid_max` and `valid_range` attributes; if \c
  valid_range is found, sets the pair `valid_min` and `valid_max` instead.
*/
void read_valid_range(const File &file, const std::string &name, VariableMetadata &variable) {
  try {
    // Never reset valid_min/max if they were set internally
    if (variable.has_attribute("valid_min") or
        variable.has_attribute("valid_max")) {
      return;
    }

    // Read the units.
    std::string file_units = file.read_text_attribute(name, "units");

    if (file_units.empty()) {
      // If the variable in the file does not have the units attribute we assume that
      // units in the file match internal (PISM) units.
      file_units = variable.get_string("units");
    }

    units::Converter c(variable.unit_system(), file_units, variable["units"]);

    std::vector<double> bounds = file.read_double_attribute(name, "valid_range");
    if (bounds.size() == 2) {             // valid_range is present
      variable["valid_min"] = {c(bounds[0])};
      variable["valid_max"] = {c(bounds[1])};
    } else {                      // valid_range has the wrong length or is missing
      bounds = file.read_double_attribute(name, "valid_min");
      if (bounds.size() == 1) {           // valid_min is present
        variable["valid_min"] = {c(bounds[0])};
      }

      bounds = file.read_double_attribute(name, "valid_max");
      if (bounds.size() == 1) {           // valid_max is present
        variable["valid_max"] = {c(bounds[0])};
      }
    }
  } catch (RuntimeError &e) {
    e.add_context("reading valid range of variable '%s' from '%s'", name.c_str(),
                  file.filename().c_str());
    throw;
  }
}

/*!
 * Return the name of the time dimension corresponding to a NetCDF variable.
 *
 * Returns an empty string if this variable is time-independent.
 */
std::string time_dimension(units::System::Ptr unit_system,
                           const File &file,
                           const std::string &variable_name) {

  auto dims = file.dimensions(variable_name);

  for (const auto &d : dims) {
    if (file.dimension_type(d, unit_system) == T_AXIS) {
      return d;
    }
  }

  return "";
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
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "PISM ERROR: can't move '%s' to '%s'",
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
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "PISM ERROR: can't remove '%s'", file_to_remove.c_str());
  }
}

} // end of namespace io
} // end of namespace pism

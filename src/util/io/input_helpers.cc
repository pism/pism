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

#include <cmath>

#include "pism/util/Grid.hh"
#include "pism/util/io/LocalInterpCtx.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/Context.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Logger.hh"

namespace pism {
namespace io {

/*!
 * Check if units are set in an input file and warn if they are not.
 */
static std::string check_units(const VariableMetadata &variable, const std::string &input_units,
                               const Logger &log) {
  std::string internal_units = variable["units"];
  if (input_units.empty() and not internal_units.empty()) {
    log.message(2,
                "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                "              Assuming that it is in '%s'.\n",
                variable.get_name().c_str(), variable.get_string("long_name").c_str(),
                internal_units.c_str());
    return internal_units;
  }

  return input_units;
}

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
static void interpolate(const grid::DistributedGridInfo &grid, const LocalInterpCtx &lic, double const *input_array,
                        double *output_array) {
  // We'll work with the raw storage here so that the array we are filling is
  // indexed the same way as the buffer we are pulling from (input_array)

  unsigned int nlevels = lic.z->n_output();

  int x_count = lic.count[X_AXIS], z_count = lic.count[Z_AXIS];
  auto input = [input_array, x_count, z_count](int X, int Y, int Z) {
    // the map from logical to linear indices for the input array
    int index = (Y * x_count + X) * z_count + Z;
    return input_array[index];
  };

  for (auto p : GridPoints(grid)) {
    const int i_global = p.i(), j_global = p.j();

    const int i = i_global - grid.xs, j = j_global - grid.ys;

    // Indices of neighboring points.
    const int X_m = lic.x->left(i), X_p = lic.x->right(i);
    const int Y_m = lic.y->left(j), Y_p = lic.y->right(j);

    for (unsigned int k = 0; k < nlevels; k++) {

      double a_mm = 0.0, a_mp = 0.0, a_pm = 0.0, a_pp = 0.0;

      if (nlevels > 1) {
        const int Z_m = lic.z->left(k), Z_p = lic.z->right(k);

        const double alpha_z = lic.z->alpha(k);

        // We pretend that there are always 8 neighbors (4 in the map plane,
        // 2 vertical levels). And compute the indices into the input_array for
        // those neighbors.
        const double
          mmm = input(X_m, Y_m, Z_m),
          mmp = input(X_m, Y_m, Z_p),
          pmm = input(X_p, Y_m, Z_m),
          pmp = input(X_p, Y_m, Z_p),
          mpm = input(X_m, Y_p, Z_m),
          mpp = input(X_m, Y_p, Z_p),
          ppm = input(X_p, Y_p, Z_m),
          ppp = input(X_p, Y_p, Z_p);

        // linear interpolation in the z-direction
        a_mm = mmm * (1.0 - alpha_z) + mmp * alpha_z;
        a_mp = pmm * (1.0 - alpha_z) + pmp * alpha_z;
        a_pm = mpm * (1.0 - alpha_z) + mpp * alpha_z;
        a_pp = ppm * (1.0 - alpha_z) + ppp * alpha_z;
      } else {
        // no interpolation in Z in the 2-D case
        a_mm = input(X_m, Y_m, 0);
        a_mp = input(X_p, Y_m, 0);
        a_pm = input(X_m, Y_p, 0);
        a_pp = input(X_p, Y_p, 0);
      }

      // interpolation coefficient in the x direction
      const double x_alpha = lic.x->alpha(i);
      // interpolation coefficient in the y direction
      const double y_alpha = lic.y->alpha(j);

      // interpolate in x direction
      const double a_m = a_mm * (1.0 - x_alpha) + a_mp * x_alpha;
      const double a_p = a_pm * (1.0 - x_alpha) + a_pp * x_alpha;

      auto index = (j * grid.xm + i) * nlevels + k;

      // index into the new array and interpolate in y direction
      output_array[index] = a_m * (1.0 - y_alpha) + a_p * y_alpha;
      // done with the point at (x,y,z)
    } // end of the loop over vertical levels
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
static bool transpose(std::vector<AxisType> dimension_types) {

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

/*!
 * Transpose data in `input`, putting results in `output`.
 *
 * We assume that both `input` and `output` hold prod(`count`) elements.
 *
 * The `output` array uses the Y,X,Z order (columns in Z are contiguous).
 *
 * The `input` array uses the ordering corresponding to `input_axes` (ordering present in
 * an input file).
 *
 * The array `count` provides the size of a block in `input`, listing axes in the order of
 * values of AxisType (T,X,Y,Z).
 */
static void transpose(const double *input, const std::vector<AxisType> &input_axes,
                      const std::array<int, 4> &count, double *output) {
  // delta[X_AXIS] is the change in the linear index corresponding to incrementing x in
  // the `input` ordering. delta[Y_AXIS], delta[Z_AXIS] and delta[T_AXIS] correspond to
  // changes in y, z, t.
  std::vector<unsigned> delta = {1, 1, 1, 1}; // 4 to store steps for T,Y,X,Z axes
  {
    int N = (int)input_axes.size();
    // compute changes in the linear index corresponding to incrementing one of the
    // "spatial" indexes, in the order used in `input`:
    std::vector<unsigned> tmp(N, 1);
    for (int k = 0; k < N; ++k) {
      for (int n = k + 1; n < N; ++n) {
        tmp[k] *= count[input_axes[n]];
      }
    }
    // re-arrange so that they are stored in the `T,X,Y,Z` order:
    for (int k = 0; k < N; ++k) {
      delta[input_axes[k]] = tmp[k];
    }
  }

  // change in the linear index corresponding to incrementing x (`output` ordering):
  unsigned delta_x = count[Z_AXIS];
  // change in the linear index corresponding to incrementing y (`output` ordering):
  unsigned delta_y = count[X_AXIS] * count[Z_AXIS];

  // traverse in the memory storage order:
  for (int y = 0; y < count[Y_AXIS]; ++y) {
    for (int x = 0; x < count[X_AXIS]; ++x) {
      for (int z = 0; z < count[Z_AXIS]; ++z) {
        auto OUT = x * delta_x + y * delta_y + z * 1;
        auto IN =  x * delta[X_AXIS] + y * delta[Y_AXIS] + z * delta[Z_AXIS];

        output[OUT] = input[IN];
      }
    }
  }
}

/*!
 * Check if some values in `buffer` match the _FillValue attribute and stop with an error
 * message if such values are found.
 */
static void check_for_missing_values(const File &file, const std::string &variable_name,
                                     double tolerance, const double *buffer, size_t buffer_length) {
  auto attribute = file.read_double_attribute(variable_name, "_FillValue");
  if (attribute.size() == 1) {
    double fill_value = attribute[0];

    for (size_t k = 0; k < buffer_length; ++k) {
      if (fabs(buffer[k] - fill_value) < tolerance) {
        throw RuntimeError::formatted(
            PISM_ERROR_LOCATION,
            "Variable '%s' in '%s' contains values matching the _FillValue attribute",
            variable_name.c_str(), file.name().c_str());
      }
      if (not std::isfinite(buffer[k])) {
        throw RuntimeError::formatted(
            PISM_ERROR_LOCATION,
            "Variable '%s' in '%s' contains values that are not finite (NaN or infinity)",
            variable_name.c_str(), file.name().c_str());
      }
    }
  }
}

struct StartCountInfo {
  std::vector<unsigned int> start;
  std::vector<unsigned int> count;
};


/*!
 * Assemble start and count arrays for use with I/O function calls.
 *
 * This function re-arranges provided `start` and `count` (listed in the order T,X,Y,Z) to
 * get start and count in the order matching the one in `dim_types`.
 */
static StartCountInfo compute_start_and_count(std::vector<AxisType> &dim_types,
                                              const std::array<int, 4> &start,
                                              const std::array<int, 4> &count) {
  auto ndims = dim_types.size();

  // Resize output vectors:
  StartCountInfo result;
  result.start.resize(ndims);
  result.count.resize(ndims);

  // Assemble start and count:
  for (unsigned int j = 0; j < ndims; j++) {
    switch (dim_types[j]) {
    case T_AXIS:
      result.start[j] = start[T_AXIS];
      result.count[j] = count[T_AXIS];
      break;
    case Y_AXIS:
      result.start[j] = start[Y_AXIS];
      result.count[j] = count[Y_AXIS];
      break;
    case X_AXIS:
      result.start[j] = start[X_AXIS];
      result.count[j] = count[X_AXIS];
      break;
    case EXP_ID_AXIS:
      result.start[j] = 0;
      result.count[j] = 1;
      break;
    default:
      // Note: the "default" case is used to handle "3D" variables where the third axis is
      // not a "Z" axis, i.e. dim_types[j] == UNKNOWN_AXIS. We use this to write
      // deposition times in the isochrone tracking model, latitude and longitude bounds,
      // etc. In all these cases the data for the "unknown" axis is input in the slot for
      // the "Z" axis. (At least this matches the in-memory storage order.)
    case Z_AXIS:
      result.start[j] = start[Z_AXIS];
      result.count[j] = count[Z_AXIS];
      break;
    }
  }

  return result;
}

//! \brief Read an array distributed according to the grid.
static void read_distributed_array(const File &file, const std::string &variable_name,
                                   std::shared_ptr<units::System> unit_system,
                                   const std::array<int,4> &start,
                                   const std::array<int,4> &count,
                                   double *output) {
  try {
    auto dim_types = dimension_types(file, variable_name, unit_system);

    auto sc = compute_start_and_count(dim_types, start, count);

    auto size = count[X_AXIS] * count[Y_AXIS] * count[Z_AXIS];

    if (transpose(dim_types)) {
      std::vector<double> tmp(size);
      file.read_variable(variable_name, sc.start, sc.count, tmp.data());
      transpose(tmp.data(), dim_types, count, output);
    } else {
      file.read_variable(variable_name, sc.start, sc.count, output);
    }

    // Stop with an error message if some values match the _FillValue attribute:
    check_for_missing_values(file, variable_name, 1e-12, output, size);

  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' from '%s'", variable_name.c_str(), file.name().c_str());
    throw;
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
void regrid_spatial_variable(const SpatialVariableMetadata &variable,
                             const Grid &target_grid,
                             const LocalInterpCtx &interp_context, const File &file,
                             const Logger &log,
                             double *output) {

  auto var = file.find_variable(variable.get_name(), variable["standard_name"]);

  const auto &profiling = target_grid.ctx()->profiling();

  profiling.begin("io.regridding.read");
  std::vector<double> buffer(interp_context.buffer_size());
  read_distributed_array(file, var.name, variable.unit_system(), interp_context.start,
                         interp_context.count, buffer.data());
  profiling.end("io.regridding.read");

  // interpolate
  profiling.begin("io.regridding.interpolate");
  interpolate(target_grid.info(), interp_context, buffer.data(), output);
  profiling.end("io.regridding.interpolate");

  // Get the units string from the file and convert the units:
  {
    std::string input_units    = file.read_text_attribute(var.name, "units");
    std::string internal_units = variable["units"];

    input_units = check_units(variable, input_units, log);

    const size_t data_size = target_grid.xm() * target_grid.ym() * interp_context.z->n_output();

    // Convert data:
    units::Converter(variable.unit_system(), input_units, internal_units)
        .convert_doubles(output, data_size);
  }
}

std::vector<double> read_1d_variable(const File &file, const std::string &variable_name,
                                     const std::string &units,
                                     std::shared_ptr<units::System> unit_system) {


  try {
    if (not file.variable_exists(variable_name)) {
      throw RuntimeError(PISM_ERROR_LOCATION, "variable " + variable_name + " is missing");
    }

    auto dims = file.dimensions(variable_name);
    if (dims.size() != 1) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "variable '%s' in '%s' should to have 1 dimension (got %d)",
                                    variable_name.c_str(), file.name().c_str(), (int)dims.size());
    }

    const auto &dimension_name = dims[0];

    unsigned int length = file.dimension_length(dimension_name);
    if (length == 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "dimension " + dimension_name + " has length zero");
    }

    units::Unit internal_units(unit_system, units), input_units(unit_system, "1");

    auto input_units_string = file.read_text_attribute(variable_name, "units");

    if (not input_units_string.empty()) {
      input_units = units::Unit(unit_system, input_units_string);
    } else {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "variable '%s' does not have the units attribute", variable_name.c_str());
    }

    std::vector<double> result(length); // memory allocation happens here

    file.read_variable(variable_name, { 0 }, { length }, result.data());

    units::Converter(input_units, internal_units).convert_doubles(result.data(), result.size());

    return result;
  } catch (RuntimeError &e) {
    e.add_context("reading 1D variable '%s' from '%s'", variable_name.c_str(),
                  file.name().c_str());
    throw;
  }
}

std::vector<double> read_timeseries_variable(const File &file, const std::string &variable_name,
                                             const std::string &units,
                                             std::shared_ptr<units::System> unit_system,
                                             size_t start, size_t count) {

  auto input_units_string = file.read_text_attribute(variable_name, "units");

  if (input_units_string.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "variable '%s' does not have the units attribute",
                                  variable_name.c_str());
  }

  std::vector<double> result(count); // memory allocation happens here
  // read from a file
  {
    std::vector<unsigned int> start_{}, count_{};

    for (const auto &dim : file.dimensions(variable_name)) {
      if (file.dimension_type(dim, unit_system) == T_AXIS) {
        start_.push_back(start);
        count_.push_back(count);
      } else {
        start_.push_back(0);
        count_.push_back(1);
      }
    }
    file.read_variable(variable_name, start_, count_, result.data());
  }

  // convert units
  {
    auto input_units = units::Unit(unit_system, input_units_string);

    units::Unit internal_units(unit_system, units);
    units::Converter(input_units, internal_units).convert_doubles(result.data(), result.size());
  }

  return result;
}


std::vector<double> read_bounds(const File &file, const std::string &bounds_variable_name,
                                const std::string &internal_units,
                                std::shared_ptr<units::System> unit_system) {

  try {
    units::Unit internal(unit_system, internal_units);

    // Find the variable:
    if (not file.variable_exists(bounds_variable_name)) {
      throw RuntimeError(PISM_ERROR_LOCATION, "variable " + bounds_variable_name + " is missing");
    }

    auto dims = file.dimensions(bounds_variable_name);

    if (dims.size() != 2) {
      throw RuntimeError(PISM_ERROR_LOCATION, "variable " + bounds_variable_name + " has to have two dimensions");
    }

    const auto &dimension_name        = dims[0];
    const auto &bounds_dimension_name = dims[1];

    // Check that we have 2 vertices (interval end-points) per record.
    size_t length = file.dimension_length(bounds_dimension_name);
    if (length != 2) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "time-bounds variable " + bounds_variable_name + " has to have exactly 2 bounds per time record");
    }

    // Get the number of time records.
    length = file.dimension_length(dimension_name);
    if (length <= 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "dimension " + dimension_name + " has length zero");
    }

    // Find the corresponding coordinate variable. (We get units from the 'time'
    // variable, because according to CF-1.5 section 7.1 a "boundary variable"
    // may not have metadata set.)
    if (not file.variable_exists(dimension_name)) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "coordinate variable " + dimension_name + " is missing");
    }

    units::Unit input_units(unit_system, "1");

    std::string input_units_string = file.read_text_attribute(dimension_name, "units");
    if (input_units_string.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "variable '%s' does not have the units attribute",
                                    dimension_name.c_str());
    }

    input_units = units::Unit(unit_system, input_units_string);


    std::vector<double> result(length * 2); // memory allocation happens here

    file.read_variable(bounds_variable_name, {0, 0}, {(unsigned)length, 2}, result.data());

    units::Converter(input_units, internal).convert_doubles(result.data(), result.size());

    return result;
  } catch (RuntimeError &e) {
    e.add_context("reading bounds variable '%s' from '%s'", bounds_variable_name.c_str(),
                  file.name().c_str());
    throw;
  }
}

/*!
 * Reads and validates times and time bounds.
 */
void read_time_info(std::shared_ptr<units::System> unit_system, const File &file,
                    const std::string &time_name, const std::string &time_units,
                    std::vector<double> &times, std::vector<double> &bounds) {

  size_t N = 0;
  {
    times = io::read_1d_variable(file, time_name, time_units, unit_system);

    if (not is_increasing(times)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "times have to be strictly increasing");
    }
    N = times.size();
  }

  // Read time bounds
  {
    std::string time_bounds_name = file.read_text_attribute(time_name, "bounds");

    if (time_bounds_name.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "please provide time bounds for '%s'",
                                    time_name.c_str());
    }

    bounds = io::read_bounds(file, time_bounds_name, time_units, unit_system);

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

//! Read a variable from a file into an array `output`.
/*! This also converts data from input units to internal units if needed.
 */
void read_spatial_variable(const SpatialVariableMetadata &variable, const Grid &grid,
                           const File &file, unsigned int time, double *output) {

  const auto &log = *grid.ctx()->log();

  // Find the variable:
  auto var = file.find_variable(variable.get_name(), variable["standard_name"]);

  if (not var.exists) {
    throw RuntimeError::formatted(
        PISM_ERROR_LOCATION, "Can't find '%s' (%s) in '%s'.", variable.get_name().c_str(),
        variable.get_string("standard_name").c_str(), file.name().c_str());
  }

  // Sanity check: the variable in an input file should have the expected
  // number of spatial dimensions.
  {
    // Set of spatial dimensions for this field:
    std::set<int> axes{ X_AXIS, Y_AXIS };

    for (const auto &dim : variable.dimensions()) {
      if (axis_type_from_string(dim["axis"]) == Z_AXIS) {
        axes.insert(Z_AXIS);
        break;
      }
    }

    int input_spatial_dim_count = 0; // number of spatial dimensions (input file)
    size_t matching_dim_count   = 0; // number of matching dimensions

    auto input_dims = file.dimensions(var.name);
    for (const auto &d : input_dims) {
      auto dim_type = file.dimension_type(d, variable.unit_system());

      if (dim_type == X_AXIS or dim_type == Y_AXIS or dim_type == Z_AXIS) {
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
          file.name().c_str(), variable.get_name().c_str(),
          variable.get_string("long_name").c_str(), static_cast<int>(axes.size()));
    }
  }

  // make sure we have at least one level
  size_t nlevels{1};
  if (variable.levels() != nullptr) {
    nlevels = std::max(variable.levels()->size(), (size_t)1);
  }

  read_distributed_array(file, var.name, variable.unit_system(),
                         {(int)time, grid.xs(), grid.ys(), 0},
                         {1, grid.xm(), grid.ym(), (int)nlevels},
                         output);

  auto input_units           = file.read_text_attribute(var.name, "units");
  const std::string &internal_units = variable["units"];

  input_units = check_units(variable, input_units, log);

  // Convert data:
  size_t size = grid.xm() * grid.ym() * nlevels;

  units::Converter(variable.unit_system(), input_units, internal_units)
      .convert_doubles(output, size);
}

/*!
 * Check the overlap of the input grid and the internal grid.
 *
 * Set `allow_extrapolation` to `true` to "extend" the vertical grid during "bootstrapping".
 */
static void check_grid_overlap(const grid::InputGridInfo &input,
                               const grid::DistributedGridInfo &internal,
                               const std::vector<double> &z_internal) {

  // Grid spacing (assume that the grid is equally-spaced) and the
  // extent of the domain. To compute the extent of the domain, assume
  // that the grid is cell-centered, so edge of the domain is half of
  // the grid spacing away from grid points at the edge.
  //
  // Note that x_min is not the same as internal.x(0).
  const double x_min = internal.x0 - internal.Lx, x_max = internal.x0 + internal.Lx,
               y_min = internal.y0 - internal.Ly, y_max = internal.y0 + internal.Ly,
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
    const double input_z_min = vector_min(input.z), input_z_max = vector_max(input.z),
                 z_min = z_internal.front(), z_max = z_internal.back();

    if (not(z_min >= input_z_min - eps and z_max <= input_z_max + eps)) {
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
                      const grid::DistributedGridInfo &internal_grid,
                      const std::vector<double> &internal_z_levels, bool allow_extrapolation) {
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

  if (not allow_extrapolation) {
    check_grid_overlap(input_grid, internal_grid, internal_z_levels);
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

VariableMetadata read_attributes(const File &file,
                                 const std::string &variable_name,
                                 std::shared_ptr<units::System> unit_system) {

  VariableMetadata variable(variable_name, unit_system);

  try {

    if (not file.variable_exists(variable_name)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "variable '%s' is missing", variable_name.c_str());
    }

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
                  variable_name.c_str(), file.name().c_str());
    throw;
  }
  return variable;
}

} // namespace io
} // namespace pism

/* Copyright (C) 2015, 2016 PISM Authors
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

#include <memory>

#include <gsl/gsl_interp.h>

#include "io_helpers.hh"
#include "PIO.hh"
#include "base/util/IceGrid.hh"
#include "base/util/VariableMetadata.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/io/LocalInterpCtx.hh"
#include "base/util/PISMTime.hh"
#include "base/util/Logger.hh"
#include "base/util/projection.hh"

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
static void regrid(const IceGrid& grid, const std::vector<double> &zlevels_out,
                   LocalInterpCtx *lic, double *output_array) {
  // We'll work with the raw storage here so that the array we are filling is
  // indexed the same way as the buffer we are pulling from (input_array)

  const int X = 1, Z = 3; // indices, just for clarity

  std::vector<double> &zlevels_in = lic->zlevels;
  unsigned int nlevels = zlevels_out.size();
  double *input_array = &(lic->buffer[0]);

  std::vector<unsigned int> kbelow(nlevels, 0);

  if (nlevels > 1) {
    gsl_interp_accel *accel = gsl_interp_accel_alloc();
    if (accel == NULL) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Failed to allocate a GSL interpolation accelerator");
    }
    for (unsigned int k = 0; k < nlevels; ++k) {
      kbelow[k] = gsl_interp_accel_find(accel, &zlevels_in[0], zlevels_in.size(), zlevels_out[k]);
    }
    gsl_interp_accel_free(accel);
  }

  // array sizes for mapping from logical to "flat" indices
  int
    x_count = lic->count[X],
    z_count = lic->count[Z];

  for (Points p(grid); p; p.next()) {
    const int i_global = p.i(), j_global = p.j();

    const int i = i_global - grid.xs(), j = j_global - grid.ys();

    // Indices of neighboring points.
    const int
      X_m = lic->x_left[i],
      X_p = lic->x_right[i],
      Y_m = lic->y_left[j],
      Y_p = lic->y_right[j];

    for (unsigned int k = 0; k < nlevels; k++) {

      double a_mm, a_mp, a_pm, a_pp;  // filled differently in 2d and 3d cases

      if (nlevels > 1) {
        // get the index into the source grid, for just below the level z
        const int kc = kbelow[k];

        // We pretend that there are always 8 neighbors (4 in the map plane,
        // 2 vertical levels). And compute the indices into the input_array for
        // those neighbors.
        const int mmm = (Y_m * x_count + X_m) * z_count + kc;
        const int mmp = (Y_m * x_count + X_m) * z_count + kc + 1;
        const int mpm = (Y_m * x_count + X_p) * z_count + kc;
        const int mpp = (Y_m * x_count + X_p) * z_count + kc + 1;
        const int pmm = (Y_p * x_count + X_m) * z_count + kc;
        const int pmp = (Y_p * x_count + X_m) * z_count + kc + 1;
        const int ppm = (Y_p * x_count + X_p) * z_count + kc;
        const int ppp = (Y_p * x_count + X_p) * z_count + kc + 1;

        // We know how to index the neighbors, but we don't yet know where the
        // point lies within this box.  This is represented by alpha_z in [0,1].
        const double zkc = zlevels_in[kc];
        double dz;
        if (kc == z_count - 1) {
          dz = zlevels_in[kc] - zlevels_in[kc-1];
        } else {
          dz = zlevels_in[kc+1] - zlevels_in[kc];
        }
        // location (x,y,z) is in target computational domain
        const double z = zlevels_out[k];
        const double alpha_z = (z - zkc) / dz;

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
      const double x_alpha = lic->x_alpha[i];
      // interpolation coefficient in the y direction
      const double y_alpha = lic->y_alpha[j];

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

static void compute_start_and_count(const PIO& nc,
                                    units::System::Ptr unit_system,
                                    const std::string &short_name,
                                    unsigned int t_start, unsigned int t_count,
                                    unsigned int x_start, unsigned int x_count,
                                    unsigned int y_start, unsigned int y_count,
                                    unsigned int z_start, unsigned int z_count,
                                    std::vector<unsigned int> &start,
                                    std::vector<unsigned int> &count,
                                    std::vector<unsigned int> &imap) {
  std::vector<std::string> dims = nc.inq_vardims(short_name);
  unsigned int ndims = dims.size();

  assert(ndims > 0);
  if (ndims == 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Cannot compute start and count: variable %s is a scalar.",
                                  short_name.c_str());
  }

  // Resize output vectors:
  start.resize(ndims);
  count.resize(ndims);
  imap.resize(ndims);

  // Assemble start, count and imap:
  for (unsigned int j = 0; j < ndims; j++) {
    std::string dimname = dims[j];

    AxisType dimtype = nc.inq_dimtype(dimname, unit_system);

    switch (dimtype) {
    case T_AXIS:
      start[j] = t_start;
      count[j] = t_count;
      imap[j]  = x_count * y_count * z_count;
      break;
    case Y_AXIS:
      start[j] = y_start;
      count[j] = y_count;
      imap[j]  = x_count * z_count;
      break;
    case X_AXIS:
      start[j] = x_start;
      count[j] = x_count;
      imap[j]  = z_count;
      break;
    default:
    case Z_AXIS:
      start[j] = z_start;
      count[j] = z_count;
      imap[j]  = 1;
      break;
    }
  }
}

//! \brief Define a dimension \b and the associated coordinate variable. Set attributes.
void define_dimension(const PIO &nc, unsigned long int length,
                      const VariableMetadata &metadata) {
  std::string name = metadata.get_name();
  try {
    nc.redef();

    nc.def_dim(name, length);

    std::vector<std::string> dims(1, name);
    nc.def_var(name, PISM_DOUBLE, dims);

    // false means "don't use glaciological_units"
    write_attributes(nc, metadata, PISM_DOUBLE, false);

  } catch (RuntimeError &e) {
    e.add_context("defining dimension '%s' in '%s'", name.c_str(), nc.inq_filename().c_str());
    throw;
  }
}

void define_time(const PIO &nc, const std::string &name, const std::string &calendar,
                 const std::string &units, units::System::Ptr unit_system) {
  try {
    if (nc.inq_var(name)) {
      return;
    }

    // time
    VariableMetadata time(name, unit_system);
    time.set_string("long_name", "time");
    time.set_string("calendar", calendar);
    time.set_string("units", units);
    time.set_string("axis", "T");

    define_dimension(nc, PISM_UNLIMITED, time);
  } catch (RuntimeError &e) {
    e.add_context("defining the time dimension in \"" + nc.inq_filename() + "\"");
    throw;
  }
}

//! \brief Append to the time dimension.
void append_time(const PIO &nc, const std::string &name, double value) {
  try {
    unsigned int start = nc.inq_dimlen(name);
    std::vector<double> data(1, value);

    nc.put_1d_var(name, start, 1, data);
  } catch (RuntimeError &e) {
    e.add_context("appending to the time dimension in \"" + nc.inq_filename() + "\"");
    throw;
  }
}

//! \brief Define dimensions a variable depends on.
static void define_dimensions(const SpatialVariableMetadata& var,
                              const IceGrid& grid, const PIO &nc) {

  // x
  std::string x_name = var.get_x().get_name();
  if (not nc.inq_dim(x_name)) {
    define_dimension(nc, grid.Mx(), var.get_x());
    nc.put_1d_var(x_name, 0, grid.x().size(), grid.x());
  }

  // y
  std::string y_name = var.get_y().get_name();
  if (not nc.inq_dim(y_name)) {
    define_dimension(nc, grid.My(), var.get_y());
    nc.put_1d_var(y_name, 0, grid.y().size(), grid.y());
  }

  // z
  std::string z_name = var.get_z().get_name();
  if (not z_name.empty()) {
    if (not nc.inq_dim(z_name)) {
      const std::vector<double>& levels = var.get_levels();
      // make sure we have at least one level
      unsigned int nlevels = std::max(levels.size(), (size_t)1);
      define_dimension(nc, nlevels, var.get_z());
      nc.put_1d_var(z_name, 0, levels.size(), levels);
    }
  }
}

/**
 * Check if the storage order of a variable in the current file
 * matches the memory storage order used by PISM.
 *
 * @param var_name name of the variable to check
 * @returns false if storage orders match, true otherwise
 */
static bool use_mapped_io(const PIO &nc,
                          units::System::Ptr unit_system,
                          const std::string &var_name) {

  std::vector<std::string> dimnames = nc.inq_vardims(var_name);

  std::vector<AxisType> storage, memory = {Y_AXIS, X_AXIS};

  for (unsigned int j = 0; j < dimnames.size(); ++j) {
    AxisType dimtype = nc.inq_dimtype(dimnames[j], unit_system);

    if (j == 0 && dimtype == T_AXIS) {
      // ignore the time dimension, but only if it is the first
      // dimension in the list
      continue;
    }

    if (dimtype == X_AXIS || dimtype == Y_AXIS) {
      storage.push_back(dimtype);
    } else if (dimtype == Z_AXIS) {
      memory.push_back(dimtype); // now memory = {Y_AXIS, X_AXIS, Z_AXIS}
      // assume that this variable has only one Z_AXIS in the file
      storage.push_back(dimtype);
    } else {
      // an UNKNOWN_AXIS or T_AXIS at index != 0 was found, use mapped I/O
      return true;
    }
  }

  // we support 2D and 3D in-memory arrays, but not 4D
  assert(memory.size() <= 3);

  if (storage == memory) {
    // same storage order, do not use mapped I/O
    return false;
  } else {
    // different storage orders, use mapped I/O
    return true;
  }
}

//! \brief Read an array distributed according to the grid.
static void get_vec(const PIO &nc, const IceGrid &grid, const std::string &var_name,
                    unsigned int z_count, unsigned int t_start, double *output) {
  try {
    std::vector<unsigned int> start, count, imap;
    const unsigned int t_count = 1;
    compute_start_and_count(nc,
                            grid.ctx()->unit_system(),
                            var_name,
                            t_start, t_count,
                            grid.xs(), grid.xm(),
                            grid.ys(), grid.ym(),
                            0, z_count,
                            start, count, imap);

    bool mapped_io = use_mapped_io(nc, grid.ctx()->unit_system(), var_name);
    if (mapped_io == true) {
      nc.get_varm_double(var_name, start, count, imap, output);
    } else {
      nc.get_vara_double(var_name, start, count, output);
    }

  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' from '%s'", var_name.c_str(),
                  nc.inq_filename().c_str());
    throw;
  }
}

//! \brief Write an array distributed according to the grid.
/*!
 * This method always writes to the last record in the file.
 */
static void put_vec(const PIO &nc, const IceGrid &grid, const std::string &var_name,
                    unsigned int z_count, const double *input) {
  try {
    // switch to data mode and perform all delayed write operations
    nc.enddef();

    unsigned int t_length = nc.inq_dimlen(grid.ctx()->config()->get_string("time.dimension_name"));

    assert(t_length >= 1);

    std::vector<unsigned int> start, count, imap;
    const unsigned int t_count = 1;
    compute_start_and_count(nc,
                            grid.ctx()->unit_system(),
                            var_name,
                            t_length - 1, t_count,
                            grid.xs(), grid.xm(),
                            grid.ys(), grid.ym(),
                            0, z_count,
                            start, count, imap);

    if (grid.ctx()->config()->get_string("output.variable_order") == "yxz") {
      // Use the faster and safer (avoids a NetCDF bug) call if the array storage
      // orders in the memory and in NetCDF files are the same.
      nc.put_vara_double(var_name, start, count, input);
    } else {
      // Revert to "mapped" I/O otherwise.
      nc.put_varm_double(var_name, start, count, imap, input);
    }

  } catch (RuntimeError &e) {
    e.add_context("writing variable '%s' to '%s' in PIO::put_vec()",
                  var_name.c_str(), nc.inq_filename().c_str());
    throw;
  }
}


//! @brief Get the interpolation context (grid information) for an input file.
/*!
 * @note The *caller* is in charge of destroying lic
 */
static LocalInterpCtx* get_interp_context(const PIO& file,
                                          const std::string &variable_name,
                                          const IceGrid &grid,
                                          const std::vector<double> &zlevels) {
  grid_info gi(file, variable_name, grid.ctx()->unit_system(), grid.periodicity());

  return new LocalInterpCtx(gi, grid, zlevels.front(), zlevels.back());
}

//! \brief Read a PETSc Vec from a file, using bilinear (or trilinear)
//! interpolation to put it on the grid defined by "grid" and zlevels_out.
static void regrid_vec(const PIO &nc, const IceGrid &grid, const std::string &var_name,
                       const std::vector<double> &zlevels_out,
                       unsigned int t_start, double *output) {
  try {
    const int X = 1, Y = 2, Z = 3; // indices, just for clarity
    std::vector<unsigned int> start, count, imap;

    std::shared_ptr<LocalInterpCtx> lic(get_interp_context(nc, var_name, grid, zlevels_out));
    assert((bool)lic);

    double *buffer = &(lic->buffer[0]);

    const unsigned int t_count = 1;
    compute_start_and_count(nc,
                            grid.ctx()->unit_system(),
                            var_name,
                            t_start, t_count,
                            lic->start[X], lic->count[X],
                            lic->start[Y], lic->count[Y],
                            lic->start[Z], lic->count[Z],
                            start, count, imap);

    bool mapped_io = use_mapped_io(nc, grid.ctx()->unit_system(), var_name);
    if (mapped_io) {
      nc.get_varm_double(var_name, start, count, imap, buffer);
    } else {
      nc.get_vara_double(var_name, start, count, buffer);
    }

    // interpolate
    regrid(grid, zlevels_out, lic.get(), output);
  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' (using linear interpolation) from '%s'",
                  var_name.c_str(), nc.inq_filename().c_str());
    throw;
  }
}

/** Regrid `var_name` from a file, replacing missing values with `default_value`.
 *
 * @param grid computational grid; used to initialize interpolation
 * @param var_name variable to regrid
 * @param zlevels_out vertical levels of the resulting grid
 * @param t_start time index of the record to regrid
 * @param default_value default value to replace `_FillValue` with
 * @param[out] output resulting interpolated field
 */
static void regrid_vec_fill_missing(const PIO &nc, const IceGrid &grid,
                                    const std::string &var_name,
                                    const std::vector<double> &zlevels_out,
                                    unsigned int t_start,
                                    double default_value,
                                    double *output) {
  try {
    const int X = 1, Y = 2, Z = 3; // indices, just for clarity
    std::vector<unsigned int> start, count, imap;

    std::shared_ptr<LocalInterpCtx> lic(get_interp_context(nc, var_name, grid, zlevels_out));
    assert((bool)lic);

    double *buffer = &(lic->buffer[0]);

    const unsigned int t_count = 1;
    compute_start_and_count(nc,
                            grid.ctx()->unit_system(),
                            var_name,
                            t_start, t_count,
                            lic->start[X], lic->count[X],
                            lic->start[Y], lic->count[Y],
                            lic->start[Z], lic->count[Z],
                            start, count, imap);

    bool mapped_io = use_mapped_io(nc, grid.ctx()->unit_system(), var_name);
    if (mapped_io == true) {
      nc.get_varm_double(var_name, start, count, imap, buffer);
    } else {
      nc.get_vara_double(var_name, start, count, buffer);
    }

    // Replace missing values if the _FillValue attribute is present,
    // and if we have missing values to replace.
    {
      std::vector<double> attribute = nc.get_att_double(var_name, "_FillValue");
      if (attribute.size() == 1) {
        const double fill_value = attribute[0],
          epsilon = 1e-12;
        for (unsigned int i = 0; i < lic->buffer.size(); ++i) {
          if (fabs(buffer[i] - fill_value) < epsilon) {
            buffer[i] = default_value;
          }
        }
      }
    }

    // interpolate
    regrid(grid, zlevels_out, lic.get(), output);
  } catch (RuntimeError &e) {
    e.add_context("reading variable '%s' (using linear interpolation) from '%s'",
                  var_name.c_str(), nc.inq_filename().c_str());
    throw;
  }
}

//! Define a NetCDF variable corresponding to a VariableMetadata object.
void define_spatial_variable(const SpatialVariableMetadata &var,
                             const IceGrid &grid, const PIO &nc,
                             IO_Type nctype,
                             const std::string &variable_order,
                             bool use_glaciological_units) {
  std::vector<std::string> dims;
  std::string name = var.get_name();

  if (nc.inq_var(name)) {
    return;
  }

  nc.set_local_extent(grid.xs(), grid.xm(), grid.ys(), grid.ym());

  define_dimensions(var, grid, nc);

  std::string order = variable_order;
  // "..._bounds" should be stored with grid corners (corresponding to
  // the "z" dimension here) last, so we override the variable storage
  // order here
  if (ends_with(name, "_bounds") and order == "zyx") {
    order = "yxz";
  }

  std::string
    x = var.get_x().get_name(),
    y = var.get_y().get_name(),
    z = var.get_z().get_name(),
    t = "";

  if (not var.get_time_independent()) {
    t = grid.ctx()->config()->get_string("time.dimension_name");
  }

  nc.redef();

  if (not t.empty()) {
    dims.push_back(t);
  }

  // Use t,x,y,z(zb) variable order: it is weird, but matches the in-memory
  // storage order and so is *a lot* faster.
  if (order == "xyz") {
    dims.push_back(x);
    dims.push_back(y);
  }

  // Use the t,y,x,z variable order: also weird, somewhat slower, but 2D fields
  // are stored in the "natural" order.
  if (order == "yxz") {
    dims.push_back(y);
    dims.push_back(x);
  }

  if (not z.empty()) {
    dims.push_back(z);
  }

  // Use the t,z(zb),y,x variables order: more natural for plotting and post-processing,
  // but requires transposing data while writing and is *a lot* slower.
  if (order == "zyx") {
    dims.push_back(y);
    dims.push_back(x);
  }

  assert(dims.size() > 1);

  nc.def_var(name, nctype, dims);

  write_attributes(nc, var, nctype, use_glaciological_units);

  // add the "grid_mapping" attribute if the grid has an associated mapping.
  const VariableMetadata &mapping = grid.get_mapping_info().mapping;
  if (mapping.has_attributes()) {
    nc.put_att_text(var.get_name(), "grid_mapping",
                    mapping.get_name());
  }
}

//! Read a variable from a file into an array `output`.
/*! This also converts data from input units to internal units if needed.
 */
void read_spatial_variable(const SpatialVariableMetadata &var,
                           const IceGrid& grid, const PIO &nc,
                           unsigned int time, double *output) {

  const Logger &log = *grid.ctx()->log();

  nc.set_local_extent(grid.xs(), grid.xm(), grid.ys(), grid.ym());

  // Find the variable:
  std::string name_found;
  bool found_by_standard_name = false, variable_exists = false;
  nc.inq_var(var.get_name(), var.get_string("standard_name"),
             variable_exists, name_found, found_by_standard_name);

  if (not variable_exists) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' (%s) in '%s'.",
                                  var.get_name().c_str(),
                                  var.get_string("standard_name").c_str(), nc.inq_filename().c_str());
  }

  // Sanity check: the variable in an input file should have the expected
  // number of spatial dimensions.
  {
    // Set of spatial dimensions this field has.
    std::set<int> axes;
    axes.insert(X_AXIS);
    axes.insert(Y_AXIS);
    if (not var.get_z().get_name().empty()) {
      axes.insert(Z_AXIS);
    }

    std::vector<std::string> input_dims;
    int input_ndims = 0;                 // number of spatial dimensions (input file)
    size_t matching_dim_count = 0; // number of matching dimensions

    input_dims = nc.inq_vardims(name_found);
    for (auto d : input_dims) {
      AxisType tmp = nc.inq_dimtype(d, var.unit_system());

      if (tmp != T_AXIS) {
        ++input_ndims;
      }

      if (axes.find(tmp) != axes.end()) {
        ++matching_dim_count;
      }
    }

    if (axes.size() != matching_dim_count) {

      // Print the error message and stop:
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "found the %dD variable %s (%s) in '%s' while trying to read\n"
                                    "'%s' ('%s'), which is %d-dimensional.",
                                    input_ndims, name_found.c_str(),
                                    join(input_dims, ",").c_str(),
                                    nc.inq_filename().c_str(),
                                    var.get_name().c_str(), var.get_string("long_name").c_str(),
                                    static_cast<int>(axes.size()));
    }
  }

  // make sure we have at least one level
  const std::vector<double>& zlevels = var.get_levels();
  unsigned int nlevels = std::max(zlevels.size(), (size_t)1);

  get_vec(nc, grid, name_found, nlevels, time, output);

  std::string input_units = nc.get_att_text(name_found, "units");

  if (var.has_attribute("units") && input_units.empty()) {
    const std::string &units_string = var.get_string("units"),
      &long_name = var.get_string("long_name");
    log.message(2,
               "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
               "              Assuming that it is in '%s'.\n",
               var.get_name().c_str(), long_name.c_str(),
               units_string.c_str());
    input_units = units_string;
  }

  // Convert data:
  size_t size = grid.xm() * grid.ym() * nlevels;

  units::Converter(var.unit_system(),
                input_units,
                var.get_string("units")).convert_doubles(output, size);
}

//! \brief Write a double array to a file.
/*!
  Converts the units if `use_glaciological_units` is `true`.
 */
void write_spatial_variable(const SpatialVariableMetadata &var,
                            const IceGrid& grid,
                            const PIO &nc, bool use_glaciological_units,
                            const double *input) {

  nc.set_local_extent(grid.xs(), grid.xm(), grid.ys(), grid.ym());

  // find or define the variable
  std::string name_found;
  bool exists, found_by_standard_name;
  nc.inq_var(var.get_name(),
             var.get_string("standard_name"),
             exists, name_found, found_by_standard_name);

  if (not exists) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' in '%s'.",
                                  var.get_name().c_str(),
                                  nc.inq_filename().c_str());
  }

  // make sure we have at least one level
  const std::vector<double>& zlevels = var.get_levels();
  unsigned int nlevels = std::max(zlevels.size(), (size_t)1);

  if (use_glaciological_units) {
    size_t data_size = grid.xm() * grid.ym() * nlevels;

    // create a temporary array, convert to glaciological units, and
    // save
    std::vector<double> tmp(data_size);
    for (size_t k = 0; k < data_size; ++k) {
      tmp[k] = input[k];
    }

    units::Converter(var.unit_system(),
                  var.get_string("units"),
                  var.get_string("glaciological_units")).convert_doubles(&tmp[0], tmp.size());

    put_vec(nc, grid, name_found, nlevels, &tmp[0]);
  } else {
    put_vec(nc, grid, name_found, nlevels, input);
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
void regrid_spatial_variable(SpatialVariableMetadata &var,
                             const IceGrid& grid, const PIO &nc,
                             RegriddingFlag flag, bool do_report_range,
                             double default_value, double *output) {
  unsigned int t_length = nc.inq_nrecords(var.get_name(),
                                          var.get_string("standard_name"),
                                          var.unit_system());

  regrid_spatial_variable(var, grid, nc, t_length - 1, flag, do_report_range,
                          default_value, output);
}

static void compute_range(MPI_Comm com, double *data, size_t data_size, double *min, double *max) {
  double
    min_result = data[0],
    max_result = data[0];
  for (size_t k = 0; k < data_size; ++k) {
    min_result = std::min(min_result, data[k]);
    max_result = std::max(max_result, data[k]);
  }

  if (min) {
    *min = GlobalMin(com, min_result);
  }

  if (max) {
    *max = GlobalMax(com, max_result);
  }
}

void regrid_spatial_variable(SpatialVariableMetadata &var,
                             const IceGrid& grid, const PIO &nc,
                             unsigned int t_start, RegriddingFlag flag,
                             bool do_report_range, double default_value,
                             double *output) {
  const Logger &log = *grid.ctx()->log();

  nc.set_local_extent(grid.xs(), grid.xm(), grid.ys(), grid.ym());

  units::System::Ptr sys = var.unit_system();
  const std::vector<double>& levels = var.get_levels();
  const size_t data_size = grid.xm() * grid.ym() * levels.size();

  // Find the variable
  bool exists, found_by_standard_name;
  std::string name_found;
  nc.inq_var(var.get_name(), var.get_string("standard_name"),
             exists, name_found, found_by_standard_name);

  if (exists) {                      // the variable was found successfully

    if (flag == OPTIONAL_FILL_MISSING or flag == CRITICAL_FILL_MISSING) {
      log.message(2,
                 "PISM WARNING: Replacing missing values with %f [%s] in variable '%s' read from '%s'.\n",
                 default_value, var.get_string("units").c_str(), var.get_name().c_str(),
                 nc.inq_filename().c_str());

      regrid_vec_fill_missing(nc, grid, name_found, levels,
                              t_start, default_value, output);
    } else {
      regrid_vec(nc, grid, name_found, levels, t_start, output);
    }

    // Now we need to get the units string from the file and convert
    // the units, because check_range and report_range expect data to
    // be in PISM (MKS) units.

    std::string input_units = nc.get_att_text(name_found, "units");

    if (input_units.empty()) {
      std::string internal_units = var.get_string("units");
      input_units = internal_units;
      if (not internal_units.empty()) {
        log.message(2,
                   "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                   "              Assuming that it is in '%s'.\n",
                   var.get_name().c_str(),
                   var.get_string("long_name").c_str(),
                   internal_units.c_str());
      }
    }

    // Convert data:
    units::Converter(sys,
                  input_units,
                  var.get_string("units")).convert_doubles(output, data_size);

    // Read the valid range info:
    read_valid_range(nc, name_found, var);

    double min = 0.0, max = 0.0;
    compute_range(nc.com(), output, data_size, &min, &max);

    // Check the range and warn the user if needed:
    var.check_range(nc.inq_filename(), min, max);

    if (do_report_range) {
      // We can report the success, and the range now:
      log.message(2, "  FOUND ");

      var.report_range(log, min, max, found_by_standard_name);
    }
  } else {                // couldn't find the variable
    if (flag == CRITICAL or flag == CRITICAL_FILL_MISSING) {
      // if it's critical, print an error message and stop
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' in the regridding file '%s'.",
                                    var.get_name().c_str(), nc.inq_filename().c_str());
    }

    // If it is optional, fill with the provided default value.
    // units::Converter constructor will make sure that units are compatible.
    units::Converter c(sys,
                    var.get_string("units"),
                    var.get_string("glaciological_units"));

    std::string spacer(var.get_name().size(), ' ');
    log.message(2,
               "  absent %s / %-10s\n"
               "         %s \\ not found; using default constant %7.2f (%s)\n",
               var.get_name().c_str(),
               var.get_string("long_name").c_str(),
               spacer.c_str(), c(default_value),
               var.get_string("glaciological_units").c_str());

    for (size_t k = 0; k < data_size; ++k) {
      output[k] = default_value;
    }
  } // end of if (exists)
}



//! Define a NetCDF variable corresponding to a time-series.
void define_timeseries(const TimeseriesMetadata& var,
                       const PIO &nc, IO_Type nctype, bool) {

  std::string name = var.get_name();
  std::string dimension_name = var.get_dimension_name();

  if (nc.inq_var(name)) {
    return;
  }

  if (not nc.inq_dim(dimension_name)) {
    define_dimension(nc, PISM_UNLIMITED,
                     VariableMetadata(dimension_name, var.unit_system()));
  }

  if (not nc.inq_var(name)) {
    std::vector<std::string> dims(1, dimension_name);
    nc.redef();
    nc.def_var(name, nctype, dims);
  }

  write_attributes(nc, var, PISM_FLOAT, true);
}

//! Read a time-series variable from a NetCDF file to a vector of doubles.
void read_timeseries(const PIO &nc, const TimeseriesMetadata &metadata,
                     const Time &time, const Logger &log, std::vector<double> &data) {

  std::string name = metadata.get_name();

  try {
    bool variable_exists = false;

    // Find the variable:
    std::string name_found,
      long_name      = metadata.get_string("long_name"),
      standard_name  = metadata.get_string("standard_name"),
      dimension_name = metadata.get_dimension_name();

    bool found_by_standard_name = false;
    nc.inq_var(name, standard_name, variable_exists, name_found, found_by_standard_name);

    if (not variable_exists) {
      throw RuntimeError(PISM_ERROR_LOCATION, "variable " + name + " is missing");
    }

    std::vector<std::string> dims = nc.inq_vardims(name_found);
    if (dims.size() != 1) {
      throw RuntimeError(PISM_ERROR_LOCATION, "a time-series variable has to be one-dimensional");
    }

    unsigned int length = nc.inq_dimlen(dimension_name);
    if (length <= 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "dimension " + dimension_name + " has length zero");
    }

    data.resize(length);          // memory allocation happens here

    nc.get_1d_var(name_found, 0, length, data);

    units::System::Ptr system = metadata.unit_system();
    bool input_has_units = false;
    units::Unit internal_units(system, metadata.get_string("units")),
      input_units(system, "1");

    std::string input_units_string = nc.get_att_text(name_found, "units");

    if (input_units_string.empty()) {
      input_has_units = false;
    } else {
      input_units_string = time.CF_units_to_PISM_units(input_units_string);

      input_units = units::Unit(system, input_units_string);
      input_has_units = true;
    }

    if (metadata.has_attribute("units") && not input_has_units) {
      std::string units_string = internal_units.format();
      log.message(2,
                 "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                 "              Assuming that it is in '%s'.\n",
                 name.c_str(), long_name.c_str(),
                 units_string.c_str());
      input_units = internal_units;
    }

    units::Converter(input_units, internal_units).convert_doubles(&data[0], data.size());

  } catch (RuntimeError &e) {
    e.add_context("reading time-series variable '%s' from '%s'", name.c_str(),
                  nc.inq_filename().c_str());
    throw;
  }
}

void write_timeseries(const PIO &nc, const TimeseriesMetadata &metadata, size_t t_start,
                      double data, IO_Type nctype) {
  std::vector<double> vector_data(1, data);

  // this call will handle errors
  write_timeseries(nc, metadata, t_start, vector_data, nctype);
}

/** @brief Write a time-series `data` to a file.
 *
 * Always use glaciological units when saving time-series.
 */
void write_timeseries(const PIO &nc, const TimeseriesMetadata &metadata, size_t t_start,
                      const std::vector<double> &data,
                      IO_Type nctype) {

  std::string name = metadata.get_name();
  try {
    if (not nc.inq_var(name)) {
      define_timeseries(metadata, nc, nctype, true);
    }

    // create a copy of "data":
    std::vector<double> tmp = data;

    units::System::Ptr system = metadata.unit_system();
    // convert to glaciological units:
    units::Converter(system,
                  metadata.get_string("units"),
                  metadata.get_string("glaciological_units")).convert_doubles(&tmp[0], tmp.size());

    nc.put_1d_var(name,
                  static_cast<unsigned int>(t_start),
                  static_cast<unsigned int>(tmp.size()), tmp);

  } catch (RuntimeError &e) {
    e.add_context("writing time-series variable '%s' to '%s'", name.c_str(),
                  nc.inq_filename().c_str());
    throw;
  }
}

void define_time_bounds(const TimeBoundsMetadata& var,
                        const PIO &nc, IO_Type nctype, bool) {
  std::string name = var.get_name();
  std::string dimension_name = var.get_dimension_name();
  std::string bounds_name = var.get_bounds_name();

  if (nc.inq_var(name)) {
    return;
  }

  nc.redef();

  if (not nc.inq_dim(dimension_name)) {
    nc.def_dim(dimension_name, PISM_UNLIMITED);
  }

  if (not nc.inq_dim(bounds_name)) {
    nc.def_dim(bounds_name, 2);
  }

  nc.redef();

  nc.def_var(name, nctype, {dimension_name, bounds_name});

  write_attributes(nc, var, nctype, true);
}

void read_time_bounds(const PIO &nc,
                      const TimeBoundsMetadata &metadata,
                      const Time &time, const Logger &log,
                      std::vector<double> &data) {

  std::string name = metadata.get_name();

  try {
    units::System::Ptr system = metadata.unit_system();
    units::Unit internal_units(system, metadata.get_string("units"));

    // Find the variable:
    if (not nc.inq_var(name)) {
      throw RuntimeError(PISM_ERROR_LOCATION, "variable " + name + " is missing");
    }

    std::vector<std::string> dims = nc.inq_vardims(name);

    if (dims.size() != 2) {
      throw RuntimeError(PISM_ERROR_LOCATION, "variable " + name + " has to has two dimensions");
    }

    std::string
      &dimension_name = dims[0],
      &bounds_name    = dims[1];

    // Check that we have 2 vertices (interval end-points) per time record.
    unsigned int length = nc.inq_dimlen(bounds_name);
    if (length != 2) {
      throw RuntimeError(PISM_ERROR_LOCATION, "time-bounds variable " + name + " has to have exactly 2 bounds per time record");
    }

    // Get the number of time records.
    length = nc.inq_dimlen(dimension_name);
    if (length <= 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "dimension " + dimension_name + " has length zero");
    }

    data.resize(2*length);                // memory allocation happens here

    std::vector<unsigned int> start(2), count(2);
    start[0] = 0;
    start[1] = 0;
    count[0] = length;
    count[1] = 2;

    nc.get_vara_double(name, start, count, &data[0]);

    // Find the corresponding 'time' variable. (We get units from the 'time'
    // variable, because according to CF-1.5 section 7.1 a "boundary variable"
    // may not have metadata set.)
    if (not nc.inq_var(dimension_name)) {
      throw RuntimeError(PISM_ERROR_LOCATION, "time coordinate variable " + dimension_name + " is missing");
    }

    bool input_has_units = false;
    units::Unit input_units(internal_units.get_system(), "1");

    std::string input_units_string = nc.get_att_text(dimension_name, "units");
    input_units_string = time.CF_units_to_PISM_units(input_units_string);

    if (input_units_string.empty() == true) {
      input_has_units = false;
    } else {
      input_units = units::Unit(internal_units.get_system(), input_units_string);
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

    units::Converter(input_units, internal_units).convert_doubles(&data[0], data.size());

    // FIXME: check that time intervals described by the time bounds
    // variable are contiguous (without gaps) and stop if they are not.
  } catch (RuntimeError &e) {
    e.add_context("reading time bounds variable '%s' from '%s'", name.c_str(),
                  nc.inq_filename().c_str());
    throw;
  }
}

void write_time_bounds(const PIO &nc, const TimeBoundsMetadata &metadata,
                       size_t t_start, const std::vector<double> &data, IO_Type nctype) {
  std::string name = metadata.get_name();
  try {
    bool variable_exists = nc.inq_var(name);
    if (not variable_exists) {
      define_time_bounds(metadata, nc, nctype, true);
    }

    // make a copy of "data"
    std::vector<double> tmp = data;

    // convert to glaciological units:
    units::System::Ptr system = metadata.unit_system();
    units::Converter(system,
                  metadata.get_string("units"),
                  metadata.get_string("glaciological_units")).convert_doubles(&tmp[0], tmp.size());

    std::vector<unsigned int> start(2), count(2);
    start[0] = static_cast<unsigned int>(t_start);
    start[1] = 0;
    count[0] = static_cast<unsigned int>(tmp.size()) / 2;
    count[1] = 2;

    nc.enddef();
    nc.put_vara_double(name, start, count, &tmp[0]);

  } catch (RuntimeError &e) {
    e.add_context("writing time-bounds variable '%s' to '%s'", name.c_str(),
                  nc.inq_filename().c_str());
    throw;
  }
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

  if (file_exists_flag == 1) {
    return true;
  } else {
    return false;
  }
}

void read_attributes(const PIO &nc, const std::string &variable_name, VariableMetadata &variable) {
  try {
    bool variable_exists = nc.inq_var(variable_name);

    if (not variable_exists) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "variable \"%s\" is missing", variable_name.c_str());
    }

    variable.clear_all_strings();
    variable.clear_all_doubles();

    unsigned int nattrs = nc.inq_nattrs(variable_name);

    for (unsigned int j = 0; j < nattrs; ++j) {
      std::string attribute_name = nc.inq_attname(variable_name, j);
      IO_Type nctype = nc.inq_atttype(variable_name, attribute_name);

      if (nctype == PISM_CHAR) {
        std::string value = nc.get_att_text(variable_name, attribute_name);

        if (attribute_name == "units") {
          variable.set_string("units", value);
        } else {
          variable.set_string(attribute_name, value);
        }
      } else {
        variable.set_doubles(attribute_name,
                             nc.get_att_double(variable_name, attribute_name));
      }
    } // end of for (int j = 0; j < nattrs; ++j)
  } catch (RuntimeError &e) {
    e.add_context("reading attributes of variable '%s' from '%s'",
                  variable_name.c_str(), nc.inq_filename().c_str());
    throw;
  }
}

//! Write variable attributes to a NetCDF file.
/*!
  - If use_glaciological_units == true, "glaciological_units" are
    written under the name "units" plus the valid range is written in
    glaciological units.

  - If both valid_min and valid_max are set, then valid_range is written
    instead of the valid_min, valid_max pair.

  - Skips empty text attributes.
 */
void write_attributes(const PIO &nc, const VariableMetadata &variable, IO_Type nctype,
                      bool use_glaciological_units) {
  std::string var_name = variable.get_name();
  try {
    // units, valid_min, valid_max and valid_range need special treatment:
    if (variable.has_attribute("units")) {
      std::string output_units = variable.get_string("units");

      if (use_glaciological_units) {
        output_units = variable.get_string("glaciological_units");
      }

      nc.put_att_text(var_name, "units", output_units);
    }

    std::vector<double> bounds(2);
    if (variable.has_attribute("valid_min")) {
      bounds[0]  = variable.get_double("valid_min");
    }
    if (variable.has_attribute("valid_max")) {
      bounds[1]  = variable.get_double("valid_max");
    }

    double fill_value = 0.0;
    if (variable.has_attribute("_FillValue")) {
      fill_value = variable.get_double("_FillValue");
    }

    // We need to save valid_min, valid_max and valid_range in the units
    // matching the ones in the output.
    if (use_glaciological_units) {

      std::string
        units               = variable.get_string("units"),
        glaciological_units = variable.get_string("glaciological_units");
      units::Converter c(variable.unit_system(), units, glaciological_units);

      bounds[0]  = c(bounds[0]);
      bounds[1]  = c(bounds[1]);
      fill_value = c(fill_value);
    }

    if (variable.has_attribute("_FillValue")) {
      nc.put_att_double(var_name, "_FillValue", nctype, fill_value);
    }

    if (variable.has_attribute("valid_min") && variable.has_attribute("valid_max")) {
      nc.put_att_double(var_name, "valid_range", nctype, bounds);
    } else if (variable.has_attribute("valid_min")) {
      nc.put_att_double(var_name, "valid_min",   nctype, bounds[0]);
    } else if (variable.has_attribute("valid_max")) {
      nc.put_att_double(var_name, "valid_max",   nctype, bounds[1]);
    }

    // Write text attributes:
    for (auto s : variable.get_all_strings()) {
      std::string
        name  = s.first,
        value = s.second;

      if (name == "units" or
          name == "glaciological_units" or
          value.empty()) {
        continue;
      }

      nc.put_att_text(var_name, name, value);
    }

    // Write double attributes:
    for (auto d : variable.get_all_doubles()) {
      std::string name  = d.first;
      std::vector<double> values = d.second;

      if (name == "valid_min"   or
          name == "valid_max"   or
          name == "valid_range" or
          name == "_FillValue"  or
          values.empty()) {
        continue;
      }

      nc.put_att_double(var_name, name, nctype, values);
    }

  } catch (RuntimeError &e) {
    e.add_context("writing attributes of variable '%s' to '%s'",
                  var_name.c_str(), nc.inq_filename().c_str());
    throw;
  }
}

/** Write global attributes to a file.
 *
 * Same as `write_attributes(nc, var, PISM_DOUBLE, false)`, but
 * prepends the history string.
 *
 * @param nc file to write to
 * @param var metadata object containing attributes
 *
 * @return 0 on success
 */
void write_global_attributes(const PIO &nc, const VariableMetadata &var) {
  try {
    VariableMetadata tmp = var;

    std::string old_history = nc.get_att_text("PISM_GLOBAL", "history");

    tmp.set_name("PISM_GLOBAL");
    tmp.set_string("history", tmp.get_string("history") + old_history);

    write_attributes(nc, tmp, PISM_DOUBLE, false);

  } catch (RuntimeError &e) {
    e.add_context("writing global attributes to \"" + nc.inq_filename());
    throw;
  }
}

//! Read the valid range information from a file.
/*! Reads `valid_min`, `valid_max` and `valid_range` attributes; if \c
    valid_range is found, sets the pair `valid_min` and `valid_max` instead.
 */
void read_valid_range(const PIO &nc, const std::string &name, VariableMetadata &variable) {
  try {
    // Never reset valid_min/max if they were set internally
    if (variable.has_attribute("valid_min") or
        variable.has_attribute("valid_max")) {
      return;
    }

    // Read the units.
    units::Converter c(variable.unit_system(),
                    nc.get_att_text(name, "units"),
                    variable.get_string("units"));

    std::vector<double> bounds = nc.get_att_double(name, "valid_range");
    if (bounds.size() == 2) {             // valid_range is present
      variable.set_double("valid_min", c(bounds[0]));
      variable.set_double("valid_max", c(bounds[1]));
    } else {                      // valid_range has the wrong length or is missing
      bounds = nc.get_att_double(name, "valid_min");
      if (bounds.size() == 1) {           // valid_min is present
        variable.set_double("valid_min", c(bounds[0]));
      }

      bounds = nc.get_att_double(name, "valid_max");
      if (bounds.size() == 1) {           // valid_max is present
        variable.set_double("valid_max", c(bounds[0]));
      }
    }
  } catch (RuntimeError &e) {
    e.add_context("reading valid range of variable '%s' from '%s'", name.c_str(),
                  nc.inq_filename().c_str());
    throw;
  }
}

} // end of namespace io
} // end of namespace pism

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

#include <cstddef>
#include <map>
#include <memory>
#include <mpi.h>
#include <vector>
#include <cassert>

#include "pism/util/Config.hh"
#include "pism/util/GridInfo.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/OutputWriter.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {

/*! Pre-process variable attributes for writing.
*/
static VariableMetadata format_attributes(const VariableMetadata &metadata) {

  // make a copy so we can edit metadata
  auto variable = metadata;

  std::string units = variable["units"], output_units = variable["output_units"];

  // output units should never be written to a file
  variable["output_units"] = "";

  bool use_output_units =
      (not units.empty() and not output_units.empty() and units != output_units);

  if (use_output_units) {
    // Replace "units" with "output_units"
    if (variable.has_attribute("units")) {
      variable["units"] = output_units;
    }

    units::Converter c(variable.unit_system(), units, output_units);

    // We need to convert units of valid_min, valid_max and valid_range:
    {
      if (variable.has_attribute("valid_range")) {
        auto bounds = variable.get_numbers("valid_range");

        variable["valid_range"] = { c(bounds[0]), c(bounds[1]) };
      } else {
        if (variable.has_attribute("valid_min")) {
          auto min              = variable.get_number("valid_min");
          variable["valid_min"] = { c(min) };
        }
        if (variable.has_attribute("valid_max")) {
          auto max              = variable.get_number("valid_max");
          variable["valid_max"] = { c(max) };
        }
      }

      if (variable.has_attribute("_FillValue")) {
        auto fill              = variable.get_number("_FillValue");
        variable["_FillValue"] = { c(fill) };
      }
    }
  }
  return variable;
}

struct OutputWriter::Impl {
  Impl(MPI_Comm comm_, const Config &config)
      : comm(comm_) {
    time_name          = config.get_string("time.dimension_name");
    use_internal_units = config.get_flag("output.use_MKS");
  }

  std::string time_name;
  MPI_Comm comm;
  std::map<std::string, std::map<std::string, std::string> > extra_attributes;
  bool use_internal_units;
  std::map<std::tuple<std::string, std::string>, bool> written;
  std::map<std::string, grid::DistributedGridInfo> grids;
};

bool &OutputWriter::already_written(const std::string &file_name,
                                    const std::string &variable_name) {
  return m_impl->written[{ file_name, variable_name }];
}

OutputWriter::OutputWriter(MPI_Comm comm, const Config &config)
    : m_impl(new Impl(comm, config)) {
}

void OutputWriter::add_extra_attributes(const std::string &file_name,
                                        const std::map<std::string, std::string> &attributes) {
  for (const auto &attr : attributes) {
    m_impl->extra_attributes[file_name][attr.first] = attr.second;
  }
}

OutputWriter::~OutputWriter() {
  delete m_impl;
}

MPI_Comm OutputWriter::comm() const {
  return m_impl->comm;
}

const std::string &OutputWriter::time_name() const {
  return m_impl->time_name;
}

const grid::DistributedGridInfo &OutputWriter::grid_info(const std::string &variable_name) const {
  return m_impl->grids[variable_name];
}

void OutputWriter::define_dimension(const std::string &file_name, const std::string &dimension_name,
                                    unsigned int length) {
  define_dimension_impl(file_name, dimension_name, length);
}

void OutputWriter::define_variable(const std::string &file_name, const VariableMetadata &metadata,
                                   const std::vector<std::string> &dims) {
  define_variable_impl(file_name, format_attributes(metadata), dims);
}

void OutputWriter::define_spatial_variable(const std::string &file_name,
                                           const SpatialVariableMetadata &metadata,
                                           const grid::DistributedGridInfo &grid) {

  // Make a copy of `metadata` so we can modify it:
  auto var               = metadata;
  const auto &name       = var.get_name();
  const auto &extra_attributes = m_impl->extra_attributes[file_name];

  m_impl->grids[name] = grid;

  if (m_impl->use_internal_units) {
    var.output_units(var["units"]);
  }

  // add extra attributes such as "grid_mapping" and "coordinates". Variables lat, lon,
  // lat_bnds, and lon_bnds should not have these attributes to support CDO (see issue
  // #384).
  if (metadata.x().get_name() == "x" and metadata.y().get_name() == "y") {
    if (not member(name, { "lat_bnds", "lon_bnds", "lat", "lon" })) {
      for (const auto &attr : extra_attributes) {
        const auto &attr_name  = attr.first;
        const auto &attr_value = attr.second;
        var[attr_name]         = attr_value;
      }
    }
  }

  std::vector<std::string> dims;

  if (not var.get_time_independent()) {
    dims = { m_impl->time_name };
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
    define_dimension(file_name, dimension_name, dimension.length());
    define_variable(file_name, dimension, { dimension_name });
  }

  assert(dims.size() > 1);

  // define the variable itself:
  define_variable(file_name, var, dims);
}

void OutputWriter::set_global_attributes(
    const std::string &file_name, const std::map<std::string, std::string> &strings,
    const std::map<std::string, std::vector<double> > &numbers) {
  set_global_attributes_impl(file_name, strings, numbers);
}

void OutputWriter::append_time(const std::string &file_name, double time_seconds) {
  append_time_impl(file_name, time_seconds);
}

void OutputWriter::append_history(const std::string &file_name, const std::string &text) {
  append_history_impl(file_name, text);
}

void OutputWriter::write_array(const std::string &file_name, const std::string &variable_name,
                               const std::vector<unsigned int> &start,
                               const std::vector<unsigned int> &count,
                               const std::vector<double> &input) {
  write_array_impl(file_name, variable_name, start, count, input.data());
}

void OutputWriter::write_array(const std::string &file_name, const VariableMetadata &metadata,
                               const std::vector<unsigned int> &start,
                               const std::vector<unsigned int> &count,
                               const std::vector<double> &input) {
  // create a copy of "data" to change units
  std::vector<double> data = input;

  // convert units "in place":
  units::Converter(metadata.unit_system(), metadata["units"], metadata["output_units"])
      .convert_doubles(data.data(), data.size());

  write_array(file_name, metadata.get_name(), start, count, data);
}

void OutputWriter::write_spatial_variable(const std::string &file_name,
                                          const SpatialVariableMetadata &metadata,
                                          const double *input) {
  const auto &variable_name = metadata.get_name();
  const auto &grid          = grid_info(variable_name);

  // check if we need to write this variable
  bool time_independent = metadata.get_time_independent();
  // avoid writing time-independent variables more than once (saves time when writing to
  // extra_files)
  if (time_independent and already_written(file_name, variable_name)) {
    return;
  }

  // write dimensions:
  {
    std::map<std::string, const std::vector<double> &> data = { { metadata.x().get_name(), grid.x },
                                                                { metadata.y().get_name(), grid.y },
                                                                { metadata.z().get_name(),
                                                                  metadata.levels() } };
    for (const auto &p : data) {
      const auto &dimension_name = p.first;
      const auto &coordinates    = p.second;

      if (dimension_name.empty() or coordinates.empty()) {
        continue;
      }

      if (not already_written(file_name, dimension_name)) {
        write_array(file_name, dimension_name, { 0 }, { (unsigned int)coordinates.size() },
                    coordinates);
        already_written(file_name, dimension_name) = true;
      }
    }
  }

  // make sure we have at least one level
  unsigned int nlevels = std::max(metadata.levels().size(), (std::size_t)1);

  std::string units = metadata["units"], output_units = metadata["output_units"];

  // override `output_units` if "output.use_MKS" is set
  if (m_impl->use_internal_units) {
    output_units = units;
  }

  std::vector<unsigned int> start, count;

  if (time_independent) {
    start = { grid.ys, grid.xs, 0 };
    count = { grid.ym, grid.xm, nlevels };
  } else {
    auto t_length = time_dimension_length(file_name);
    auto t_start  = t_length > 0 ? t_length - 1 : 0;

    start = { t_start, grid.ys, grid.xs, 0 };
    count = {       1, grid.ym, grid.xm, nlevels };
  }

  if (units != output_units) {
    auto data_size = grid.xm * grid.ym * nlevels;

    // create a temporary array, convert to output units, and
    // save
    std::vector<double> tmp(data_size);
    for (unsigned int k = 0; k < data_size; ++k) {
      tmp[k] = input[k];
    }

    // convert units "in place"
    units::Converter(metadata.unit_system(), units, output_units)
        .convert_doubles(tmp.data(), tmp.size());

    write_distributed_array_impl(file_name, variable_name, start, count, tmp.data());
  } else {
    write_distributed_array_impl(file_name, variable_name, start, count, input);
  }
  already_written(file_name, variable_name) = true;
}

void OutputWriter::append(const std::string &file_name) {
  append_impl(file_name);
}

void OutputWriter::sync(const std::string &file_name) {
  sync_impl(file_name);
}

void OutputWriter::close(const std::string &file_name) {
  close_impl(file_name);
}

unsigned int OutputWriter::time_dimension_length(const std::string &file_name) {
  return time_dimension_length_impl(file_name);
}

double OutputWriter::last_time_value(const std::string &file_name) {
  return last_time_value_impl(file_name);
}

} // namespace pism

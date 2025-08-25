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
#include <string>
#include <vector>
#include <cassert>

#include "pism/util/io/IO_Flags.hh"
#include "pism/util/Config.hh"
#include "pism/util/GridInfo.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/OutputWriter.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/error_handling.hh"

namespace pism {

/*!
 * Pre-process variable attributes for writing.
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
    time_name            = config.get_string("time.dimension_name");
    use_internal_units   = config.get_flag("output.use_MKS");
    experiment_id        = config.get_string("output.experiment_id");
    experiment_id_name   = config.get_string("output.experiment_id_dimension");
    experiment_id_length = (int)config.get_number("output.experiment_id_max_length");

    if (not experiment_id.empty()) {
      auto format = config.get_string("output.format");
      if (format == "netcdf3" or format == "pnetcdf") {
        throw RuntimeError::formatted(
            PISM_ERROR_LOCATION,
            "cannot save experiment ID \"%s\" to NetCDF-3 output files ('output.format' == \"%s\")",
            experiment_id.c_str(), format.c_str());
      }
    }
  }

  std::string time_name;
  MPI_Comm comm;
  std::map<std::string, std::map<std::string, std::string> > extra_attributes;
  bool use_internal_units;
  std::map<std::tuple<std::string, std::string>, bool> written_time_independent;
  std::map<std::tuple<std::string, std::string>, bool> written_time_dependent;
  std::map<std::string, grid::DistributedGridInfo> grids;
  std::map<std::string, SpatialVariableMetadata> variables;
  std::string experiment_id_name;
  std::string experiment_id;
  int experiment_id_length;
};

bool &OutputWriter::already_written(const std::string &file_name,
                                    const std::string &variable_name,
                                    bool time_dependent) {
  if (time_dependent) {
    return m_impl->written_time_dependent[{ file_name, variable_name }];
  }

  return m_impl->written_time_independent[{ file_name, variable_name }];
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

const SpatialVariableMetadata&
OutputWriter::spatial_variable_info(const std::string &variable_name) const {
  auto i = m_impl->variables.find(variable_name);
  if (i != m_impl->variables.end()) {
    return i->second;
  }
  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "variable '%s' was not added using add_spatial_variable()",
                                variable_name.c_str());
}

void OutputWriter::define_dimension(const std::string &file_name, const std::string &dimension_name,
                                    unsigned int length) {
  define_dimension_impl(file_name, dimension_name, length);
}

void OutputWriter::define_variable(const std::string &file_name, const VariableMetadata &metadata,
                                   const std::vector<std::string> &dims) {
  define_variable_impl(file_name, format_attributes(metadata), dims);
}

void OutputWriter::add_spatial_variable(const SpatialVariableMetadata &metadata,
                                        const grid::DistributedGridInfo &grid) {
  const auto &name = metadata.get_name();

  if (m_impl->variables.find(name) == m_impl->variables.end()) {
    m_impl->variables.insert({ name, metadata });
  }

  if (m_impl->grids.find(name) == m_impl->grids.end()) {
    m_impl->grids[name] = grid;
  }
}

void OutputWriter::define_spatial_variable(const std::string &file_name,
                                           const std::string &variable_name) {

  // Make a copy of `metadata` so we can modify it:
  auto metadata                = spatial_variable_info(variable_name);
  const auto &name             = metadata.get_name();
  const auto &extra_attributes = m_impl->extra_attributes[file_name];

  if (m_impl->use_internal_units) {
    metadata.output_units(metadata["units"]);
  }

  // add extra attributes such as "grid_mapping" and "coordinates". Variables lat, lon,
  // lat_bnds, and lon_bnds should not have these attributes to support CDO (see issue
  // #384).
  //
  // We check names of x and y dimensions to avoid setting extra attributes for variables
  // that use a different grid (e.g. viscous_bed_displacement written by the Lingle-Clark
  // bed deformation model).
  if (metadata.x().get_name() == "x" and metadata.y().get_name() == "y") {
    if (not member(name, { "lat_bnds", "lon_bnds", "lat", "lon" })) {
      for (const auto &attr : extra_attributes) {
        const auto &attr_name  = attr.first;
        const auto &attr_value = attr.second;
        metadata[attr_name]    = attr_value;
      }
    }
  }

  std::vector<std::string> dims{};

  if (not experiment_id().empty()) {
    define_experiment_id(file_name, metadata.unit_system());
    // add the "experiment_id" dimension to the beginning of the list of dimensions
    dims = { m_impl->experiment_id_name };
  }

  if (not metadata.get_time_independent()) {
    dims.push_back(m_impl->time_name);
  }

  // define dimensions and coordinate variables; assemble the list of dimension names:
  //
  // Note the order: y,x,z
  for (const auto &dimension : { metadata.y(), metadata.x(), metadata.z() }) {

    auto dimension_name = dimension.get_name();

    if (dimension_name.empty()) {
      // var.z().dimension_name() is empty if var is a 2D variable
      continue;
    }

    dims.push_back(dimension_name);
    define_dimension(file_name, dimension_name, dimension.length());
    define_variable(file_name, dimension, { dimension_name });
  }

  assert(dims.size() > 1);

  // define the variable itself:
  define_variable(file_name, metadata, dims);
}

void OutputWriter::define_timeseries_variable(const std::string &file_name,
                                              const VariableMetadata &metadata) {

  std::vector<std::string> dims{};

  if (not experiment_id().empty()) {
    define_experiment_id(file_name, metadata.unit_system());
    // add the "experiment_id" dimension to the beginning of the list of dimensions
    dims = { m_impl->experiment_id_name };
  }

  dims.push_back(m_impl->time_name);

  define_variable(file_name, metadata, dims);
}


void OutputWriter::set_global_attributes(
    const std::string &file_name, const std::map<std::string, std::string> &strings,
    const std::map<std::string, std::vector<double> > &numbers) {
  set_global_attributes_impl(file_name, strings, numbers);
}

void OutputWriter::append_time(const std::string &file_name, double time_seconds) {
  append_time_impl(file_name, time_seconds);
  // mark all time-dependent variables as "not written yet"
  m_impl->written_time_dependent.clear();
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

void OutputWriter::write_text(const std::string &file_name, const std::string &variable_name,
                               const std::vector<unsigned int> &start,
                               const std::vector<unsigned int> &count,
                              const std::string &input) {
  write_text_impl(file_name, variable_name, start, count, input);
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

  write_array_impl(file_name, metadata.get_name(), start, count, data.data());
}

void OutputWriter::write_spatial_variable(const std::string &file_name,
                                          const std::string &variable_name,
                                          const double *input) {
  const auto &metadata = spatial_variable_info(variable_name);
  const auto &grid          = grid_info(variable_name);

  // check if we need to write this variable
  bool time_dependent = not metadata.get_time_independent();

  // Avoid writing time-independent variables more than once (saves time when writing to
  // extra_files) and also avoid writing time-dependent variables more than once per time
  // record
  if (already_written(file_name, variable_name, time_dependent)) {
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

      if (not already_written(file_name, dimension_name, false)) {
        write_array(file_name, dimension_name, { 0 }, { (unsigned int)coordinates.size() },
                    coordinates);
        already_written(file_name, dimension_name, false) = true;
      }
    }
  }

  // write experiment ID
  if (not experiment_id().empty()) {
    write_experiment_id(file_name);
  }

  // make sure we have at least one level
  unsigned int n_levels = std::max(metadata.levels().size(), (std::size_t)1);

  std::string units = metadata["units"];
  std::string output_units = m_impl->use_internal_units ? units : metadata["output_units"];

  if (units != output_units) {
    auto data_size = grid.xm * grid.ym * n_levels;

    // create a temporary array, convert to output units, and
    // save
    std::vector<double> tmp(data_size);
    for (unsigned int k = 0; k < data_size; ++k) {
      tmp[k] = input[k];
    }

    // convert units "in place"
    units::Converter(metadata.unit_system(), units, output_units)
        .convert_doubles(tmp.data(), tmp.size());

    write_spatial_variable_impl(file_name, variable_name, tmp.data());
  } else {
    write_spatial_variable_impl(file_name, variable_name, input);
  }
  already_written(file_name, variable_name, time_dependent) = true;
}

void OutputWriter::write_timeseries_variable(const std::string &file_name,
                                             const VariableMetadata &metadata,
                                             const std::vector<unsigned int> &start,
                                             const std::vector<unsigned int> &count,
                                             const std::vector<double> &input) {
  auto S = start;
  auto C = count;

  if (not experiment_id().empty()) {
    write_experiment_id(file_name);
    S.insert(S.cbegin(), 0);
    C.insert(C.cbegin(), 1);
  }

  write_array(file_name, metadata, S, C, input);
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

void OutputWriter::define_experiment_id(const std::string &file_name,
                                        std::shared_ptr<units::System> unit_system) {
  auto &dim_name = m_impl->experiment_id_name;

  VariableMetadata exp_id(dim_name, unit_system);
  // NOTE: this long name is significant: we use it to recognize the experiment ID
  // dimension in File::dimension_type(). This is needed to compute start and count arrays
  // correctly when re-starting from a file containing this dimension.
  exp_id.set_output_type(io::PISM_CHAR).long_name("experiment ID");

  define_dimension(file_name, dim_name, 1);
  define_dimension(file_name, "nc", m_impl->experiment_id_length);
  define_variable(file_name, exp_id, { dim_name, "nc" });
}

void OutputWriter::write_experiment_id(const std::string &file_name) {
  const auto &variable_name = m_impl->experiment_id_name;
  if (already_written(file_name, variable_name, false)) {
    return;
  }

  const auto &exp_id = experiment_id();
  write_text(file_name, variable_name, { 0, 0 }, { 1, (unsigned)exp_id.size() + 1 }, exp_id);

  // mark experiment ID as "already written"
  already_written(file_name, m_impl->experiment_id_name, false) = true;
}

const std::string &OutputWriter::experiment_id() const {
  return m_impl->experiment_id;
}


} // namespace pism

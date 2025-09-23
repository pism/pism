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

#include "io_helpers.hh"
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
static VariableAttributes format_attributes(const VariableAttributes &metadata) {

  // make a copy so we can edit metadata
  auto variable = metadata;

  auto get = [&variable](const std::string &name) {
    auto j = variable.strings.find(name);
    if (j != variable.strings.end()) {
      return j->second;
    }

    return std::string{};
  };

  std::string units        = get("units");
  std::string output_units = get("output_units");

  // output units should never be written to a file
  variable.strings["output_units"] = "";

  bool use_output_units =
      (not units.empty() and not output_units.empty() and units != output_units);

  if (use_output_units) {
    // Replace "units" with "output_units"
    if (variable.is_set("units")) {
      variable.strings["units"] = output_units;
    }

    units::Converter c(variable.unit_system, units, output_units);

    // We need to convert units of valid_min, valid_max and valid_range:
    {
      if (variable.is_set("valid_range")) {
        auto bounds = variable.numbers["valid_range"];

        variable.numbers["valid_range"] = { c(bounds[0]), c(bounds[1]) };
      } else {
        if (variable.is_set("valid_min")) {
          auto min                      = variable.numbers["valid_min"][0];
          variable.numbers["valid_min"] = { c(min) };
        }
        if (variable.is_set("valid_max")) {
          auto max                      = variable.numbers["valid_max"][0];
          variable.numbers["valid_max"] = { c(max) };
        }
      }

      if (variable.is_set("_FillValue")) {
        auto fill                      = variable.numbers["_FillValue"][0];
        variable.numbers["_FillValue"] = { c(fill) };
      }
    }
  }
  return variable;
}

struct OutputWriter::Impl {
  Impl(MPI_Comm comm_, const Config &config)
      : comm(comm_) {
    time_name            = config.get_string("time.dimension_name");
    experiment_id        = config.get_string("output.experiment_id");
    experiment_id_name   = config.get_string("output.experiment_id_dimension");

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

  std::map<std::tuple<std::string, std::string>, bool> written_time_independent;
  std::map<std::tuple<std::string, std::string>, bool> written_time_dependent;
  std::map<std::string, VariableMetadata> variables;
  std::string experiment_id_name;
  std::string experiment_id;
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

OutputWriter::~OutputWriter() {
  delete m_impl;
}

MPI_Comm OutputWriter::comm() const {
  return m_impl->comm;
}

const std::string &OutputWriter::time_name() const {
  return m_impl->time_name;
}

const VariableMetadata &OutputWriter::variable_info(const std::string &variable_name) const {
  auto i = m_impl->variables.find(variable_name);
  if (i != m_impl->variables.end()) {
    return i->second;
  }
  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "variable '%s' was not added using add_variable()",
                                variable_name.c_str());
}

void OutputWriter::define_dimension(const std::string &file_name, const std::string &dimension_name,
                                    unsigned int length) {
  define_dimension_impl(file_name, dimension_name, length);
}

void OutputWriter::define_variable(const std::string &file_name, const std::string &variable_name,
                                   const std::vector<std::string> &dims, io::Type type,
                                   const VariableAttributes &attributes) {
  define_variable_impl(file_name, variable_name, dims, type, format_attributes(attributes));
}

void OutputWriter::define_variable(const std::string &file_name, const VariableMetadata &variable) {

  add_variable(variable);

  // define dimensions and corresponding coordinate variables
  for (const auto &dimension : variable.dimensions()) {
    define_dimension(file_name, dimension.get_name(), dimension.length());
    if (dimension.coordinate_variable()) {
      define_variable(file_name, dimension.get_name(), dimension.dimension_names(),
                      dimension.get_output_type(), dimension.attributes());
    }
  }

  // build the list of dimension names, adding the name of the time dimension (for
  // time-dependent variables) and the experiment ID dimension (if requested)
  auto dimensions = variable.dimension_names();
  {
    if (variable.get_time_dependent()) {
      dimensions.insert(dimensions.begin(), m_impl->time_name);
    }

    if (not m_impl->experiment_id.empty()) {
      dimensions.insert(dimensions.begin(), m_impl->experiment_id_name);
    }
  }

  // define the variable
  define_variable(file_name, variable.get_name(), dimensions, variable.get_output_type(),
                  variable.attributes());
}

void OutputWriter::add_variable(const VariableMetadata &metadata) {
  const auto &name = metadata.get_name();

  if (m_impl->variables.find(name) == m_impl->variables.end()) {
    m_impl->variables.insert({ name, metadata });
  }
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
  const auto &metadata = variable_info(variable_name);

  if (metadata.grid_info() == nullptr) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "variable '%s' has no grid info",
                                  variable_name.c_str());
  }

  const auto &grid = *metadata.grid_info();

  // check if we need to write this variable
  bool time_dependent = metadata.get_time_dependent();

  // Avoid writing time-independent variables more than once (saves time when writing to
  // extra_files) and also avoid writing time-dependent variables more than once per time
  // record
  if (already_written(file_name, variable_name, time_dependent)) {
    return;
  }

  // write dimensions:
  for (const auto &dim : metadata.dimensions()) {
    auto axis_type = axis_type_from_string(dim["axis"]);

    const std::vector<double> *coordinates = nullptr;
    switch (axis_type) {
    case X_AXIS:
      coordinates = &grid.x;
      break;
    case Y_AXIS:
      coordinates = &grid.y;
      break;
    default:
      coordinates = metadata.levels();
    }

    if (coordinates == nullptr) {
      continue;
    }

    auto dimension_name = dim.get_name();
    if (not already_written(file_name, dimension_name, false)) {
      write_array(file_name, dimension_name, { 0 }, { (unsigned int)coordinates->size() },
                  *coordinates);
      already_written(file_name, dimension_name, false) = true;
    }
  }

  // write experiment ID
  if (not experiment_id().empty()) {
    write_experiment_id(file_name);
  }

  // make sure we have at least one level
  size_t n_levels{1};
  if (metadata.levels() != nullptr) {
    n_levels = std::max(metadata.levels()->size(), (std::size_t)1);
  }

  std::string units = metadata["units"];
  std::string output_units = metadata["output_units"];

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

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

#include <memory>
#include <mpi.h>
#include <map>

#include "IO_Flags.hh"
#include "pism/util/io/File.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/OutputWriter.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/error_handling.hh"

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
  Impl(MPI_Comm comm_, const Config &config, const grid::DistributedGridInfo &grid_,
       const VariableMetadata &mapping_)
      : comm(comm_), grid(grid_), mapping(mapping_) {
    time_name          = config.get_string("time.dimension_name");
    use_internal_units = config.get_flag("output.use_MKS");
  }

  const File &get_file(const std::string &filename) {
    if (files[filename] == nullptr) {
      files[filename] =
          std::make_shared<File>(comm, filename, io::PISM_GUESS, io::PISM_READWRITE_MOVE);
    }

    return *files[filename];
  }

  std::string time_name;
  MPI_Comm comm;
  grid::DistributedGridInfo grid;
  VariableMetadata mapping;
  std::map<std::string, std::shared_ptr<File> > files;
  bool use_internal_units;
};

OutputWriter::OutputWriter(MPI_Comm comm, const Config &config,
                           const grid::DistributedGridInfo &grid,
                           const VariableMetadata &cf_mapping)
    : m_impl(new Impl(comm, config, grid, cf_mapping)) {
}

OutputWriter::~OutputWriter() {
  delete m_impl;
}

void OutputWriter::define_dimension(const std::string &filename, const std::string &name,
                                    size_t length) {
  define_dimension_impl(filename, name, length);
}

void OutputWriter::define_variable(const std::string &filename, const VariableMetadata &metadata,
                                   const std::vector<std::string> &dims) {
  define_variable_impl(filename, metadata, dims);
}

void OutputWriter::define_spatial_variable(const std::string &filename,
                                           const SpatialVariableMetadata &metadata) {

  // Make a copy of `metadata` so we can modify it:
  auto var  = metadata;
  const auto &name = var.get_name();
  const auto &cf_mapping = m_impl->mapping;

  if (m_impl->use_internal_units) {
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
    define_dimension(filename, dimension_name, dimension.length());
    define_variable(filename, dimension, { dimension_name });
  }

  assert(dims.size() > 1);

  // define the variable itself:
  define_variable(filename, var, dims);
}

void OutputWriter::write_attributes(const std::string &filename, const VariableMetadata &metadata) {
  write_attributes_impl(filename, metadata);
}

void OutputWriter::append_time(const std::string &filename, double time_seconds) {
  append_time_impl(filename, time_seconds);
}

void OutputWriter::write_array(const std::string &filename, const std::string &name,
                               unsigned int start, unsigned int M, unsigned int N,
                               const std::vector<double> &data) {
  write_array_impl(filename, name, start, M, N, data);
}

void OutputWriter::write_array(const std::string &filename, const VariableMetadata &metadata,
                               unsigned int start, unsigned int M, unsigned int N,
                               const std::vector<double> &input) {
  // create a copy of "data" to change units
  std::vector<double> data = input;

  units::Converter(metadata.unit_system(), metadata["units"], metadata["output_units"])
      .convert_doubles(data.data(), data.size());

  write_array(filename, metadata.get_name(), start, M, N, data);
}

void OutputWriter::write_spatial_variable(const SpatialVariableMetadata &metadata,
                                          const std::string &filename, const double *input) {
  write_spatial_variable_impl(metadata, filename, input);
}

void OutputWriter::close(const std::string &filename) {
  close_impl(filename);
}

void OutputWriter::define_dimension_impl(const std::string &filename, const std::string &name,
                                         size_t length) {
  const auto &file = m_impl->get_file(filename);

  if (file.dimension_exists(name)) {
    return;
  }

  file.define_dimension(name, length);
}

void OutputWriter::define_variable_impl(const std::string &filename,
                                        const VariableMetadata &metadata,
                                        const std::vector<std::string> &dims) {
  const auto &file = m_impl->get_file(filename);

  if (file.variable_exists(metadata.get_name())) {
    return;
  }

  file.define_variable(metadata.get_name(), metadata.get_output_type(), dims);

  write_attributes(filename, metadata);
}

void OutputWriter::write_attributes_impl(const std::string &filename, const std::string &var_name,
                                         const std::map<std::string, std::string> &strings,
                                         const std::map<std::string, std::vector<double> > &numbers,
                                         io::Type output_type) {
  const auto &file = m_impl->get_file(filename);

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

void OutputWriter::append_time_impl(const std::string &filename, double time_seconds) {
  const auto &file = m_impl->get_file(filename);
  auto name = m_impl->time_name;

  write_array_impl(filename, m_impl->time_name, file.dimension_length(name), 1, 1, { time_seconds });
}

void OutputWriter::write_array_impl(const std::string &filename, const std::string &name,
                                    unsigned int start, unsigned int M, unsigned int N,
                                    const std::vector<double> &data) {
  const auto &file = m_impl->get_file(filename);

  file.write_variable(name, { start, 0 }, { M, N }, data.data());
}

void OutputWriter::write_spatial_variable_impl(const SpatialVariableMetadata &metadata,
                                               const std::string &filename, const double *input) {
  const auto &file = m_impl->get_file(filename);

  const auto &name = metadata.get_name();

  if (not file.variable_exists(name)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' in '%s'.", name.c_str(),
                                  file.name().c_str());
  }

  // check if we need to write this variable
  bool time_independent = metadata.get_time_independent();
  // avoid writing time-independent variables more than once (saves time when writing to
  // extra_files)
  if (time_independent and file.get_variable_was_written(name)) {
    return;
  }

  // write dimensions:
  {
    std::map<std::string, const std::vector<double> &> data = {
      { metadata.x().get_name(), m_impl->grid.x },
      { metadata.y().get_name(), m_impl->grid.y },
      { metadata.z().get_name(), metadata.levels() }
    };
    for (const auto &p : data) {
      const auto &dimension   = p.first;
      const auto &coordinates = p.second;
      bool exists             = file.dimension_exists(dimension);
      bool written            = file.get_variable_was_written(dimension);
      if (exists and not written) {
        write_array(filename, dimension, 0, coordinates.size(), 1, coordinates);
        file.set_variable_was_written(dimension);
      }
    }
  }
  
  // make sure we have at least one level
  unsigned int nlevels = std::max(metadata.levels().size(), (size_t)1);

  std::string units = metadata["units"], output_units = metadata["output_units"];

  // override `output_units` if "output.use_MKS" is set
  if (m_impl->use_internal_units) {
    output_units = units;
  }

  // FIXME: use put_vara_double(...) instead and remove write_darray()
  if (units != output_units) {
    size_t data_size = m_impl->grid.xm * m_impl->grid.ym * nlevels;

    // create a temporary array, convert to output units, and
    // save
    std::vector<double> tmp(data_size);
    for (size_t k = 0; k < data_size; ++k) {
      tmp[k] = input[k];
    }

    units::Converter(metadata.unit_system(), units, output_units)
        .convert_doubles(tmp.data(), tmp.size());

    file.write_distributed_array(name, m_impl->grid, nlevels, not time_independent, tmp.data());
  } else {
    file.write_distributed_array(name, m_impl->grid, nlevels, not time_independent, input);
  }
  file.set_variable_was_written(name);
}

void OutputWriter::close_impl(const std::string &filename) {
  m_impl->files[filename].reset();
}

} // namespace pism

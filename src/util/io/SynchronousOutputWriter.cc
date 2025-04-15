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

#include <string>

#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"

namespace pism {

struct Impl;

const File &SynchronousOutputWriter::get_file(const std::string &filename) {
  if (m_files[filename] == nullptr) {
    m_files[filename] =
      std::make_shared<File>(comm(), filename, io::PISM_GUESS, io::PISM_READWRITE_MOVE);
  }

  return *m_files[filename];
}

SynchronousOutputWriter::SynchronousOutputWriter(MPI_Comm comm, const Config &config,
                                                 const grid::DistributedGridInfo &grid,
                                                 const VariableMetadata &mapping)
    : OutputWriter(comm, config, grid, mapping) {
}

void SynchronousOutputWriter::define_dimension_impl(const std::string &filename, const std::string &name,
                                         size_t length) {
  const auto &file = get_file(filename);

  if (file.dimension_exists(name)) {
    return;
  }

  file.define_dimension(name, length);
}

void SynchronousOutputWriter::define_variable_impl(const std::string &filename,
                                        const VariableMetadata &metadata,
                                        const std::vector<std::string> &dims) {
  const auto &file = get_file(filename);

  if (file.variable_exists(metadata.get_name())) {
    return;
  }

  file.define_variable(metadata.get_name(), metadata.get_output_type(), dims);

  write_attributes(filename, metadata);
}

void SynchronousOutputWriter::write_attributes_impl(const std::string &filename, const std::string &var_name,
                                         const std::map<std::string, std::string> &strings,
                                         const std::map<std::string, std::vector<double> > &numbers,
                                         io::Type output_type) {
  const auto &file = get_file(filename);

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

void SynchronousOutputWriter::append_time_impl(const std::string &filename, double time_seconds) {
  const auto &file = get_file(filename);
  auto name = m_impl->time_name;

  write_array_impl(filename, m_impl->time_name, file.dimension_length(name), 1, 1, { time_seconds });
}

void SynchronousOutputWriter::write_array_impl(const std::string &filename, const std::string &name,
                                    unsigned int start, unsigned int M, unsigned int N,
                                    const std::vector<double> &data) {
  const auto &file = get_file(filename);

  file.write_variable(name, { start, 0 }, { M, N }, data.data());
}

void SynchronousOutputWriter::write_spatial_variable_impl(const SpatialVariableMetadata &metadata,
                                               const std::string &filename, const double *input) {
}

void SynchronousOutputWriter::append_impl(const std::string &filename) {
}

void SynchronousOutputWriter::close_impl(const std::string &filename) {
  m_impl->files[filename].reset();
}

} // namespace pism

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

#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/Config.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/GridInfo.hh"

namespace pism {

const File &SynchronousOutputWriter::file(const std::string &file_name) {
  if (m_files[file_name] == nullptr) {
    auto file = std::make_shared<File>(comm(), file_name, m_backend, io::PISM_READWRITE_MOVE);

    file->set_compression_level(m_compression_level);

    m_files[file_name] = file;
  }

  return *m_files[file_name];
}

SynchronousOutputWriter::SynchronousOutputWriter(MPI_Comm comm, const Config &config)
    : OutputWriter(comm, config) {
  m_compression_level = static_cast<int>(config.get_number("output.compression_level"));
  m_backend = string_to_backend(config.get_string("output.format"));
}

void SynchronousOutputWriter::define_variable_impl(const std::string &file_name,
                                                   const VariableMetadata &metadata,
                                                   const std::vector<std::string> &dims) {
  const auto &output_file = file(file_name);

  if (output_file.variable_exists(metadata.get_name())) {
    return;
  }

  const auto &variable_name = metadata.get_name();

  auto type = metadata.get_output_type();

  output_file.define_variable(variable_name, type, dims);

  write_attributes(file_name, variable_name, metadata.all_strings(), metadata.all_doubles(), type);
}

void SynchronousOutputWriter::append_time_impl(const std::string &file_name, double time_seconds) {
  write_array(file_name, time_name(), { time_dimension_length(file_name) }, { 1 },
              { time_seconds });
}

void SynchronousOutputWriter::append_history_impl(const std::string &file_name,
                                                  const std::string &text) {
  const auto &output_file = file(file_name);

  auto old_history = output_file.read_text_attribute("PISM_GLOBAL", "history");

  output_file.write_attribute("PISM_GLOBAL", "history", old_history + text);
}

void SynchronousOutputWriter::append_impl(const std::string &file_name) {
  if (m_files[file_name] == nullptr) {
    m_files[file_name] =
        std::make_shared<File>(comm(), file_name, io::PISM_GUESS, io::PISM_READWRITE);
  }
}

void SynchronousOutputWriter::sync_impl(const std::string &file_name) {
  m_files[file_name]->sync();
}

void SynchronousOutputWriter::close_impl(const std::string &file_name) {
  if (m_files[file_name] != nullptr) {
    m_files[file_name]->close();
    m_files[file_name].reset();
  }
}

void SynchronousOutputWriter::define_dimension_impl(const std::string &file_name,
                                                    const std::string &name, unsigned int length) {
  const auto &output_file = file(file_name);

  if (output_file.dimension_exists(name)) {
    return;
  }

  output_file.define_dimension(name, length);
}

void SynchronousOutputWriter::set_global_attributes_impl(
    const std::string &file_name, const std::map<std::string, std::string> &strings,
    const std::map<std::string, std::vector<double> > &numbers) {
  write_attributes(file_name, "PISM_GLOBAL", strings, numbers, io::PISM_DOUBLE);
}

void SynchronousOutputWriter::write_attributes(
    const std::string &file_name, const std::string &var_name,
    const std::map<std::string, std::string> &strings,
    const std::map<std::string, std::vector<double> > &numbers, io::Type output_type) {
  const auto &output_file = file(file_name);

  // Write text attributes:
  for (const auto &s : strings) {
    const auto &name  = s.first;
    const auto &value = s.second;

    if (value.empty()) {
      continue;
    }

    output_file.write_attribute(var_name, name, value);
  }

  // Write double attributes:
  for (const auto &d : numbers) {
    const auto &name   = d.first;
    const auto &values = d.second;

    if (values.empty()) {
      continue;
    }

    if (output_type == io::PISM_CHAR) {
      // save attributes of a character variable as "double"
      output_type = io::PISM_DOUBLE;
    }

    output_file.write_attribute(var_name, name, output_type, values);
  }
}

unsigned int SynchronousOutputWriter::time_dimension_length_impl(const std::string &file_name) {
  const auto &output_file = file(file_name);
  return output_file.dimension_length(time_name());
}

double SynchronousOutputWriter::last_time_value_impl(const std::string &file_name) {
  const auto &output_file = file(file_name);

  auto t_length = output_file.dimension_length(time_name());

  if (t_length > 0) {
    double time = 0;
    output_file.read_variable("time", { t_length - 1 }, { 1 }, &time);
    return time;
  }
  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "time dimension in '%s' is absent or has length zero",
                                file_name.c_str());
}

void SynchronousOutputWriter::write_array_impl(const std::string &file_name,
                                               const std::string &variable_name,
                                               const std::vector<unsigned int> &start,
                                               const std::vector<unsigned int> &count,
                                               const double *data) {
  const auto &output_file = file(file_name);

  output_file.write_variable(variable_name, start, count, data);
}

void SynchronousOutputWriter::write_text_impl(const std::string &file_name,
                                              const std::string &variable_name,
                                              const std::vector<unsigned int> &start,
                                              const std::vector<unsigned int> &count,
                                              const std::string &input) {
  const auto &output_file = file(file_name);

  output_file.write_text_variable(variable_name, start, count, input);
}

void SynchronousOutputWriter::write_spatial_variable_impl(const std::string &file_name,
                                                          const SpatialVariableMetadata &metadata,
                                                          const double *data) {

  const auto &output_file = file(file_name);

  const auto &variable_name = metadata.get_name();
  const auto &grid = grid_info(variable_name);
  unsigned int n_levels = std::max(metadata.levels().size(), (std::size_t)1);
  
  std::vector<unsigned int> start = { grid.ys, grid.xs, 0 };
  std::vector<unsigned int> count = { grid.ym, grid.xm, n_levels };

  if (not metadata.get_time_independent()) {
    auto t_length = time_dimension_length(file_name);
    auto t_start  = t_length > 0 ? t_length - 1 : 0;

    start.insert(start.cbegin(), t_start);
    count.insert(count.cbegin(), 1);
  }

  if (not experiment_id().empty()) {
    start.insert(start.cbegin(), 0);
    count.insert(count.cbegin(), 1);
  }
  
  output_file.write_variable(variable_name, start, count, data);
}

} // namespace pism

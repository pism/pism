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

#include "pism/util/io/OutputWriter.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism {

OutputFile::OutputFile(std::shared_ptr<OutputWriter> writer, const std::string &file_name)
    : m_file_name(file_name), m_writer(writer) {
  // empty
}

void OutputFile::add_extra_attributes(const std::map<std::string, std::string> &attributes) const {
  m_writer->add_extra_attributes(m_file_name, attributes);
}

void OutputFile::define_dimension(const std::string &dimension_name, unsigned int length) const {
  m_writer->define_dimension(m_file_name, dimension_name, length);
}

void OutputFile::define_variable(const std::string &variable_name,
                                 const std::vector<std::string> &dims, io::Type type,
                                 const VariableAttributes &attributes) const {
  m_writer->define_variable(m_file_name, variable_name, dims, type, attributes);
}

void OutputFile::define_spatial_variable(const SpatialVariableMetadata &metadata) const {
  const auto *grid_info = metadata.grid_info();
  if (grid_info != nullptr) {
    m_writer->add_spatial_variable(metadata, *grid_info);
  }
  m_writer->define_spatial_variable(m_file_name, metadata.get_name());
}

void OutputFile::define_timeseries_variable(const VariableMetadata &metadata) const {
  m_writer->define_timeseries_variable(m_file_name, metadata);
}

void OutputFile::set_global_attributes(
    const std::map<std::string, std::string> &strings,
    const std::map<std::string, std::vector<double> > &numbers) const {
  m_writer->set_global_attributes(m_file_name, strings, numbers);
}

void OutputFile::append_time(double time_seconds) const {
  m_writer->append_time(m_file_name, time_seconds);
}

void OutputFile::append_history(const std::string &text) const {
  m_writer->append_history(m_file_name, text);
}

void OutputFile::write_array(const std::string &variable_name,
                             const std::vector<unsigned int> &start,
                             const std::vector<unsigned int> &count,
                             const std::vector<double> &input) const {
  m_writer->write_array(m_file_name, variable_name, start, count, input);
}

void OutputFile::write_array(const VariableMetadata &metadata,
                             const std::vector<unsigned int> &start,
                             const std::vector<unsigned int> &count,
                             const std::vector<double> &input) const {
  m_writer->write_array(m_file_name, metadata, start, count, input);
}

void OutputFile::write_spatial_variable(const std::string &variable_name,
                                        const double *input) const {
  m_writer->write_spatial_variable(m_file_name, variable_name, input);
}

void OutputFile::write_timeseries_variable(const VariableMetadata &metadata,
                                           const std::vector<unsigned int> &start,
                                           const std::vector<unsigned int> &count,
                                           const std::vector<double> &input) const {
  m_writer->write_timeseries_variable(m_file_name, metadata, start, count, input);
}

void OutputFile::write_text(const std::string &variable_name,
                            const std::vector<unsigned int> &start,
                            const std::vector<unsigned int> &count, const std::string &input) const {
  m_writer->write_text(m_file_name, variable_name, start, count, input);
}

void OutputFile::append() {
  m_writer->append(m_file_name);
}

void OutputFile::sync() {
  m_writer->sync(m_file_name);
}

void OutputFile::close() {
  m_writer->close(m_file_name);
}

unsigned int OutputFile::time_dimension_length() const {
  return m_writer->time_dimension_length(m_file_name);
}

double OutputFile::last_time_value() const {
  return m_writer->last_time_value(m_file_name);
}

std::string OutputFile::name() const {
  return m_file_name;
}

} // namespace pism

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

#include "pism/util/io/OutputWriter.hh"

namespace pism {

OutputFile::OutputFile(std::shared_ptr<OutputWriter> writer, const std::string &file_name)
    : m_file_name(file_name), m_writer(writer) {
  // empty
}

void OutputFile::define_dimension(const std::string &dimension_name, unsigned int length) {
  m_writer->define_dimension(m_file_name, dimension_name, length);
}

void OutputFile::define_variable(const VariableMetadata &metadata,
                                 const std::vector<std::string> &dims) {
  m_writer->define_variable(m_file_name, metadata, dims);
}

void OutputFile::define_spatial_variable(const SpatialVariableMetadata &metadata,
                                         const grid::DistributedGridInfo &grid) {
  m_writer->define_spatial_variable(m_file_name, metadata, grid);
}

void OutputFile::write_attributes(const VariableMetadata &variable) {
  m_writer->write_attributes(m_file_name, variable);
}

void OutputFile::append_time(double time_seconds) {
  m_writer->append_time(m_file_name, time_seconds);
}

void OutputFile::write_array(const std::string &variable_name,
                             const std::vector<unsigned int> &start,
                             const std::vector<unsigned int> &count,
                             const std::vector<double> &input) {
  m_writer->write_array(m_file_name, variable_name, start, count, input);
}

void OutputFile::write_array(const VariableMetadata &metadata,
                             const std::vector<unsigned int> &start,
                             const std::vector<unsigned int> &count,
                             const std::vector<double> &input) {
  m_writer->write_array(m_file_name, metadata, start, count, input);
}

void OutputFile::write_spatial_variable(const SpatialVariableMetadata &metadata,
                                        const double *input) {
  m_writer->write_spatial_variable(m_file_name, metadata, input);
}

void OutputFile::append() {
  m_writer->append(m_file_name);
}

void OutputFile::close() {
  m_writer->close(m_file_name);
}

unsigned int OutputFile::time_dimension_length() {
  return m_writer->time_dimension_length(m_file_name);
}

std::string OutputFile::name() {
  return m_file_name;
}

} // namespace pism

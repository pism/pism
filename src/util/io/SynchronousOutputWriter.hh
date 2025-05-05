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

#ifndef PISM_SYNCHRONOUSOUTPUTWRITER_H
#define PISM_SYNCHRONOUSOUTPUTWRITER_H

#include "pism/util/io/OutputWriter.hh"

#include <memory>

namespace pism {

class File;

namespace io {
enum Backend : int;
}

class SynchronousOutputWriter : public OutputWriter {
public:
  SynchronousOutputWriter(MPI_Comm comm, const Config &config);
  virtual ~SynchronousOutputWriter() = default;

private:
  std::map<std::string, std::shared_ptr<File> > m_files;
  int m_compression_level;
  io::Backend m_backend;

  const File &file(const std::string &file_name);

  void define_dimension_impl(const std::string &file_name, const std::string &name,
                             unsigned int length);

  void define_variable_impl(const std::string &file_name, const VariableMetadata &metadata,
                            const std::vector<std::string> &dims);

  void set_global_attributes_impl(const std::string &file_name,
                                  const std::map<std::string, std::string> &strings,
                                  const std::map<std::string, std::vector<double> > &numbers);

  void write_attributes(const std::string &file_name, const std::string &var_name,
                        const std::map<std::string, std::string> &strings,
                        const std::map<std::string, std::vector<double> > &numbers,
                        io::Type output_type);

  void append_time_impl(const std::string &file_name, double time_seconds);

  void append_history_impl(const std::string &file_name, const std::string &text);

  unsigned int time_dimension_length_impl(const std::string &file_name);

  double last_time_value_impl(const std::string &file_name);

  void write_array_impl(const std::string &file_name, const std::string &variable_name,
                        const std::vector<unsigned int> &start,
                        const std::vector<unsigned int> &count, const double *data);

  void write_distributed_array_impl(const std::string &file_name, const std::string &variable_name,
                                    const std::vector<unsigned int> &start,
                                    const std::vector<unsigned int> &count, const double *data);

  void append_impl(const std::string &file_name);
  void sync_impl(const std::string &file_name);
  void close_impl(const std::string &file_name);
};

} // namespace pism

#endif /* PISM_SYNCHRONOUSOUTPUTWRITER_H */

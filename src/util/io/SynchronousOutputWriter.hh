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

#include "pism/util/io/File.hh"
#include "pism/util/io/OutputWriter.hh"
#include <memory>

namespace pism {

class SynchronousOutputWriter : public OutputWriter {
public:
  SynchronousOutputWriter(MPI_Comm comm, const Config &config,
                          const grid::DistributedGridInfo &grid,
                          const VariableMetadata &mapping);
  virtual ~SynchronousOutputWriter() = default;

private:
  std::map<std::string, std::shared_ptr<File> > m_files;
  
  void define_dimension_impl(const std::string &filename, const std::string &name, size_t length);

  void define_variable_impl(const std::string &filename, const VariableMetadata &metadata,
                            const std::vector<std::string> &dims);

  void write_attributes_impl(const std::string &filename, const std::string &var_name,
                             const std::map<std::string, std::string> &strings,
                             const std::map<std::string, std::vector<double> > &numbers,
                             io::Type output_type);

  void append_time_impl(const std::string &filename, double time_seconds);

  void write_array_impl(const std::string &filename, const std::string &name, unsigned int start,
                        unsigned int M, unsigned int N, const std::vector<double> &data);

  void write_spatial_variable_impl(const SpatialVariableMetadata &metadata,
                                   const std::string &filename, const double *input);

  void append_impl(const std::string &filename);

  void close_impl(const std::string &filename);

  const File &get_file(const std::string &filename);
};

}

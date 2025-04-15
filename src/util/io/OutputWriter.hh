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

#include <map>
#include <string>
#include <vector>

#include <mpi.h>

namespace pism {

class Config;
class SpatialVariableMetadata;
class VariableMetadata;

namespace grid {
class DistributedGridInfo;
}

namespace io {
enum Type : int;
}

class OutputWriter {
public:
  OutputWriter(MPI_Comm comm, const Config &config, const grid::DistributedGridInfo &grid,
               const VariableMetadata &mapping);
  virtual ~OutputWriter();

  void define_dimension(const std::string &filename, const std::string &name, size_t length);

  void define_variable(const std::string &filename, const VariableMetadata &metadata,
                       const std::vector<std::string> &dims);

  void define_spatial_variable(const std::string &filename,
                               const SpatialVariableMetadata &metadata);

  void write_attributes(const std::string &filename, const VariableMetadata &variable);

  void append_time(const std::string &filename, double time_seconds);

  //! Write an one- or two-dimensional array to an output file
  void write_array(const std::string &filename, const std::string &name, unsigned int start,
                   unsigned int M, unsigned int N, const std::vector<double> &data);

  void write_array(const std::string &filename, const VariableMetadata &metadata,
                   unsigned int start, unsigned int M, unsigned int N,
                   const std::vector<double> &input);

  void write_spatial_variable(const SpatialVariableMetadata &metadata, const std::string &filename,
                              const double *input);

  void append(const std::string &filename);

  void close(const std::string &filename);

protected:
  virtual void define_dimension_impl(const std::string &filename, const std::string &name,
                                     size_t length) = 0;

  virtual void define_variable_impl(const std::string &filename,
                                    const VariableMetadata &metadata,
                                    const std::vector<std::string> &dims) = 0;

  virtual void write_attributes(const std::string &filename, const std::string &var_name,
                                const std::map<std::string, std::string> &strings,
                                const std::map<std::string, std::vector<double> > &numbers,
                                io::Type output_type) = 0;

  virtual void append_time_impl(const std::string &filename, double time_seconds) = 0;

  virtual void write_array(const std::string &filename, const std::string &name,
                           const std::vector<unsigned int> &start,
                           const std::vector<unsigned int> &count,
                           const std::vector<double> &data) = 0;

  virtual void write_distributed_array(const std::string &filename, const std::string &name,
                                       const std::vector<unsigned int> &start,
                                       const std::vector<unsigned int> &count,
                                       const double *data) = 0;

  virtual void append_impl(const std::string &filename) = 0;

  virtual void close_impl(const std::string &filename) = 0;

  MPI_Comm comm() const;

private:
  struct Impl;
  Impl *m_impl;
};

} // namespace pism

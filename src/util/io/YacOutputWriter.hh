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

#ifndef PISM_YACOUTPUTWRITER_H
#define PISM_YACOUTPUTWRITER_H

#include <memory>
#include <mpi.h>

#include "OutputWriter.hh"
#include "pism/util/io/OutputWriter.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/json.hpp"

namespace pism {

class File;

namespace io {
enum Backend : int;
}

/*!
 * Synchronous implementation of OutputWriter.
 */
class YacOutputWriter : public OutputWriter {
public:
  YacOutputWriter(MPI_Comm comm, const Config &config, const Geometry & geometry);
  virtual ~YacOutputWriter() = default;

private:
  std::map<std::string, std::shared_ptr<File> > m_files;
  int m_compression_level;
  io::Backend m_backend;
  MPI_Comm intercomm;
  bool yac_initialized = false;
  bool yac_init_finished = false;
  int x_size;
  int y_size;
  int grid_size;
  const Geometry& m_geometry;
  std::map<std::string, int> field_ids;
  std::map<std::string, int> dim_sizes;
  std::map<std::string, unsigned int> variable_tags;
  nlohmann::json non_spatial_variables_metadata, global_attributes;

  std::vector<MPI_Request> field_reqs;
  unsigned int sent_fields_count = 0;

  //YAC variables
  int grid_id;
  int vertex_points_id;

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

  void write_spatial_variable_impl(const std::string &file_name,
                                   const SpatialVariableMetadata &metadata, const double *data);

  void append_impl(const std::string &file_name);
  void sync_impl(const std::string &file_name);
  void close_impl(const std::string &file_name);
  void initialize_yac();
  void define_yac_field(const VariableMetadata &metadata,
                        const std::vector<std::string> &dims);
  void finalize_yac_initialization();
};

} // namespace pism

#endif /* PISM_YACOUTPUTWRITER_H */

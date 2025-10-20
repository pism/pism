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

enum ServerActions {
    FINISH,
    CREATE_FILE,
    SET_FILE_DIMENSION,
    INIT_YAC_GRID,
    SET_FILE_ATTRIBUTES,
    DEFINE_NON_SPATIAL_VARIABLE,
    SEND_NON_SPATIAL_VARIABLE,
    FINISH_YAC_INITIALIZATION,
    DEFINE_SPATIAL_VARIABLE,
    SEND_SPATIAL_VARIABLE,
    UPDATE_TIME_LENGTH
};

/*!
 * Synchronous implementation of OutputWriter.
 */
class YacOutputWriter : public OutputWriter {
public:
  YacOutputWriter(MPI_Comm comm, const Config &config, const Geometry & geometry);
  ~YacOutputWriter();

private:
  std::map<std::string, std::shared_ptr<File> > m_files;
  int m_compression_level;
  io::Backend m_backend;
  MPI_Comm intercomm;
  bool yac_init_finished = false;
  bool yac_grid_initialized = false;
  int grid_size;
  int my_rank = -1;
  int local_x_size;
  int local_y_size;
  const Geometry& m_geometry;
  std::map<std::string, int> field_ids;
  std::map<std::string, std::map<std::string, int>> dim_sizes;
  std::map<std::string, unsigned int> variable_tags;
  std::vector<std::string> text_field_buffers;
  nlohmann::json non_spatial_variables_metadata;
  std::vector<MPI_Request> mpi_requests; 
  std::vector<double *> array_data;

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

  void write_text_impl(const std::string &file_name, const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       const std::string &input);

  void write_spatial_variable_impl(const std::string &file_name,
                                   const SpatialVariableMetadata &metadata, const double *data);

  void append_impl(const std::string &file_name);
  void sync_impl(const std::string &file_name);
  void close_impl(const std::string &file_name);
  void initialize_yac();
  void define_yac_field(const std::string file_name,
                        const VariableMetadata &metadata,
                        const std::vector<std::string> &dims);
  void initialize_grid();
  void finalize_yac_initialization();
  void server_send_action(int server_action_id, 
                          const std::string &server_action_metadata = "");

  // Utility: Given grid size and patch bounds, return global indices of patch vertices
  static std::vector<int> compute_patch_global_indices(int x_global_size, int x_start, int x_size, int y_start, int y_size);
};

} // namespace pism

#endif /* PISM_YACOUTPUTWRITER_H */

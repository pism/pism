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

// These have to match the actions which are defined on the server
enum ServerActions {
    CREATE_FILE,
    SET_FILE_DIMENSION,
    SET_FILE_ATTRIBUTES,
    START_YAC_INITIALIZATION,
    FINISH_YAC_INITIALIZATION,
    DEFINE_NON_SPATIAL_VARIABLE,
    DEFINE_SPATIAL_VARIABLE,
    SEND_NON_SPATIAL_VARIABLE,
    SEND_SPATIAL_VARIABLE,
    UPDATE_TIME_LENGTH,
    FINISH
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
  MPI_Comm m_intercomm;
  bool m_yac_init_finished = false;
  bool m_yac_grid_initialized = false;

  // This flag defines whether the client should also write the files written by the server
  // It might be useful for debugging purposes
  bool m_suppress_client_file_operations = true;
  int m_grid_size;
  bool m_leader;
  int m_max_collection_size = 0;
  const Geometry& m_geometry;
  std::map<std::string, int> m_field_ids;
  std::map<std::string, std::map<std::string, int>> m_dim_sizes;
  std::map<std::string, unsigned int> m_variable_tags;
  std::vector<std::string> m_text_field_buffers;
  std::map<std::string, std::map<std::string, bool> > m_defined_variable, m_defined_dimension;

  std::vector<MPI_Request> m_mpi_requests; 
  std::vector<double *> m_array_data;
  double *** m_yac_raw_send_array = nullptr;

  std::map<std::string, bool> m_server_allowed_files;

  //YAC variables
  int m_grid_id;
  int m_vertex_points_id;

  // --- Server-related subroutines ---
  void create_intercomm();

  void initialize_yac_grid(const std::string &variable_name);

  void define_yac_field(const std::string &file_name, const std::string &variable_name,
                        const std::vector<std::string> &dims, io::Type type,
                        const VariableAttributes &attributes);
  void end_yac_definitions();
  void server_send_action(int server_action_id,
                          const std::string &server_action_metadata = "");
  void server_ensure_file_exists(const std::string &file_name);

  // Utility: Given grid size and patch bounds, return global indices of patch vertices
  static std::vector<int> compute_patch_global_indices(unsigned int x_global_size,
                                                       unsigned int x_start, unsigned int x_size,
                                                       unsigned int y_start, unsigned int y_size);

  // --- Interface implementation and utilities ---
  const File &file(const std::string &file_name);
  void initialize_impl(const std::set<VariableMetadata> &array_variables);

  void define_dimension_impl(const std::string &file_name, const std::string &name,
                             unsigned int length);

  void define_variable_impl(const std::string &file_name, const std::string &variable_name,
                            const std::vector<std::string> &dims, io::Type type,
                            const VariableAttributes &attributes);

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

  void write_distributed_array_impl(const std::string &file_name, const std::string &variable_name,
                                    const double *data);

  void append_impl(const std::string &file_name);
  void sync_impl(const std::string &file_name);
  void close_impl(const std::string &file_name);
};

} // namespace pism

#endif /* PISM_YACOUTPUTWRITER_H */

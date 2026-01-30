/* Copyright (C) 2025, 2026 PISM Authors
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
#include <string>

#include "OutputWriter.hh"
#include "pism/util/io/OutputWriter.hh"
#include "pism/util/json.hpp"

namespace pism {

class File;

namespace io {
enum Backend : int;
}

// These have to match the actions which are defined on the server
enum ServerActions : int {
    OPEN_FILE = 1,
    CLOSE_FILE = 2,
    DEFINE_DIMENSION = 3,
    SET_FILE_ATTRIBUTES = 4,
    APPEND_HISTORY = 5,
    DEFINE_YAC_GRID = 6,
    DEFINE_YAC_FIELD = 7,
    FINISH_YAC_INITIALIZATION = 8,
    DEFINE_VARIABLE = 9,
    SEND_VARIABLE = 10,
    SEND_GRIDDED_VARIABLE = 11,
    APPEND_TIME = 12,
    FINISH = 13,
    SYNC = 14
};

/*!
 * Synchronous implementation of OutputWriter.
 */
class YacOutputWriter : public OutputWriter {
public:
  YacOutputWriter(MPI_Comm comm, const Config &config);
  ~YacOutputWriter();

private:
  MPI_Comm m_intercomm;

  //! True if the current MPI process is responsible for sending non-gridded data.
  bool m_leader;

  //! Length of the time dimension in a file
  std::map<std::string, int> m_time_length;

  //! last time value in an output file
  std::map<std::string, double> m_last_time;
  
  // Note (and possibly FIXME): these lists of variables and dimensions below are
  // inaccurate for files that are opened for appending. This means that when appending to
  // a file PISM will attempt to define variables *once*. The output server code is
  // responsible for checking if a dimension (variable) exists and ignoring a request to
  // define it.

  //! List of all variables defined in a given file (used to avoid defining a variable
  //! more than once)
  std::map<std::string, std::map<std::string, bool> > m_defined_variable;

  //! List of all dimensions defined in a given file (used to avoid defining a dimension
  //! more than once)
  std::map<std::string, std::map<std::string, bool> > m_defined_dimension;

  // --- YAC Grid information

  //! YAC point set ID corresponding to a grid name
  std::map<std::string, int> m_point_set_id;

  //! YAC field ID corresponding to a particular variable (by name)
  std::map<std::string, int> m_field_ids;

  //! Maximum collection size corresponding to a grid name
  std::map<std::string, int> m_max_collection_size;

  //! Size of the local grid patch corresponding to a grid name
  std::map<std::string, int> m_patch_size;

  // --- Buffers ---

  //! buffers used to send text (write_text_impl())
  std::vector<std::string> m_text_buffers;

  //! buffers used to send arrays of double
  std::vector<double *> m_buffers;

  double *** m_yac_raw_send_array = nullptr;

  // --- MPI Communication

  //! Tags for MPI messages sending non-gridded variable data

  enum TagTreatment : int {CREATE_NEW_TAG, GET_EXISTING_TAG};
  int tag(const std::string &variable_name, TagTreatment flag = GET_EXISTING_TAG);
  std::map<std::string, int> m_variable_tags;

  std::vector<MPI_Request> m_mpi_requests;

  // --- Server-related subroutines ---
  void create_intercomm();

  void define_yac_grid(const VariableMetadata &variable);

  void define_yac_field(const VariableMetadata &variable);

  void end_yac_definitions();

  void send_action(int action_id, const nlohmann::json &metadata);

  // --- Interface implementation ---
  void initialize_impl(const std::set<VariableMetadata> &array_variables);

  void define_dimension_impl(const std::string &file_name, const std::string &name,
                             unsigned int length);

  void define_variable_impl(const std::string &file_name, const std::string &variable_name,
                            const std::vector<std::string> &dims, io::Type type,
                            const VariableAttributes &attributes);

  void set_global_attributes_impl(const std::string &file_name,
                                  const std::map<std::string, std::string> &strings,
                                  const std::map<std::string, std::vector<double> > &numbers);

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

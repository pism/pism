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

#include <cstddef>              // size_t
#include <memory>
#include <cmath>
#include <cstring>
#include <mpi.h>

#include "pism/util/Config.hh"
#include "pism/util/GridInfo.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/YacOutputWriter.hh"
#include "pism/util/json.hpp"

extern "C" {
#include "yac.h"
}

namespace pism {

namespace details {
/*!
 * Store variable attributes from `attributes` in a JSON object `json`.
 */
static void to_json(const VariableAttributes &attributes, nlohmann::json &json) {
  for (const auto &attribute : attributes.strings) {
    json[attribute.first] = attribute.second;
  }

  for (const auto &attribute : attributes.numbers) {
    json[attribute.first] = attribute.second;
  }
}

/*!
 * Convert a PISM data type to a string that describes the same type as a NumPy dtype
 * object (in Python).
 */
static std::string to_python_type(pism::io::Type input) {
  static const std::map<pism::io::Type, std::string> type_map = {
    { io::PISM_BYTE, "i1" }, { io::PISM_CHAR, "S1" },  { io::PISM_SHORT, "i2" },
    { io::PISM_INT, "i4" },  { io::PISM_FLOAT, "f4" }, { io::PISM_DOUBLE, "f8" }
  };

  auto it = type_map.find(input);
  return (it != type_map.end()) ? it->second : "None"; // "None" for NC_NAT
}

/*!
 * Calculate global indices for local points given the start and size of the local patch
 * and the global x size.
 */
static std::vector<int> patch_global_indices(unsigned int x_global_size, unsigned int x_start,
                                             unsigned int x_size, unsigned int y_start,
                                             unsigned int y_size) {
  std::vector<int> indices;
  indices.reserve((size_t)x_size * y_size);
  for (unsigned int j = y_start; j < y_start + y_size; ++j) {
    for (unsigned int i = x_start; i < x_start + x_size; ++i) {
      indices.push_back((int)(j * x_global_size + i));
    }
  }
  return indices;
}

} // namespace details

// Even if we are using YAC, certain forms of interaction with the server cannot 
// be made using only YAC functionalities. Sending actions, non-gridded data 
// and metadata (after definitions) are examples of such. In order to be able 
// to send all the information to the server we need to define an intercommunicator,
// which then allows direct MPI communication between the client and the server.
void YacOutputWriter::create_intercomm() {
  // At this point YAC has already been initialized in the PetscInitializer.cc
  // and on the server side. Both client and server components have already been defined.

  // We get the local component communicator and a global communicator which 
  // contains all client and server processes
  //
  // FIXME: global_comm should be de-allocated using MPI_Comm_free(). The same is
  // true about local_group and global_group (using MPI_Group_free()).
  MPI_Comm global_comm = MPI_COMM_NULL;
  {
    const int nbr_comps               = 2;
    const char *comp_names[nbr_comps] = { "pism", "pism_output_server" };
    yac_cget_comps_comm(comp_names, nbr_comps, &global_comm);
  }

  int global_size = 0;
  MPI_Comm_size(global_comm, &global_size);
  std::vector<int> component_leaders_ranks(global_size);

  // We retrieve the process group information from both local and global communicators
  MPI_Group local_group = MPI_GROUP_NULL, global_group = MPI_GROUP_NULL;
  MPI_Comm_group(comm(), &local_group);
  MPI_Comm_group(global_comm, &global_group);

  // For the creation of the intercommunicator we need to set leaders on both groups.
  // We define process 0 of each component to be the leader of its local group.
  // We then find the corresponding rank of each leader in the global group and
  // exchange this information between the processes. 
  // The intercomm creation is then done using this information. 
  const int local_leader_rank[]  = {0};
  int local_leader_global_rank[] = {-1};
  MPI_Group_translate_ranks(local_group, 1, local_leader_rank, 
                            global_group, local_leader_global_rank);
  MPI_Allgather(local_leader_global_rank, 1, MPI_INT, 
                component_leaders_ranks.data(), 1, MPI_INT, global_comm);
  int remote_leader = component_leaders_ranks.back();

  // FIXME: we need to call MPI_Comm_free() in the destructor.
  int tag = 0;
  MPI_Intercomm_create(comm(), local_leader_rank[0], global_comm, remote_leader, tag,
                       &m_intercomm);
}

// Initializes the YAC grid and sends the geometrical information to 
// the server so that it can also initialize its own grid
void YacOutputWriter::initialize_yac_grid(const std::string &variable_name) {

  // Distributed grid containing the domain decomposition information
  const auto &grid = *variable_info(variable_name).grid_info();

  send_action(START_YAC_INITIALIZATION);

  // Sends the global domain sizes to the server
  if (m_leader) {
    // FIXME: one MPI_Isend() sending 2 integers would be enough here
    int x_size = (int)grid.Mx;
    m_mpi_requests.emplace_back();
    MPI_Isend((void *) &x_size, 1, MPI_INT, 0, 0, m_intercomm, &m_mpi_requests.back());
  
    int y_size = (int)grid.My;
    m_mpi_requests.emplace_back();
    MPI_Isend((void *) &y_size, 1, MPI_INT, 0, 0, m_intercomm, &m_mpi_requests.back());
  }

  // Will hold the amount of points for the local grid patch of this process
  int local_patch_size = (int)(grid.xm * grid.ym);
  m_grid_size = local_patch_size;

  // Gathers on the server the size of the local patch from each process
  MPI_Gather(&local_patch_size, 1, MPI_INT, NULL, 1, MPI_INT, 0, m_intercomm);

  // Generate "fake" longitudes and latitudes for the local patch of the grid. Here
  // longitudes and latitudes range from 0 to 1 (inclusive) for the "global" grid (a local
  // patch covers a part of this).
  //
  // This is sufficient for moving data from PISM to the output server and does not
  // require projection info.
  //
  // FIXME: make it possible to use projection info to use real lon,lat coordinates of
  // grid points. This will be necessary for "on the fly" post-processing in the output
  // server.
  std::vector<double> latitudes(local_patch_size);
  std::vector<double> longitudes(local_patch_size);
  {
    const auto &x = grid.x;
    const auto &y = grid.y;

    double x_min  = x.front();
    double y_min  = y.front();
    double x_span = x.back() - x_min;
    double y_span = y.back() - y_min;
    int it        = 0;
    for (auto p : GridPoints(grid)) {
      const int i = p.i(), j = p.j();
      longitudes[it] = (x[i] - x_min) / x_span;
      latitudes[it]  = (y[j] - y_min) / y_span;
      it++;
    }
  }

  // Translate local point indices to global point indices
  auto patch_global_indices =
    details::patch_global_indices(grid.Mx, grid.xs, grid.xm, grid.ys, grid.ym);

  // Sends the global indices of local points to the server, followed by local latitudes and longitudes
  MPI_Gatherv(patch_global_indices.data(), local_patch_size, MPI_INT, NULL, NULL, NULL, MPI_INT, 0,
              m_intercomm);
  MPI_Gatherv(latitudes.data(), local_patch_size, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0,
              m_intercomm);
  MPI_Gatherv(longitudes.data(), local_patch_size, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0,
              m_intercomm);

  // Defines the YAC grid and points using the local points 
  int cyclic_dims[] = {0, 0};
  int nbr_vertices[] = {(int)grid.xm, (int)grid.ym};
  yac_cdef_grid_curve2d("pism_grid", nbr_vertices, cyclic_dims,
                        longitudes.data(), latitudes.data(), &m_grid_id);
  yac_cdef_points_unstruct(m_grid_id, local_patch_size, YAC_LOCATION_CORNER,
                           longitudes.data(), latitudes.data(), &m_vertex_points_id);

  m_yac_grid_initialized = true;
}

// Subroutine to define a YAC field
void YacOutputWriter::define_yac_field(const std::string &file_name,
                                       const std::string &variable_name,
                                       const std::vector<std::string> &dims, io::Type type,
                                       const VariableAttributes &attributes) {

  const auto &variable = variable_info(variable_name);

  int collection_size = std::max((int)variable.levels().size(), 1);

  // Gather all the field metadata into the json object and send to the output server
  {
    nlohmann::json info;

    details::to_json(attributes, info);
    info["dimensions"]      = dims;
    info["dtype"]           = details::to_python_type(type);
    info["file_name"]       = file_name;
    info["variable_name"]   = variable_name;
    info["timestep"]        = "PT1M";
    info["collection_size"] = collection_size;

    send_action(DEFINE_SPATIAL_VARIABLE, info.dump());
  }

  // This allows a single allocation/deallocation of the yac_raw_send_array.
  //
  // With multiple files writing spatial variables it might become a problem
  // if the actual max collection size will be defined after a previous file
  // already finished its initialization.
  //
  // The check on the yac_raw_send_array nullity is to prevent extra
  // deallocation in the destructor.
  if (collection_size > m_max_collection_size and m_yac_raw_send_array == nullptr) {
    m_max_collection_size = collection_size;
  }

  // If the field has already been defined, return
  //
  // Note that a spatial variable can theoretically be defined on the server
  // for multiple files but only one YAC field will be created
  if (m_field_ids.find(variable_name) != m_field_ids.end()) {
    return;
  }

  int field_id = -1;
  yac_cdef_field(variable_name.c_str(), 1, &m_vertex_points_id, 1, collection_size, "PT1M",
                 YAC_TIME_UNIT_ISO_FORMAT, &field_id);
  m_field_ids[variable_name] = field_id;
}

// This subroutine ends the YAC definitions phase.
// No components, grids or fields can be defined after this.
void YacOutputWriter::end_yac_definitions() {
    send_action(FINISH_YAC_INITIALIZATION);
    yac_cenddef();

    // These arrays are deleted in the destructor
    if (m_yac_raw_send_array == nullptr) {
      m_yac_raw_send_array = new double **[m_max_collection_size];
      for (int c = 0; c < m_max_collection_size; c++) {
        m_yac_raw_send_array[c]    = new double *[1];
        m_yac_raw_send_array[c][0] = new double[m_grid_size];
      }
    }

    m_yac_init_finished = true;
}

// This subroutine sends an action to the server, so it knows what to do next.
// An action is composed of an integer action id and an optional metadata field.
// Actions like FINISH need no metadata, while actions like CREATE_FILE need additional 
// parameters in the metadata. 
// TODO: the metadata should be a json string for all actions
void YacOutputWriter::send_action(int action_id,
                                  const std::string &action_metadata) {
  // Only the leader process needs to send actions to the server
  if (not m_leader) {
    return;
  }

  // FIXME: encode action ID in JSON and then use one MPI_Isend() to send the message. The
  // Python code will then use MPI_Probe() + MPI_Get_count() + MPI_Recv() to get all of
  // the info.
  //
  // First the action id is sent to the server
  m_mpi_requests.emplace_back();
  MPI_Isend((void *)&action_id, 1, MPI_INT, 0, 0, m_intercomm, &m_mpi_requests.back());

  if (action_metadata.empty()) {
    return;
  }

  // If there are metadata fields, then the string length and the contents are sent to the server
  int action_metadata_length = (int)action_metadata.length();
  m_mpi_requests.emplace_back();
  MPI_Isend((void *)&action_metadata_length, 1, MPI_INT, 0, 0, m_intercomm, &m_mpi_requests.back());

  //FIXME: Currently there are runtime errors if this is also done asynchronously
  MPI_Send((void *)action_metadata.data(), action_metadata_length, MPI_CHAR, 0, 0,
           m_intercomm);
}

YacOutputWriter::YacOutputWriter(MPI_Comm comm, const Config &config)
    : OutputWriter(comm, config) {

  {
    int rank = -1;
    MPI_Comm_rank(comm, &rank);
    m_leader = (rank == 0);
  }

  create_intercomm();
}

YacOutputWriter::~YacOutputWriter() {
    send_action(FINISH);

    for (auto &i : m_array_data) {
      delete i;
    }

    for (int c = 0; c < m_max_collection_size; c++) {
      delete m_yac_raw_send_array[c][0];
      delete m_yac_raw_send_array[c];
    }

    delete m_yac_raw_send_array;

    yac_cfinalize();
}

/*!
 * Define all the grids and send grid information to the other side.
 */
void YacOutputWriter::initialize_impl(const std::set<VariableMetadata> &array_variables) {

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "FIXME: not implemented");
}

void YacOutputWriter::define_variable_impl(const std::string &file_name,
                                           const std::string &variable_name,
                                           const std::vector<std::string> &dims, io::Type type,
                                           const VariableAttributes &attributes) {

  // If this variable was already defined for this file, return
  if (m_defined_variable[file_name][variable_name]) {
    return;
  }

  const auto &variable = variable_info(variable_name);

  bool gridded_variable = variable.grid_info() != nullptr;

  if (gridded_variable) {
    // Initialize the YAC grid if it has not yet been initialized
    if (not m_yac_grid_initialized) {
      initialize_yac_grid(variable_name);
    }
    // If the variable is gridded, define a yac field for it
    define_yac_field(file_name, variable_name, dims, type, attributes);
  } else {
    // If the variable is not gridded, pack all its metadata and tell the server to
    // define it for this file

    nlohmann::json metadata;
    details::to_json(attributes, metadata);

    metadata["variable_name"] = variable_name;
    metadata["dimensions"]    = dims;
    metadata["dtype"]         = details::to_python_type(type);
    metadata["file_name"]     = file_name;

    if (not dims.empty()) {
      metadata["tag"]                = m_variable_tags.size();
      m_variable_tags[variable_name] = m_variable_tags.size();
    }
    send_action(DEFINE_NON_SPATIAL_VARIABLE, metadata.dump());
  }

  // Save the variable as already defined for this file
  m_defined_variable[file_name][variable_name] = true;
}

void YacOutputWriter::append_time_impl(const std::string &file_name, double time_seconds) {
  // FIXME: we need to increment the value of the time length for this file so
  // time_dimension_length() returns the correct value (and the same with
  // last_time_value()).
  //
  // FIXME: this (the code block below) should not be necessary. Writing to the variable
  // `time_name()` will automatically update the length of the time dimension.
  {
    // Gathers time_dimension_length metadata and sends it to the server
    nlohmann::json file_metadata;
    file_metadata["file_name"] = file_name;
    file_metadata["time_dimension_length"] = time_dimension_length(file_name);
    send_action(UPDATE_TIME_LENGTH, file_metadata.dump());
  }
  
  write_array(file_name, time_name(), { time_dimension_length(file_name) }, { 1 },
              { time_seconds });
}

void YacOutputWriter::append_history_impl(const std::string &file_name, const std::string &text) {
  // FIXME: create the APPEND_HISTORY action and re-implement this

  // Gathers history metadata and sends it to the server
  nlohmann::json info;
  info["file_name"]             = file_name;
  info["attributes"]["history"] = text;
  send_action(SET_FILE_ATTRIBUTES, info.dump());
}

void YacOutputWriter::append_impl(const std::string &file_name) {
  // FIXME: send the "action" to open the file in "append" mode
  //
  // The output server will open the file and send back the time dimension length and the
  // last value of the "time" variable.
  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "FIXME: not implemented");
}

void YacOutputWriter::sync_impl(const std::string &file_name) {
  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "FIXME: not implemented");
}

void YacOutputWriter::close_impl(const std::string &file_name) {
  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "FIXME: not implemented");
}

void YacOutputWriter::define_dimension_impl(const std::string &file_name,
                                            const std::string &name, unsigned int length) {

  // If this dimension has already been defined for this file, return
  if (m_defined_dimension[file_name][name]) {
    return;
  }

  // save the length of the time dimension
  if (name == time_name()) {
    m_time_length[file_name] = (int)length;
  }

  // Gathers the dimension metadata and sends it to the server
  {
    nlohmann::json info;
    info["file_name"]        = file_name;
    info["dimension_name"]   = name;
    info["dimension_length"] = length;
    send_action(SET_FILE_DIMENSION, info.dump());
  }

  m_defined_dimension[file_name][name] = true;
}

void YacOutputWriter::set_global_attributes_impl(
    const std::string &file_name, const std::map<std::string, std::string> &strings,
    const std::map<std::string, std::vector<double> > &numbers) {

  // Gather the global_attributes into the JSON object and send it to the server.
  nlohmann::json attributes_json;
  for (const auto &attribute : strings) {
    attributes_json[attribute.first] = attribute.second;
  }

  for (const auto &attribute : numbers) {
    attributes_json[attribute.first] = attribute.second;
  }

  nlohmann::json file_attributes_json;
  file_attributes_json["file_name"]  = file_name;
  file_attributes_json["attributes"] = attributes_json;
  send_action(SET_FILE_ATTRIBUTES, file_attributes_json.dump());
}

unsigned int YacOutputWriter::time_dimension_length_impl(const std::string &file_name) {
  return m_time_length[file_name];
}

double YacOutputWriter::last_time_value_impl(const std::string &file_name) {
  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "FIXME: not implemented");
  return 0.0;
}

void YacOutputWriter::write_array_impl(const std::string &file_name,
                                       const std::string &variable_name,
                                       const std::vector<unsigned int> &start,
                                       const std::vector<unsigned int> &count, const double *data) {

  // Gathers the variable name into the json object and sends it
  // to the server for identification of which variable to receive
  nlohmann::json info;
  info["file_name"]     = file_name;
  info["variable_name"] = variable_name;
  send_action(SEND_NON_SPATIAL_VARIABLE, info.dump());

  // Non-gridded variables are sent by the leader process
  if (m_leader) {
    // FIXME: arguments `start` and `count` are used incorrectly

    // Buffers the argument array so that the asynchronous operation can finish after its lifetime
    // These arrays are deleted in the destructor
    m_array_data.push_back(new double[count[0]]);
    memcpy(m_array_data.back(), data + start[0], count[0] * sizeof(double));
    m_mpi_requests.emplace_back();
    MPI_Isend((void *)(m_array_data.back()), (int)count[0], MPI_DOUBLE, 0,
              (int)m_variable_tags[variable_name], m_intercomm, &m_mpi_requests.back());
  }
}

void YacOutputWriter::write_text_impl(const std::string &file_name,
                                      const std::string &variable_name,
                                      const std::vector<unsigned int> &start,
                                      const std::vector<unsigned int> &count,
                                      const std::string &input) {

  // Gathers the variable name into the json object and sends it
  // to the server for identification of which variable to receive
  nlohmann::json info;
  info["file_name"]     = file_name;
  info["variable_name"] = variable_name;
  send_action(SEND_NON_SPATIAL_VARIABLE, info.dump());

  // Text variables are sent by the leader process
  if (m_leader) {
    // FIXME: arguments `start` and `count` are used incorrectly

    // Text fields are buffered so that the asynchronous send can finish after the
    // arguments lifetime. Since it is buffered inside of a vector, the de-allocation
    // happens automatically at the destructor
    m_text_field_buffers.push_back(input);
    m_mpi_requests.emplace_back();
    MPI_Isend((void *)(m_text_field_buffers.back().data() + start[0]), (int)count[0], MPI_CHAR, 0,
              (int)m_variable_tags[variable_name], m_intercomm, &m_mpi_requests.back());
  }
}

void YacOutputWriter::write_distributed_array_impl(const std::string &file_name,
                                                   const std::string &variable_name,
                                                   const double *data) {


  auto variable = variable_info(variable_name);

  const auto *grid = variable.grid_info();

  // If the YAC end of definitions has not yet been done, perfom it now
  if (not m_yac_init_finished) {
    end_yac_definitions();
  }

  // Gather the variable name into the JSON object and send it to the server for
  // identification of which variable to receive.
  {
    nlohmann::json info;
    info["file_name"]     = file_name;
    info["variable_name"] = variable_name;
    send_action(SEND_SPATIAL_VARIABLE, info.dump());
  }

  // Copies the data from the argument array to the yac_raw_send_array
  // YAC will automatically buffer the data it is passed to
  // Since the output interface is only called when the output is done,
  // all calls to yac_cput should result in an actual data exchange
  int local_x_size    = (int)grid->xm;
  int local_y_size    = (int)grid->ym;
  int collection_size = (int)variable.levels().size();
  for (int c = 0; c < collection_size; c++) {
    for (int x = 0; x < local_x_size; x++) {
      for (int y = 0; y < local_y_size; y++) {
        int delta_x = collection_size;
        int delta_y = collection_size * local_x_size;

        int pism_index = y * delta_y + x * delta_x + c;

        int yac_index = x + y * local_x_size;

        m_yac_raw_send_array[c][0][yac_index] = data[pism_index];
      }
    }
  }

  int info, error;
  // TODO: we can add a check to verify that the time is still below the simulation end
  // Since the snapshot output calls are normally equal or smaller than the number of time
  // steps this should work fine nonetheless
  yac_cput(m_field_ids[variable_name], collection_size, m_yac_raw_send_array, &info, &error);
}

} // namespace pism

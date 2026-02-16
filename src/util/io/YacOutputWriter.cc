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

#include <cstddef>              // size_t
#include <memory>
#include <cmath>
#include <cstring>
#include <mpi.h>
#include <vector>

#include "pism/util/Config.hh"
#include "pism/util/GridInfo.hh"
#include "pism/util/Grid.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/YacOutputWriter.hh"
#include "pism/external/nlohmann/json.hpp"

extern "C" {
#include "yac.h"
}

namespace pism {

namespace details {
/*!
 * Store variable attributes from `attributes` in a JSON object `json`.
 */
static void to_json(const VariableAttributes &attributes, nlohmann::json &output) {
  for (const auto &attribute : attributes.strings) {
    output[attribute.first] = attribute.second;
  }

  for (const auto &attribute : attributes.numbers) {
    output[attribute.first] = attribute.second;
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

/*!
 * Return the string identifying the grid used by a variable.
 *
 * Uses names of the first two dimensions. (The third dimension, if present, will
 * correspond to a vertical or some other coordinate.)
 */
static std::string grid_name(const VariableMetadata &variable) {
  const auto &dims = variable.dimensions();
  return dims[0].get_name() + "-" + dims[1].get_name();
}

/*!
 * Free buffers used to cache numeric data sent using MPI_Isend().
 */
static void free_array_buffers(std::vector<double *> &buffers) {
  for (auto *ptr : buffers) {
    delete ptr;
  }
  buffers.clear();
}

/*!
 * Compute coordinates of grid points owned by the current MPI process.
 */
void compute_point_coordinates(const grid::DistributedGridInfo &grid,
                               const std::string & proj_string, std::vector<double> &longitudes,
                               std::vector<double> &latitudes) {
  int patch_size = (int)(grid.xm * grid.ym);

  longitudes.resize(patch_size);
  latitudes.resize(patch_size);

  if (proj_string.empty()) {
    // Generate "fake" longitudes and latitudes for the local patch of the grid. Here
    // longitudes and latitudes range from 0 to 1 (inclusive) for the "global" grid (a local
    // patch covers a part of this).
    //
    // This is sufficient for moving data from PISM to the output server and does not
    // require projection info.
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
  } else {
    // FIXME: make it possible to use projection info to use real lon,lat coordinates of
    // grid points. This will be necessary for "on the fly" post-processing in the output
    // server.
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "output writer: non-trivial projections are not implemented yet");
  }
}

} // namespace details

// Even if we are using YAC, certain forms of interaction with the server cannot
// be made using only YAC functionalities. Sending actions, non-gridded data
// and metadata (after definitions) are examples of such. In order to be able
// to send all the information to the server we need to define an intercommunicator,
// which then allows direct MPI communication between the client and the server.
void YacOutputWriter::create_intercomm() {
  // At this point YAC has been initialized in the pism::initialize() and on the server
  // side. Both client and server components have been defined.

  // We get the local component communicator and a global communicator which
  // contains all client and server processes
  MPI_Comm global_comm = MPI_COMM_NULL;
  {
    const int nbr_comps               = 2;
    const char *comp_names[nbr_comps] = { "pism", "pism_output" };
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

  int tag = 0;
  MPI_Intercomm_create(comm(), local_leader_rank[0], global_comm, remote_leader, tag,
                       &m_intercomm);

  MPI_Group_free(&local_group);
  MPI_Group_free(&global_group);
  MPI_Comm_free(&global_comm);
}

int YacOutputWriter::tag(const std::string &variable_name, TagTreatment flag) {
  auto i = m_variable_tags.find(variable_name);
  if (i != m_variable_tags.end()) {
    return i->second;
  }

  if (flag == CREATE_NEW_TAG) {
    int new_tag = (int)m_variable_tags.size();
    m_variable_tags[variable_name] = new_tag;
    return new_tag;
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Bug: no tag for variable '%s'",
                                variable_name.c_str());
}

// Initializes the YAC grid and sends the geometrical information to
// the server so that it can also initialize its own grid
void YacOutputWriter::define_yac_grid(const VariableMetadata &variable) {

  const auto grid_name = details::grid_name(variable);

  if (m_point_set_id.find(grid_name) != m_point_set_id.end()) {
    // this grid was defined already
    return;
  }

  // Distributed grid containing the domain decomposition information
  const auto &grid = *variable.grid_info();

  {
    nlohmann::json info;
    info["grid_name"] = grid_name;

    send_action(DEFINE_YAC_GRID, info);
  }

  // Sends the global domain sizes to the server
  if (m_leader) {
    int grid_size[2] = {(int)grid.Mx, (int)grid.My};
    m_mpi_requests.emplace_back();
    MPI_Isend((void *) &grid_size, 2, MPI_INT, 0, 0, m_intercomm, &m_mpi_requests.back());
  }

  std::vector<double> latitudes;
  std::vector<double> longitudes;
  details::compute_point_coordinates(grid, "", longitudes, latitudes);

  int local_patch_size = (int)latitudes.size();
  // Gathers on the server the size of the local patch from each process
  MPI_Gather(&local_patch_size, 1, MPI_INT, NULL, 1, MPI_INT, 0, m_intercomm);

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
  int grid_id = -1;
  int point_set_id = -1;
  yac_cdef_grid_curve2d(grid_name.c_str(), nbr_vertices, cyclic_dims, longitudes.data(),
                        latitudes.data(), &grid_id);
  yac_cdef_points_unstruct(grid_id, local_patch_size, YAC_LOCATION_CORNER, longitudes.data(),
                           latitudes.data(), &point_set_id);

  m_point_set_id[grid_name] = point_set_id;
}

// Subroutine to define a YAC field
void YacOutputWriter::define_yac_field(const VariableMetadata &variable) {

  const auto &variable_name = variable.get_name();

  // If the field has been defined already, return
  if (m_field_ids.find(variable_name) != m_field_ids.end()) {
    return;
  }

  int collection_size = std::max((int)variable.levels().size(), 1);

  // define the field
  {
    int point_set_id = m_point_set_id[details::grid_name(variable)];

    int field_id = -1;
    yac_cdef_field(variable_name.c_str(), 1, &point_set_id, 1, collection_size, "PT1M",
                   YAC_TIME_UNIT_ISO_FORMAT, &field_id);

    m_field_ids[variable_name] = field_id;
  }

  // tell the output server to define the field
  {
    nlohmann::json info;

    info["variable_name"]   = variable_name;
    info["timestep"]        = "PT1M";
    info["collection_size"] = collection_size;
    info["grid_name"]       = details::grid_name(variable);

    send_action(DEFINE_YAC_FIELD, info);
  }
}

// This subroutine ends the YAC definitions phase.
// No components, grids or fields can be defined after this.
void YacOutputWriter::end_yac_definitions() {
  send_action(FINISH_YAC_INITIALIZATION, {});
  yac_cenddef();

  m_field_buffer = new double[m_field_buffer_size];
}

void YacOutputWriter::send_action(int action_id,
                                  const nlohmann::json &metadata) {
  // Only the leader process needs to send actions to the server
  if (not m_leader) {
    return;
  }

  nlohmann::json message{};
  message["action"] = action_id;
  message["info"]   = metadata;

  m_text_buffers.push_back(message.dump());
  const auto &message_string = m_text_buffers.back();

  // Send the metadata string to the output server:
  m_mpi_requests.emplace_back();
  MPI_Isend((void *)message_string.data(), (int)message_string.length(), MPI_CHAR, 0, 0,
            m_intercomm, &m_mpi_requests.back());
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
  send_action(FINISH, {});

  try {
    waitall();
  } catch (RuntimeError &e) {
    // ignore failures: we should not let exceptions escape a destructor
  }

  details::free_array_buffers(m_buffers);

  delete m_field_buffer;

  MPI_Comm_free(&m_intercomm);

  yac_cfinalize();
}

/*!
 * Define all the grids and send grid information to the other side.
 */
void YacOutputWriter::initialize_impl(const std::set<VariableMetadata> &array_variables) {

  for (const auto &variable : array_variables) {
    if (variable.grid_info() == nullptr) {
      continue;
    }

    auto grid_name = details::grid_name(variable);

    const auto &grid = *variable.grid_info();

    int collection_size = std::max((int)variable.levels().size(), 1);
    int array_size      = (int)(grid.xm * grid.ym * collection_size);
    m_field_buffer_size = std::max(array_size, m_field_buffer_size);

    // define the grid (if necessary)
    define_yac_grid(variable);

    // define the YAC field
    define_yac_field(variable);
  }

  end_yac_definitions();
}

void YacOutputWriter::define_variable_impl(const std::string &file_name,
                                           const std::string &variable_name,
                                           const std::vector<std::string> &dims, io::Type type,
                                           const VariableAttributes &attributes) {

  // If this variable was already defined for this file, return
  if (m_defined_variable[file_name][variable_name]) {
    return;
  }

  {
    nlohmann::json info, nc_attributes;
    details::to_json(attributes, nc_attributes);
    info["attributes"]    = nc_attributes;
    info["dimensions"]    = dims;
    info["dtype"]         = details::to_python_type(type);
    info["file_name"]     = file_name;
    info["variable_name"] = variable_name;
    info["tag"]           = tag(variable_name, CREATE_NEW_TAG);

    send_action(DEFINE_VARIABLE, info);
  }

  // Save the variable as already defined for this file
  m_defined_variable[file_name][variable_name] = true;
}

void YacOutputWriter::append_time_impl(const std::string &file_name, double time_seconds) {
  // Note: these values are updated *without* communication with the output server.
  {
    m_time_length[file_name] += 1;
    m_last_time[file_name] = time_seconds;
  }

  {
    nlohmann::json info;
    info["file_name"] = file_name;
    info["time"]      = time_seconds;
    send_action(APPEND_TIME, info);
  }
}

void YacOutputWriter::append_history_impl(const std::string &file_name, const std::string &text) {
  nlohmann::json info;
  info["file_name"] = file_name;
  info["history"]   = text;
  send_action(APPEND_HISTORY, info);
}

void YacOutputWriter::append_impl(const std::string &file_name) {

  nlohmann::json info;
  info["file_name"] = file_name;
  send_action(OPEN_FILE, info);

  // get file information from the output server
  //
  // FIXME: The Recv + Bcast pair can be replaced with a single Bcast using m_intercomm.
  if (m_leader) {
    MPI_Status status;
    int time_length = -1;
    MPI_Recv(&time_length, 1, MPI_INT, 0, 0, m_intercomm, &status);
    double last_time = -1;
    MPI_Recv(&last_time, 1, MPI_DOUBLE, 0, 0, m_intercomm, &status);

    m_time_length[file_name] = time_length;
    m_last_time[file_name] = last_time;
  }

  // scatter to other ranks in `comm()`:
  MPI_Bcast(&m_time_length[file_name], 1, MPI_INT, 0, comm());
  MPI_Bcast(&m_last_time[file_name], 1, MPI_DOUBLE, 0, comm());
}

void YacOutputWriter::waitall() {
  int error_code =
      MPI_Waitall((int)m_mpi_requests.size(), m_mpi_requests.data(), MPI_STATUSES_IGNORE);
  if (error_code != MPI_SUCCESS) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Fatal error in MPI_Waitall()");
  }
}

void YacOutputWriter::sync_impl(const std::string & /*file_name*/) {
  send_action(SYNC, {});

  waitall();

  m_mpi_requests.clear();
  m_text_buffers.clear();
  details::free_array_buffers(m_buffers);

  // FIXME: consider adding MPI_Barrier(m_intercomm)
}

void YacOutputWriter::close_impl(const std::string &file_name) {
  nlohmann::json info;
  info["file_name"] = file_name;
  send_action(CLOSE_FILE, info);
}

void YacOutputWriter::define_dimension_impl(const std::string &file_name,
                                            const std::string &name, unsigned int length) {

  // If this dimension has been defined already for this file, return
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
    info["file_name"] = file_name;
    info["name"]      = name;
    info["length"]    = length;
    send_action(DEFINE_DIMENSION, info);
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
  send_action(SET_FILE_ATTRIBUTES, file_attributes_json);
}

unsigned int YacOutputWriter::time_dimension_length_impl(const std::string &file_name) {
  return m_time_length[file_name];
}

double YacOutputWriter::last_time_value_impl(const std::string &file_name) {
  return m_last_time[file_name];
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
  info["start"]         = start;
  info["count"]         = count;
  info["tag"]           = tag(variable_name);
  send_action(SEND_VARIABLE, info);

  // Non-gridded variables are sent by the leader process
  if (m_leader) {
    int data_size = 1;
    for (const auto &c : count) {
      data_size *= (int)c;
    }
    // Buffers the argument array so that the asynchronous operation can finish after its
    // lifetime. These arrays are deleted using details::free_array_buffers() in the
    // destructor and sync_impl().
    double *buffer = new double[data_size];
    m_buffers.push_back(buffer);
    memcpy(buffer, data, data_size * sizeof(double));
    m_mpi_requests.emplace_back();
    MPI_Isend((void *)(buffer), data_size, MPI_DOUBLE, 0, tag(variable_name),
              m_intercomm, &m_mpi_requests.back());
  }
}

void YacOutputWriter::write_text_impl(const std::string &file_name,
                                      const std::string &variable_name,
                                      const std::vector<unsigned int> &start,
                                      const std::vector<unsigned int> &count,
                                      const std::string &input) {

  // info["text"] = true indicates that this is a text variable and data should be
  // received using the MPI_CHAR type instead of MPI_DOUBLE
  nlohmann::json info;
  info["file_name"]     = file_name;
  info["variable_name"] = variable_name;
  info["start"]         = start;
  info["count"]         = count;
  info["tag"]           = tag(variable_name);
  info["text"]          = true;
  send_action(SEND_VARIABLE, info);

  // Text variables are sent by the leader process
  if (m_leader) {
    int data_size = 1;
    for (const auto &c : count) {
      data_size *= (int)c;
    }
    // Text fields are buffered so that the asynchronous send can finish after the
    // arguments lifetime. Since it is buffered inside of a vector, the de-allocation
    // happens automatically at the destructor
    m_text_buffers.push_back(input);
    m_mpi_requests.emplace_back();
    MPI_Isend((void *)(m_text_buffers.back().data()), data_size, MPI_CHAR, 0,
              tag(variable_name), m_intercomm, &m_mpi_requests.back());
  }
}

void YacOutputWriter::write_distributed_array_impl(const std::string &file_name,
                                                   const std::string &variable_name,
                                                   const double *data) {


  auto variable = variable_info(variable_name);

  const auto *grid = variable.grid_info();

  // Gather the variable name into the JSON object and send it to the server for
  // identification of which variable to receive.
  {
    nlohmann::json info;
    info["file_name"]      = file_name;
    info["variable_name"]  = variable_name;
    info["ndims"]          = variable.n_spatial_dimensions();
    info["time_dependent"] = variable.get_time_dependent();
    info["grid_name"]      = details::grid_name(variable);
    send_action(SEND_GRIDDED_VARIABLE, info);
  }

  int x_size          = (int)grid->xm;
  int y_size          = (int)grid->ym;
  int collection_size = std::max((int)variable.levels().size(), 1);

  // Assemble the "send_field" argument for yac_cput():
  //
  // This method comes from examples/toy_dummy/dummy_ocean_c.c in YAC's source code tree.
  std::vector<double**> collection_data(collection_size, nullptr);
  std::vector<double*> point_set_data(collection_size, nullptr);
  for (int j = 0; j < collection_size; ++j) {
    point_set_data[j]  = m_field_buffer + (int)(j * (x_size * y_size));
    collection_data[j] = &(point_set_data[j]);
  }

  // Copy data from the argument array to the send_field.
  //
  // Note: YAC will automatically buffer data that is passed to yac_cput().
  {
    // PISM indexing helpers:
    int delta_x_p = collection_size;
    int delta_y_p = collection_size * x_size;
    // YAC indexing helpers:
    int delta_x_y = 1;
    int delta_y_y = x_size;
    for (int c = 0; c < collection_size; c++) {
      for (int x = 0; x < x_size; x++) {
        for (int y = 0; y < y_size; y++) {
          int pism_index = y * delta_y_p + x * delta_x_p + c;
          int yac_index  = y * delta_y_y + x * delta_x_y;

          collection_data[c][0][yac_index] = data[pism_index];
        }
      }
    }
  }
  // Since the output interface is only called when the output is done, all calls to
  // yac_cput() should result in an actual data exchange.
  //
  // TODO: we can add a check to verify that the time is still below the simulation end
  // Since the snapshot output calls are normally equal or smaller than the number of time
  // steps this should work fine nonetheless
  int info, error;
  yac_cput(m_field_ids[variable_name], collection_size, collection_data.data(), &info, &error);
}


} // namespace pism

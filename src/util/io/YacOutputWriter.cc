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
#include <iostream>
#include <memory>
#include <cmath>
#include <cstring>

#include "pism/util/Config.hh"
#include "pism/util/GridInfo.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/YacOutputWriter.hh"

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

static std::string to_python_type(pism::io::Type input) {
  static const std::map<pism::io::Type, std::string> type_map = {
    { io::PISM_BYTE, "i1" }, { io::PISM_CHAR, "S1" },  { io::PISM_SHORT, "i2" },
    { io::PISM_INT, "i4" },  { io::PISM_FLOAT, "f4" }, { io::PISM_DOUBLE, "f8" }
  };

  auto it = type_map.find(input);
  return (it != type_map.end()) ? it->second : "None"; // "None" for NC_NAT
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
  int nbr_comps = 2;
  const char * comp_names[] = {"pism", "pism_output_server"};
  const int local_leader_rank[]  = {0};
  int local_leader_global_rank[] = {-1};
  int remote_leader;
  int global_size;

  MPI_Comm local_comm, global_comm;
  MPI_Group local_group, global_group;

  // We get the local component communicator and a global communicator which 
  // contains all client and server processes
  yac_cget_comp_comm(1, &local_comm);
  yac_cget_comps_comm(comp_names, nbr_comps, &global_comm);
  MPI_Comm_size(global_comm, &global_size);
  {
    int rank = -1;
    MPI_Comm_rank(local_comm, &rank);
    m_leader = (rank == 0);
  }
  std::vector<int> component_leaders_ranks(global_size);

  // We retrieve the process group information from both local and global communicators
  MPI_Comm_group(local_comm, &local_group);
  MPI_Comm_group(global_comm, &global_group);

  // For the creation of the intercommunicator we need to set leaders on both groups.
  // We define process 0 of each component to be the leader of its local group.
  // We then find the corresponding rank of each leader in the global group and
  // exchange this information between the processes. 
  // The intercomm creation is then done using this information. 
  MPI_Group_translate_ranks(local_group, 1, local_leader_rank, 
                            global_group, local_leader_global_rank);
  MPI_Allgather(local_leader_global_rank, 1, MPI_INT, 
                component_leaders_ranks.data(), 1, MPI_INT, global_comm);
  remote_leader = component_leaders_ranks.back();

  MPI_Intercomm_create(local_comm, 0, global_comm, remote_leader, 0, &m_intercomm);
}

// Initializes the YAC grid and sends the gometrical information to 
// the server so that it can also initialize its own grid
void YacOutputWriter::initialize_yac_grid(const std::string &variable_name) {

  // Distributed grid containing the domain decomposition information
  const auto &grid = *variable_info(variable_name).grid_info();

  server_send_action(START_YAC_INITIALIZATION);

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

  array::AccessScope list
    {
      &m_geometry.latitude,
      &m_geometry.longitude,
    };

  int cyclic_dims[] = {0, 0};
  int nbr_vertices[] = {(int)grid.xm, (int)grid.ym};
  std::vector<double> latitudes(local_patch_size);
  std::vector<double> longitudes(local_patch_size);

  // Converts the process local latitudes and longitudes from degrees to radians
  int it = 0;
  for (auto p : GridPoints(grid)) {
    const int i = p.i(), j = p.j();
    latitudes[it] = (m_geometry.latitude(i, j) * M_PI) / 180.0;
    longitudes[it] = (m_geometry.longitude(i, j) * M_PI) / 180.0;
    it++;
  }

  // Translate local point indices to global point indices
  auto patch_global_indices =
      compute_patch_global_indices(grid.Mx, grid.xs, grid.xm, grid.ys, grid.ym);

  // Sends the global indices of local points to the server, followed by local latitudes and longitudes
  MPI_Gatherv(patch_global_indices.data(), local_patch_size, MPI_INT, NULL, NULL, NULL, MPI_INT, 0,
              m_intercomm);
  MPI_Gatherv(latitudes.data(), local_patch_size, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0,
              m_intercomm);
  MPI_Gatherv(longitudes.data(), local_patch_size, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0,
              m_intercomm);

  // Defines the YAC grid and points using the local points 
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

  int collection_size = 1;
  if (dims.size() > 3) {
    collection_size = m_dim_sizes[file_name][dims[3]];
  }

  // Gathers all the field metadata into the json object
  nlohmann::json field_metadata;
  {
    details::to_json(attributes, field_metadata);

    field_metadata["dimensions"]      = dims;
    field_metadata["dtype"]           = details::to_python_type(type);
    field_metadata["file_name"]       = file_name;
    field_metadata["variable_name"]   = variable_name;
    field_metadata["timestep"]        = "PT1M";
    field_metadata["collection_size"] = collection_size;
  }

    // This allows a single allocation/deallocation of the yac_raw_send_array.
    // With multiple files writing spatial variables it might become a problem 
    // if the actual max collection size will be defined after a previous file 
    // already finished its initialization.
    // The check on the yac_raw_send_array nullity is to prevent extra 
    // deallocation in the destructor.
    if (collection_size > m_max_collection_size and m_yac_raw_send_array == nullptr) {
      m_max_collection_size = collection_size;
    }

    server_send_action(DEFINE_SPATIAL_VARIABLE, field_metadata.dump());

    // If the field has already been defined, return
    // Note that a spatial variable can theoretically be defined on the server 
    // for multiple files but only one YAC field will be created
    if (m_field_ids.find(variable_name) != m_field_ids.end()) {
      return;
    }

    int field_id = -1;
    yac_cdef_field(variable_name.c_str(), 1, &m_vertex_points_id, 1,
                   collection_size, "PT1M", YAC_TIME_UNIT_ISO_FORMAT, &field_id);
    m_field_ids[variable_name] = field_id;
}

// This subroutine ends the YAC definitions phase.
// No components, grids or fields can be defined after this.
void YacOutputWriter::end_yac_definitions() {
    server_send_action(FINISH_YAC_INITIALIZATION);
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
void YacOutputWriter::server_send_action(int server_action_id,
                                         const std::string &server_action_metadata) {
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
  MPI_Isend((void *)&server_action_id, 1, MPI_INT, 0, 0, m_intercomm, &m_mpi_requests.back());

  if (server_action_metadata.empty()) {
    return;
  }

  // If there are metadata fields, then the string length and the contents are sent to the server
  int action_metadata_length = (int)server_action_metadata.length();
  m_mpi_requests.emplace_back();
  MPI_Isend((void *)&action_metadata_length, 1, MPI_INT, 0, 0, m_intercomm, &m_mpi_requests.back());

  //FIXME: Currently there are runtime errors if this is also done asynchronously
  MPI_Send((void *)server_action_metadata.data(), action_metadata_length, MPI_CHAR, 0, 0,
           m_intercomm);
}

// This subroutine provides a similar functionality as the synchronous "file" subroutine but for the server side.
// It should be always called before any other action is performed on a server file.
// It will check whether the file already exists and if not will tell the server to create it.
void YacOutputWriter::server_ensure_file_exists(const std::string &file_name) {
    // The server_allowed_files map is used to define which files can be written by the server.
    // This map is checked for all relevant actions which operate on a file. 
    // Currently, due to the different staged and staggered definitions 
    // in YAC and PISM, only snapshot files can be written.
    // These files are however the ones which will likely provide the largest 
    // benefits in model execution time reduction when written asynchronously.
    // Future server extensions to other allowed file classes should be 
    // straightforward by simply adding new clauses to this if condition.
    if(file_name.find("snapshot") != std::string::npos and !m_server_allowed_files[file_name]) {
      m_server_allowed_files[file_name] = true;
      nlohmann::json info;
      info["file_name"] = file_name;
      server_send_action(CREATE_FILE, info.dump());
    }
}

// This subroutine takes the start and size of the local patch and 
// the global x size to calculate the corresponding global indices for local points.
std::vector<int> YacOutputWriter::compute_patch_global_indices(unsigned int x_global_size,
                                                               unsigned int x_start,
                                                               unsigned int x_size,
                                                               unsigned int y_start,
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

const File &YacOutputWriter::file(const std::string &file_name) {
  if (m_files[file_name] == nullptr) {
    auto file = std::make_shared<File>(comm(), file_name, m_backend, io::PISM_READWRITE_MOVE);

    file->set_compression_level(m_compression_level);

    m_files[file_name] = file;
  }

  return *m_files[file_name];
}

YacOutputWriter::YacOutputWriter(MPI_Comm comm, const Config &config, const Geometry& geometry)
    : OutputWriter(comm, config), m_geometry(geometry) {
  m_compression_level = static_cast<int>(config.get_number("output.compression_level"));
  m_backend = string_to_backend(config.get_string("output.format"));
  create_intercomm();
}

YacOutputWriter::~YacOutputWriter() {
    server_send_action(FINISH);

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

void YacOutputWriter::initialize_impl(const std::set<VariableMetadata> &array_variables) {

}

void YacOutputWriter::define_variable_impl(const std::string &file_name,
                                           const std::string &variable_name,
                                           const std::vector<std::string> &dims, io::Type type,
                                           const VariableAttributes &attributes) {
  server_ensure_file_exists(file_name);

  // Checks whether the file for which this action was called is allowed for the server
  if (m_server_allowed_files[file_name]) {
    // Initialize the YAC grid if it has not yet been initialized 
    if(not m_yac_grid_initialized and dims.size() > 1) {
      initialize_yac_grid(variable_name);
    }

    // If this variable was already defined for this file, return
    if (m_defined_variable[file_name][variable_name]) {
      return;
    }

    // Small heuristic to find if the variable is gridded or not
    int horizontal_dims = 0;
    for (const auto &dim : dims) {
      if (dim == "x" or dim == "y") {
        horizontal_dims++;
      }
    }

    if (horizontal_dims == 2) {
      // If the variable is gridded, define a yac field for it
      define_yac_field(file_name, variable_name, dims, type, attributes);
    } else {
      // If the variable is not gridded, pack all its metadata and tell the server to
      // define it for this file

      nlohmann::json metadata;
      details::to_json(attributes, metadata);

      metadata["variable_name"] = variable_name;
      metadata["dimensions"] = dims;
      metadata["dtype"] = details::to_python_type(type);
      metadata["file_name"] = file_name;

      if(not dims.empty()) {
        metadata["tag"] = m_variable_tags.size();
        m_variable_tags[variable_name] = m_variable_tags.size();
      }
      server_send_action(DEFINE_NON_SPATIAL_VARIABLE, metadata.dump());
    }

    // Save the variable as already defined for this file
    m_defined_variable[file_name][variable_name] = true;

    if (m_suppress_client_file_operations) {
      return;
    }
  }

  const auto &output_file = file(file_name);

  if (output_file.variable_exists(variable_name)) {
    return;
  }

  output_file.define_variable(variable_name, type, dims);

  write_attributes(file_name, variable_name, attributes.strings, attributes.numbers, type);
}

void YacOutputWriter::append_time_impl(const std::string &file_name, double time_seconds) {
  if (m_server_allowed_files[file_name]) {
    // FIXME: the argument `time_seconds` is ignored
    // Gathers time_dimension_length metadata and sends it to the server
    nlohmann::json file_metadata;
    file_metadata["file_name"] = file_name;
    file_metadata["time_dimension_length"] = time_dimension_length(file_name);
    server_send_action(UPDATE_TIME_LENGTH, file_metadata.dump());

    if (m_suppress_client_file_operations) {
      return;
    }
  }
  
  write_array(file_name, time_name(), { time_dimension_length(file_name) }, { 1 },
              { time_seconds });
}

void YacOutputWriter::append_history_impl(const std::string &file_name,
                                                  const std::string &text) {
  server_ensure_file_exists(file_name);
  if (m_server_allowed_files[file_name]) {
    // Gathers history metadata and sends it to the server
    nlohmann::json info;
    info["file_name"] = file_name;
    info["attributes"]["history"] = text;
    server_send_action(SET_FILE_ATTRIBUTES, info.dump());

    if (m_suppress_client_file_operations) {
      return;
    }
  }

  const auto &output_file = file(file_name);

  auto old_history = output_file.read_text_attribute("PISM_GLOBAL", "history");

  output_file.write_attribute("PISM_GLOBAL", "history", old_history + text);
}

void YacOutputWriter::append_impl(const std::string &file_name) {
  if (m_files[file_name] == nullptr) {
    m_files[file_name] =
        std::make_shared<File>(comm(), file_name, io::PISM_GUESS, io::PISM_READWRITE);
  }
}

void YacOutputWriter::sync_impl(const std::string &file_name) {
  m_files[file_name]->sync();
}

void YacOutputWriter::close_impl(const std::string &file_name) {
  if (m_files[file_name] != nullptr) {
    m_files[file_name]->close();
    m_files[file_name].reset();
  }
}

void YacOutputWriter::define_dimension_impl(const std::string &file_name,
                                            const std::string &name, unsigned int length) {
  server_ensure_file_exists(file_name);
  if (m_server_allowed_files[file_name]) {
    m_dim_sizes[file_name][name] = (int)length;

      // If this dimension has already been defined for this file, return
    if (m_defined_dimension[file_name][name]) {
      return;
    }

      // Gathers the dimension metadata and sends it to the server
      nlohmann::json file_dim;
      file_dim["file_name"] = file_name;
      file_dim["dimension_name"] = name;
      file_dim["dimension_length"] = length;
      server_send_action(SET_FILE_DIMENSION, file_dim.dump());

      m_defined_dimension[file_name][name] = true;

      if (m_suppress_client_file_operations) {
        return;
      }
  }
  const auto &output_file = file(file_name);

  if (output_file.dimension_exists(name)) {
    return;
  }

  output_file.define_dimension(name, length);
}

void YacOutputWriter::set_global_attributes_impl(
    const std::string &file_name, const std::map<std::string, std::string> &strings,
    const std::map<std::string, std::vector<double> > &numbers) {

  if (m_server_allowed_files[file_name]) {

    // Gathers the global_attributes into the json object
    // and sends it to the server
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
    server_send_action(SET_FILE_ATTRIBUTES, file_attributes_json.dump());

    if (m_suppress_client_file_operations) {
      return;
    }
  }

  write_attributes(file_name, "PISM_GLOBAL", strings, numbers, io::PISM_DOUBLE);
}

void YacOutputWriter::write_attributes(
    const std::string &file_name, const std::string &var_name,
    const std::map<std::string, std::string> &strings,
    const std::map<std::string, std::vector<double> > &numbers, io::Type output_type) {

  //Since this subroutine is only called through others which already have
  //checks for server files, we don't have to also explicitly call it here

  const auto &output_file = file(file_name);

  // Write text attributes:
  for (const auto &s : strings) {
    const auto &name  = s.first;
    const auto &value = s.second;

    if (value.empty()) {
      continue;
    }

    output_file.write_attribute(var_name, name, value);
  }

  // Write double attributes:
  for (const auto &d : numbers) {
    const auto &name   = d.first;
    const auto &values = d.second;

    if (values.empty()) {
      continue;
    }

    if (output_type == io::PISM_CHAR) {
      // save attributes of a character variable as "double"
      output_type = io::PISM_DOUBLE;
    }

    output_file.write_attribute(var_name, name, output_type, values);
  }
}

unsigned int YacOutputWriter::time_dimension_length_impl(const std::string &file_name) {
  server_ensure_file_exists(file_name);
  if (m_server_allowed_files[file_name] and m_suppress_client_file_operations) {
    return m_dim_sizes[file_name][time_name()];
  }

  const auto &output_file = file(file_name);
  return output_file.dimension_length(time_name());
}

double YacOutputWriter::last_time_value_impl(const std::string &file_name) {
  server_ensure_file_exists(file_name);
  const auto &output_file = file(file_name);

  auto t_length = output_file.dimension_length(time_name());

  if (t_length > 0) {
    double time = 0;
    output_file.read_variable("time", { t_length - 1 }, { 1 }, &time);
    return time;
  }
  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "time dimension in '%s' is absent or has length zero",
                                file_name.c_str());
}

void YacOutputWriter::write_array_impl(const std::string &file_name,
                                               const std::string &variable_name,
                                               const std::vector<unsigned int> &start,
                                               const std::vector<unsigned int> &count,
                                               const double *data) {
  server_ensure_file_exists(file_name);
  if(m_server_allowed_files[file_name]) {

    // Gathers the variable name into the json object and sends it 
    // to the server for identification of which variable to receive
    nlohmann::json info;
    info["file_name"] = file_name;
    info["variable_name"] = variable_name;
    server_send_action(SEND_NON_SPATIAL_VARIABLE, info.dump());

    // Non-gridded variables are sent uniquely by the process with rank 0
    if (m_leader) {
      // Buffers the argument array so that the asynchronous operation can finish after its lifetime
      // These arrays are deleted in the destructor
      m_array_data.push_back(new double[count[0]]);
      memcpy(m_array_data.back(), data + start[0], count[0] * sizeof(double));
      m_mpi_requests.emplace_back();
      MPI_Isend((void *)(m_array_data.back()), (int)count[0], MPI_DOUBLE, 0,
                (int)m_variable_tags[variable_name], m_intercomm, &m_mpi_requests.back());
    }

    if (m_suppress_client_file_operations) {
      return;
    }
  }

  const auto &output_file = file(file_name);
  output_file.write_variable(variable_name, start, count, data);
}

void YacOutputWriter::write_text_impl(const std::string &file_name,
                                      const std::string &variable_name,
                                      const std::vector<unsigned int> &start,
                                      const std::vector<unsigned int> &count,
                                      const std::string &input) {
  server_ensure_file_exists(file_name);
  if(m_server_allowed_files[file_name]) {
    // Gathers the variable name into the json object and sends it 
    // to the server for identification of which variable to receive
    nlohmann::json info;
    info["file_name"] = file_name;
    info["variable_name"] = variable_name;
    server_send_action(SEND_NON_SPATIAL_VARIABLE, info.dump());

    // Text variables are sent uniquely by the process with rank 0
    if (m_leader) {
      // Text fields are buffered so that the asynchronous send can finish after the arguments lifetime
      // Since it is buffered inside of a vector, the deallocation happens automatically at the destructor
      m_text_field_buffers.push_back(input);
      m_mpi_requests.emplace_back();
      MPI_Isend((void *)(m_text_field_buffers.back().data() + start[0]), (int)count[0], MPI_CHAR, 0,
                (int)m_variable_tags[variable_name], m_intercomm, &m_mpi_requests.back());
    }

    if (m_suppress_client_file_operations) {
      return;
    }
  }

  const auto &output_file = file(file_name);
  output_file.write_text_variable(variable_name, start, count, input);
}

void YacOutputWriter::write_distributed_array_impl(const std::string &file_name,
                                                   const std::string &variable_name,
                                                   const double *data) {


  auto variable = variable_info(variable_name);

  const auto *grid = variable.grid_info();
  
  server_ensure_file_exists(file_name);
  if (m_server_allowed_files[file_name]) {
    // If the YAC end of definitions has not yet been done, perfom it now
    if(not m_yac_init_finished) {
      end_yac_definitions();
    }

    // Gathers the variable name into the json object and sends it 
    // to the server for identification of which variable to receive
    {
      nlohmann::json info;
      info["file_name"]     = file_name;
      info["variable_name"] = variable_name;
      server_send_action(SEND_SPATIAL_VARIABLE, info.dump());
    }

    // Copies the data from the argument array to the yac_raw_send_array
    // YAC will automatically buffer the data it is passed to
    // Since the output interface is only called when the output is done,
    // all calls to yac_cput should result in an actual data exchange
    int local_x_size = (int)grid->xm;
    int local_y_size = (int)grid->ym;
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
    if (m_suppress_client_file_operations) {
      return;
    }
  }

  unsigned int n_levels = std::max(variable.levels().size(), (std::size_t)1);
  const auto &output_file = file(file_name);

  std::vector<unsigned int> start = { grid->ys, grid->xs, 0 };
  std::vector<unsigned int> count = { grid->ym, grid->xm, n_levels };

  if (variable.get_time_dependent()) {
    auto t_length = time_dimension_length(file_name);
    auto t_start  = t_length > 0 ? t_length - 1 : 0;

    start.insert(start.cbegin(), t_start);
    count.insert(count.cbegin(), 1);
  }

  if (not experiment_id().empty()) {
    start.insert(start.cbegin(), 0);
    count.insert(count.cbegin(), 1);
  }

  output_file.write_variable(variable_name, start, count, data);
}

} // namespace pism

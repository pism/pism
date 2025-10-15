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

#include <iostream>
#include <memory>
#include <cmath>
#include <cstring>

#include "pism/util/io/YacOutputWriter.hh"
#include "pism/util/io/IO_Flags.hh"
#include "YacOutputWriter.hh"
#include "pism/util/Config.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/GridInfo.hh"
#include "mpi.h"

extern "C" {
#include "yac.h"
}

namespace pism {

std::string pism_type_to_python_nc_type(pism::io::Type input) {
  static const std::map<pism::io::Type, std::string> type_map = {
    {io::PISM_BYTE,  "i1"},
    {io::PISM_CHAR,  "S1"},
    {io::PISM_SHORT, "i2"},
    {io::PISM_INT,   "i4"},
    {io::PISM_FLOAT, "f8"},
    {io::PISM_DOUBLE,"f8"}
  };

  auto it = type_map.find(input);
  return (it != type_map.end()) ? it->second : "None";  // "None" for NC_NAT
}


void YacOutputWriter::initialize_yac() {
  int nbr_comps = 2;

  const char * comp_names[] = {"pism", "pism_output_server"};
  const int local_leader_rank[]  = {0};
  int global_leader_rank[] = {-1};
  int remote_leader;
  int global_size;

  MPI_Comm local_comm, global_comm;
  MPI_Group local_group, global_group;

  yac_cget_comp_comm(1, &local_comm);
  yac_cget_comps_comm(comp_names, nbr_comps, &global_comm);
  MPI_Comm_size(global_comm, &global_size);
  MPI_Comm_rank(local_comm, &local_rank);
  std::vector<int> component_leaders_ranks(global_size);

  MPI_Comm_group(local_comm, &local_group);
  MPI_Comm_group(global_comm, &global_group);

  MPI_Group_translate_ranks(local_group, 1, local_leader_rank, global_group, global_leader_rank);
  MPI_Allgather(global_leader_rank, 1, MPI_INT, component_leaders_ranks.data(), 1, MPI_INT, global_comm);
  remote_leader = component_leaders_ranks.back();

  MPI_Intercomm_create(local_comm, 0, global_comm, remote_leader, 0, &intercomm);
  yac_initialized = true;
}


void YacOutputWriter::initialize_grid() {
  int local_patch_size = -1;
  auto global_grid = m_geometry.latitude.grid();
  auto distributed_grid = m_geometry.latitude.grid()->info();
  x_size = global_grid->Mx();
  y_size = global_grid->My();
  local_x_size = distributed_grid.xm;
  local_y_size = distributed_grid.ym;
  
  int server_action = INIT_YAC_GRID;
  field_reqs.emplace_back();
  MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

  field_reqs.emplace_back();
  MPI_Isend((void *) &x_size, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

  field_reqs.emplace_back();
  MPI_Isend((void *) &y_size, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

  std::vector<int> patch_global_indices = compute_patch_global_indices(
                                                x_size, 
                                                distributed_grid.xs, 
                                                distributed_grid.xm,
                                                distributed_grid.ys, 
                                                distributed_grid.ym);
  local_patch_size = patch_global_indices.size();

  MPI_Gather(&local_patch_size, 1, MPI_INT, NULL, 1, MPI_INT, 0, intercomm);

  array::AccessScope list
    {
      &m_geometry.latitude,
      &m_geometry.longitude,
    };

  grid_size = local_x_size * local_y_size;
  int cyclic_dims[] = {0, 0};
  int nbr_vertices[] = {local_x_size, local_y_size};
  std::vector<double> latitudes(local_patch_size);
  std::vector<double> longitudes(local_patch_size);

  int it = 0;
  for (auto p = global_grid->points(); p; p.next(), it++) {
    const int i = p.i(), j = p.j();
    latitudes[it] = (m_geometry.latitude(i, j) * M_PI) / 180.0;
    longitudes[it] = (m_geometry.longitude(i, j) * M_PI) / 180.0;
  }

  MPI_Gatherv(patch_global_indices.data(), local_patch_size, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, intercomm);
  MPI_Gatherv(latitudes.data(), local_patch_size, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, intercomm);
  MPI_Gatherv(longitudes.data(), local_patch_size, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, intercomm);

  yac_cdef_grid_curve2d("pism_grid", nbr_vertices, cyclic_dims,
                        longitudes.data(), latitudes.data(), &grid_id);

  yac_cdef_points_unstruct(grid_id, grid_size, YAC_LOCATION_CORNER,
                           longitudes.data(), latitudes.data(), &vertex_points_id);

  yac_grid_initialized = true;
}

void YacOutputWriter::define_yac_field(const std::string file_name,
                                       const VariableMetadata &metadata,
                                       const std::vector<std::string> &dims){
    int field_id, collection_size = 1;
    nlohmann::json field_metadata;

    field_metadata["dimensions"] = dims;
    field_metadata["dtype"] = pism_type_to_python_nc_type(metadata.get_output_type());

    for (auto string_attribute : metadata.all_strings())
        field_metadata[string_attribute.first] = string_attribute.second;

    for (auto double_attribute : metadata.all_doubles())
        field_metadata[double_attribute.first] = double_attribute.second;

    if(dims.size() > 3)
        collection_size = dim_sizes[file_name][dims[3]];

    field_metadata["timestep"] = "PT1M";
    field_metadata["collection_size"] = collection_size;
    field_metadata["variable_name"] = metadata.get_name();
    field_metadata["file_name"] = file_name;
    server_define_spatial_variable(file_name, field_metadata.dump());

    if(field_ids.find(metadata.get_name()) != field_ids.end()) return;

    yac_cdef_field(metadata.get_name().c_str(), 1, &vertex_points_id, 1,
                   collection_size, "PT1M", YAC_TIME_UNIT_ISO_FORMAT, &field_id);
    yac_cdef_field_metadata("pism", "pism_grid", metadata.get_name().c_str(), field_metadata.dump().c_str());
    field_ids[metadata.get_name()] = field_id;
}

const File &YacOutputWriter::file(const std::string &file_name) {
  if (m_files[file_name] == nullptr) {
    auto file = std::make_shared<File>(comm(), file_name, m_backend, io::PISM_READWRITE_MOVE);

    file->set_compression_level(m_compression_level);

    m_files[file_name] = file;
    file_time_lengths[file_name] = 1;

    if(file_name.find("snapshot") != std::string::npos or file_name.find("timeseries") != std::string::npos)
      server_create_file(file_name);
  }

  return *m_files[file_name];
}

YacOutputWriter::YacOutputWriter(MPI_Comm comm, const Config &config, const Geometry& geometry)
    : OutputWriter(comm, config), m_geometry(geometry) {
  m_compression_level = static_cast<int>(config.get_number("output.compression_level"));
  m_backend = string_to_backend(config.get_string("output.format"));
  initialize_yac();
}

YacOutputWriter::~YacOutputWriter() {
    int server_action = FINISH;
    field_reqs.emplace_back();
    MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());
}

void YacOutputWriter::server_send_action_metadata(const std::string &action_metadata) {
    int action_metadata_length = action_metadata.length();

    field_reqs.emplace_back();
    MPI_Isend((void *) &action_metadata_length, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

    text_field_buffers.push_back(action_metadata);
    field_reqs.emplace_back();
    MPI_Isend((void *) text_field_buffers.back().data(), action_metadata_length, MPI_CHAR, 0, 0, intercomm, &field_reqs.back());
}

void YacOutputWriter::server_set_file_dimension(const std::string &file_name, 
                                                const std::string &name, 
                                                unsigned int length) {
      nlohmann::json file_dim;
      file_dim["file_name"] = file_name;
      file_dim["dimension_name"] = name;
      file_dim["dimension_length"] = length;

      std::string file_dimensions = file_dim.dump();
      int file_dimensions_length = file_dimensions.length();

      int server_action = SET_FILE_DIMENSION;
      field_reqs.emplace_back();
      MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());
      server_send_action_metadata(file_dimensions);
}

void YacOutputWriter::server_create_file(const std::string &file_name) {
    int server_action = CREATE_FILE;
    field_reqs.emplace_back(); 
    MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());
    
    int file_name_length = file_name.length();
    server_send_action_metadata(file_name);
}

void YacOutputWriter::server_set_file_attributes(const std::string &file_name) {
    int server_action = SET_FILE_ATTRIBUTES;
    field_reqs.emplace_back(); 
    MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

    nlohmann::json file_attributes_json;
    file_attributes_json["file_name"] = file_name;
    file_attributes_json["attributes"] = global_attributes[file_name];

    std::string file_attributes = file_attributes_json.dump();
    int file_attributes_length = file_attributes.length();
    server_send_action_metadata(file_attributes); 
}

void YacOutputWriter::server_define_non_spatial_variable(const std::string &file_name, 
                                                         const std::string &variable_metadata) {
    int server_action = DEFINE_NON_SPATIAL_VARIABLE;
    field_reqs.emplace_back();
    MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

    int variable_metadata_length = variable_metadata.length();
    server_send_action_metadata(variable_metadata); 
}

void YacOutputWriter::server_define_spatial_variable(const std::string &file_name, 
                                                     const std::string &variable_metadata) {
    int server_action = DEFINE_SPATIAL_VARIABLE;
    field_reqs.emplace_back();
    MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

    int variable_metadata_length = variable_metadata.length();
    server_send_action_metadata(variable_metadata); 
}

void YacOutputWriter::define_variable_impl(const std::string &file_name,
                                           const VariableMetadata &metadata,
                                           const std::vector<std::string> &dims) {
  if(file_name.find("snapshot") != std::string::npos and not yac_grid_initialized and dims.size() > 1) {
    initialize_grid();
  }

  const auto &output_file = file(file_name);

  if (output_file.variable_exists(metadata.get_name())) {
    return;
  }

  if(file_name.find("snapshot") != std::string::npos and yac_grid_initialized and
     file_name != current_snapshot_file) {
      current_snapshot_file = file_name;
  }

  if (file_name.find("snapshot") != std::string::npos or
      file_name.find("timeseries") != std::string::npos) {

    int horizontal_dims = 0;
    for (auto dim : dims)
      if (dim == "x" or dim == "y")
        horizontal_dims++;

    if (horizontal_dims == 2) {
        define_yac_field(file_name, metadata, dims);
    } else {
      non_spatial_variables_metadata[metadata.get_name()]["variable_name"] = metadata.get_name();
      non_spatial_variables_metadata[metadata.get_name()]["dimensions"] = dims;
      non_spatial_variables_metadata[metadata.get_name()]["dtype"] = pism_type_to_python_nc_type(metadata.get_output_type());
      non_spatial_variables_metadata[metadata.get_name()]["file_name"] = file_name;

      for (auto string_attribute : metadata.all_strings())
        non_spatial_variables_metadata[metadata.get_name()][string_attribute.first] = string_attribute.second;

      for (auto double_attribute : metadata.all_doubles())
        non_spatial_variables_metadata[metadata.get_name()][double_attribute.first] = double_attribute.second;

      if(dims.size() > 0) {
        non_spatial_variables_metadata[metadata.get_name()]["tag"] = variable_tags.size();
        variable_tags[metadata.get_name()] = variable_tags.size();
      }
      server_define_non_spatial_variable(file_name, non_spatial_variables_metadata[metadata.get_name()].dump());
    }
  }

  const auto &variable_name = metadata.get_name();

  auto type = metadata.get_output_type();

  output_file.define_variable(variable_name, type, dims);

  write_attributes(file_name, variable_name, metadata.all_strings(), metadata.all_doubles(), type);
}

void YacOutputWriter::append_time_impl(const std::string &file_name, double time_seconds) {

  if (file_name.find("snapshot") != std::string::npos or
      file_name.find("timeseries") != std::string::npos) {
        
    int server_action = UPDATE_TIME_LENGTH;
    field_reqs.emplace_back();
    MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

    nlohmann::json file_metadata;
    file_metadata["file_name"] = file_name;
    file_metadata["time_dimension_length"] = time_dimension_length(file_name);

    std::string file_info = file_metadata.dump();
    int file_info_length = file_info.length();
    server_send_action_metadata(file_info);
  }
  
  write_array(file_name, time_name(), { time_dimension_length(file_name) }, { 1 },
              { time_seconds });
}

void YacOutputWriter::append_history_impl(const std::string &file_name,
                                                  const std::string &text) {
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
  const auto &output_file = file(file_name);

  if (output_file.dimension_exists(name)) {
    return;
  }

  dim_sizes[file_name][name] = length;

  if(file_name.find("snapshot") != std::string::npos or file_name.find("timeseries") != std::string::npos)
    server_set_file_dimension(file_name, name, length);

  output_file.define_dimension(name, length);
}

void YacOutputWriter::set_global_attributes_impl(
    const std::string &file_name, const std::map<std::string, std::string> &strings,
    const std::map<std::string, std::vector<double> > &numbers) {

    for (auto string_attribute : strings)
        global_attributes[file_name][string_attribute.first] = string_attribute.second;

    for (auto double_attribute : numbers)
        global_attributes[file_name][double_attribute.first] = double_attribute.second;

    if(file_name.find("snapshot") != std::string::npos or file_name.find("timeseries") != std::string::npos)
      server_set_file_attributes(file_name);

  write_attributes(file_name, "PISM_GLOBAL", strings, numbers, io::PISM_DOUBLE);
}

void YacOutputWriter::write_attributes(
    const std::string &file_name, const std::string &var_name,
    const std::map<std::string, std::string> &strings,
    const std::map<std::string, std::vector<double> > &numbers, io::Type output_type) {
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
  const auto &output_file = file(file_name);
  return output_file.dimension_length(time_name());
}

double YacOutputWriter::last_time_value_impl(const std::string &file_name) {
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

void YacOutputWriter::finalize_yac_initialization() {
    nlohmann::json serialized_dims, component_metadata;

    for (auto file : m_files) {
      for (auto dim : dim_sizes[file.first])
          serialized_dims[dim.first] = dim.second;

      component_metadata[file.first]["global"] = global_attributes[file.first];
      component_metadata[file.first]["non_field_variables"] = non_spatial_variables_metadata;
      component_metadata[file.first]["dimensions"] = serialized_dims;
    }

    int server_action = FINISH_YAC_INITIALIZATION;
    field_reqs.emplace_back();
    MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

    yac_cdef_component_metadata("pism", component_metadata.dump().c_str());
    yac_cdef_grid_metadata("pism_grid", serialized_dims.dump().c_str());
    yac_cenddef();
    yac_init_finished = true;
}

void YacOutputWriter::write_array_impl(const std::string &file_name,
                                               const std::string &variable_name,
                                               const std::vector<unsigned int> &start,
                                               const std::vector<unsigned int> &count,
                                               const double *data) {
  const auto &output_file = file(file_name);
  MPI_Datatype send_type;
  MPI_Request send_req_handle;

  if(not yac_init_finished and file_name.find("snapshot") != std::string::npos) {
    finalize_yac_initialization();
  }

  if(file_name.find("snapshot") != std::string::npos) {
    if( non_spatial_variables_metadata[variable_name]["dtype"] == "f8")
        send_type = MPI_DOUBLE;
    else
        send_type = MPI_INT;

    int server_action = SEND_NON_SPATIAL_VARIABLE;
    field_reqs.emplace_back();
    MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

    nlohmann::json variable_info_json;
    variable_info_json["file_name"] = file_name;
    variable_info_json["variable_name"] = variable_name;
    
    std::string variable_info = variable_info_json.dump();
    int variable_info_length = variable_info.length();
    server_send_action_metadata(variable_info);

    MPI_Isend((void *) (data + start[0]), count[0], send_type, 0, variable_tags[variable_name], intercomm, &send_req_handle);
    field_reqs.push_back(send_req_handle);
    sent_fields_count++;
  }

  output_file.write_variable(variable_name, start, count, data);
}

void YacOutputWriter::write_text_impl(const std::string &file_name,
                                      const std::string &variable_name,
                                      const std::vector<unsigned int> &start,
                                      const std::vector<unsigned int> &count,
                                      const std::string &input) {

  MPI_Request send_req_handle;

  if(file_name.find("snapshot") != std::string::npos) {

    int server_action = SEND_NON_SPATIAL_VARIABLE;
    field_reqs.emplace_back();
    MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

    nlohmann::json variable_info_json;
    variable_info_json["file_name"] = file_name;
    variable_info_json["variable_name"] = variable_name;
    
    std::string variable_info = variable_info_json.dump();
    int variable_info_length = variable_info.length();
    server_send_action_metadata(variable_info);

    text_field_buffers.push_back(input);
    MPI_Isend((void *) (text_field_buffers.back().data() + start[0]), count[0], MPI_CHAR, 0, variable_tags[variable_name], intercomm, &send_req_handle);
    field_reqs.push_back(send_req_handle);
    sent_fields_count++;
  }

  const auto &output_file = file(file_name);
  output_file.write_text_variable(variable_name, start, count, input);
}

void YacOutputWriter::write_spatial_variable_impl(const std::string &file_name,
                                           const SpatialVariableMetadata &metadata,
                                           const double *data) {
  const auto &output_file = file(file_name);

  const auto &variable_name = metadata.get_name();
  const auto &grid = grid_info(variable_name);
  unsigned int n_levels = std::max(metadata.levels().size(), (std::size_t)1);

  if(file_name.find("snapshot") != std::string::npos and time_dimension_length(file_name) != file_time_lengths[file_name]) {
      file_time_lengths[file_name] = time_dimension_length(file_name);
  }

  if(not yac_init_finished and file_name.find("snapshot") != std::string::npos) {
    finalize_yac_initialization();
  }

  if (file_name.find("snapshot") != std::string::npos) {
    int collection_size = metadata.z().length();

    int server_action = SEND_SPATIAL_VARIABLE;
    field_reqs.emplace_back();
    MPI_Isend((void *) &server_action, 1, MPI_INT, 0, 0, intercomm, &field_reqs.back());

    nlohmann::json variable_info_json;
    variable_info_json["file_name"] = file_name;
    variable_info_json["variable_name"] = variable_name;
    
    std::string variable_info = variable_info_json.dump();
    int variable_info_length = variable_info.length();
    server_send_action_metadata(variable_info);

    double ***send_field = new double **[collection_size];
    for (int c = 0; c < collection_size; c++) {
      // FIXME: memory leaks
      send_field[c]    = new double *[1];
      send_field[c][0] = new double[grid_size];
      for (int x = 0; x < local_x_size; x++) {
        for (int y = 0; y < local_y_size; y++) {
          int delta_x = collection_size;
          int delta_y = collection_size * local_x_size;

          int pism_index = y * delta_y + x * delta_x + c;

          int yac_index = x + y * local_x_size;

          send_field[c][0][yac_index] = data[pism_index];
        }
      }
    }

    int info, error;
    yac_cput(field_ids[variable_name], collection_size, send_field, &info, &error);
    written_vars[variable_name] = true;
  }

  std::vector<unsigned int> start = { grid.ys, grid.xs, 0 };
  std::vector<unsigned int> count = { grid.ym, grid.xm, n_levels };

  if (not metadata.get_time_independent()) {
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

std::vector<int> YacOutputWriter::compute_patch_global_indices(int x_global_size, int x_start, int x_size, int y_start, int y_size) {
    std::vector<int> indices;
    indices.reserve(x_size * y_size);
    for (int j = y_start; j < y_start + y_size; ++j) {
        for (int i = x_start; i < x_start + x_size; ++i) {
            indices.push_back(j * x_global_size + i);
        }
    }
    return indices;
}

} // namespace pism

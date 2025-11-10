# Copyright (C) 2025 PISM Authors
#
# This file is part of PISM.
#
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

from netCDF4 import Dataset
from yac import *
import numpy as np
import json
from mpi4py import MPI
from enum import Enum

# Actions that the server can handle
# These have to match the actions which are defined on the client
class ServerActions(Enum):
    CREATE_FILE = 0
    SET_FILE_DIMENSION = 1
    SET_FILE_ATTRIBUTES = 2
    START_YAC_INITIALIZATION = 3
    FINISH_YAC_INITIALIZATION = 4
    DEFINE_NON_SPATIAL_VARIABLE = 5
    DEFINE_SPATIAL_VARIABLE = 6
    SEND_NON_SPATIAL_VARIABLE = 7
    SEND_SPATIAL_VARIABLE = 8
    UPDATE_TIME_LENGTH = 9
    FINISH = 10

# YAC general component, grid and configuration variables
source_comp_name = "pism"
source_grid_name = "pism_grid"
target_grid_name = "pism_grid_output"
target_comp_name = "pism_output_server"
time_unit        = TimeUnit.ISO_FORMAT
time_reduction   = Reduction.TIME_NONE
       
# Class to hold all YAC variables and definitions
class YacWrapper:
    # In the constructor the YAC initialization happens, followed
    # by the creation of the intercommunicator, which is necessary for 
    # direct communication with the client and receiving its information 
    def __init__(self):
        self.interpolation_stack = None
        self.grid = None
        self.global_vertex_indices = None
        self.x_size = 0
        self.y_size = 0
        self.fields = {} 

        self.yac = YAC()
        self.component = self.yac.def_comp(target_comp_name)
        
        # Get the local server communicator and group
        local_comm  = self.component.comp_comm
        local_group  = local_comm.Get_group()

        # Get the global communicator and group (contatining all processes from client and server)
        global_comm = self.yac.get_comps_comm(["pism", "pism_output_server"])
        global_group = global_comm.Get_group()
        
        # Get the rank of the local leader in the global communicator
        local_leader_global_rank = local_group.Translate_ranks([0], global_group)[0]
        
        # Create buffers for exchanging global ranks for local component leaders
        sendbuf = np.array([local_leader_global_rank], dtype='i')
        recvbuf = np.empty(global_comm.Get_size(), dtype='i')
        
        # Exchange global ranks of local leaders
        global_comm.Allgather([sendbuf, MPI.INT], [recvbuf, MPI.INT])
        self.remote_leader = recvbuf[0]
        self.remote_size = global_comm.Get_size() - local_comm.Get_size()
        
        # Create the intercommunicator with the received information
        self.intercomm = local_comm.Create_intercomm(
            0,
            global_comm,
            self.remote_leader,
            tag=0
        )
        
    # Starts the coupler initialization
    # For YAC specifically the main task happening in 
    # this subroutine is the grid definition
    def start_initialization(self):
        # Variables for the global x and y sizes
        x_size = np.empty(1, dtype='i')
        y_size = np.empty(1, dtype='i')
        client_local_domain_sizes = np.empty(self.remote_size, dtype='i')
        displacements = np.empty(self.remote_size, dtype='i')
        
        # Receive the global x and y sizes from the leader in the client
        self.intercomm.Recv([x_size, MPI.INT], source = self.remote_leader, tag = 0)
        self.intercomm.Recv([y_size, MPI.INT], source = self.remote_leader, tag = 0)
        self.x_size = x_size[0]
        self.y_size = y_size[0] 

        # Total amount of grid points
        grid_points = x_size[0] * y_size[0]
        
        # Gather the sizes of the local domains of all client processes
        self.intercomm.Gather(None, client_local_domain_sizes, root = MPI.ROOT)
        
        # Calculate the array displacements for receiving geometrical data from client processes
        displacements[0] = 0
        for i in range(1, self.remote_size):
            displacements[i] = displacements[i-1] + client_local_domain_sizes[i-1]
        
        # Create arrays for geometrical data and global data ordering
        longitudes = np.empty(grid_points, dtype='d')
        latitudes = np.empty(grid_points, dtype='d')
        self.global_vertex_indices = np.empty(grid_points, dtype='i')

        # Receive latitudes, longitudes and the global ordering of points
        self.intercomm.Gatherv(None, (self.global_vertex_indices, (client_local_domain_sizes, displacements)), root = MPI.ROOT)
        self.intercomm.Gatherv(None, (latitudes, (client_local_domain_sizes, displacements)), root = MPI.ROOT)
        self.intercomm.Gatherv(None, (longitudes, (client_local_domain_sizes, displacements)), root = MPI.ROOT)
        
        # Create grid and define point locations for corners (vertices)
        self.grid = CloudGrid(target_grid_name, longitudes, latitudes)
        self.grid.corner_points = self.grid.def_points(longitudes, latitudes)
        
        # Create the interpolation stack and add an nearest-neighbor inteporlation
        # For the current purposes of the output server we just want the data to be transferred
        # from the client to the server without any interpolation. Therefore, we use the NNN with a 
        # single neighbor, and since the client and server grids match, the data should arrive the same.
        self.interpolation_stack = InterpolationStack()
        self.interpolation_stack.add_nnn(NNNReductionType.AVG, 1, 1., 0.)

    # Finalize the YAC definitions phase
    def finish_initialization(self):
        self.yac.enddef()

    # Defines a YAC field for a spatial variable
    # In the future, when multiple files share gridded variables, either a
    # single field definition could be reused for multiple files or
    # individual dedicated fields could be created for each variable/file combination
    def define_field(self, variable_metadata):
        field_name = variable_metadata["variable_name"]
        timestep = variable_metadata["timestep"]
        collection_size = variable_metadata["collection_size"]

        self.fields[field_name] = Field.create(
            field_name, self.component, self.grid.corner_points, collection_size, timestep, TimeUnit.ISO_FORMAT
        )

        self.yac.def_couple(
        source_comp_name, source_grid_name, field_name,
        target_comp_name, target_grid_name, field_name,
        timestep, time_unit, time_reduction,
        self.interpolation_stack)

# Class to wrap the NetCDF functionalities for writing files
class OutputFile:
    # The constructor does empty initialization of most members and 
    # creates the dataset for the NetCDF file
    def __init__(self, file_name):
        self.name = file_name
        self.nc_dataset = Dataset(file_name[:-3] + "_server.nc", 'w')
        self.nc_variables = {}
        self.variables_metadata = {}
        self.dimensions = {} 
        self.comm_reqs = []
        self.text_req_indices = {}
        self.variables_data = {}
        self.time_index = 0 

    # Sets the attributes in the NetCDF dataset
    def set_attributes(self, file_attributes):
        for attr in file_attributes:
            setattr(self.nc_dataset, attr, file_attributes[attr])

    # Updates the time index used for writing to the datase
    def update_time_length(self, time_dimension_length):
        self.time_index = time_dimension_length
            
    # Creates a dimension in the NetCDF dataset
    def set_dimension(self, dimension):
        self.dimensions[dimension["dimension_name"]] = dimension["dimension_length"]
        self.nc_dataset.createDimension(dimension["dimension_name"],
                                        dimension["dimension_length"] if dimension["dimension_length"] > 0
                                        else None)
                                    
    # Defines a new variable in the NetCDF dataset 
    def define_variable(self, variable_metadata):
        field_name = variable_metadata["variable_name"]
        self.variables_metadata[field_name] = variable_metadata
    
        fill_value = variable_metadata["_FillValue"] if '_FillValue' in variable_metadata else None
        self.nc_variables[field_name] = self.nc_dataset.createVariable(field_name,
                                                                   variable_metadata["dtype"],
                                                                   variable_metadata["dimensions"],
                                                                   fill_value=fill_value)
    
        # There are some attributes which were either already set in the variable creation
        # or that are used only for the execution of the server. These all should be ignored
        # when setting NetCDF variable attributes.
        ignored_attributes = ['_FillValue', 'dimensions', 'tag', 'dtype', 
                              'collection_size', 'file_name', 'timestep', 'variable_name']
        for attr in variable_metadata:
            if attr not in ignored_attributes and variable_metadata[attr] != "":
                self.nc_variables[field_name].setncattr(attr, variable_metadata[attr])

    # Returns the metadata stored for that variable
    def get_variable_metadata(self, variable_name):
        return self.variables_metadata[variable_name]

    # Receives a spatial field through the coupler and writes its data to the corresponding NetCDF variable
    def receive_spatial_field(self, field_name, yac_wrapper):
        # Number of vertical levels for that field in YAC
        collection_size = yac_wrapper.fields[field_name].collection_size

        # Variable for holding the raw data provided by the YAC get, which is directly issued here
        data = yac_wrapper.fields[field_name].get()[0]

        # For each vertical level, reorders the received data according to the horizontal global vertex indices
        for level in range(collection_size):
            data[level] = data[level, np.argsort(yac_wrapper.global_vertex_indices)]

        # Creates a 3D array for reshapping the raw data received from YAC
        values = np.ndarray(shape=(collection_size, yac_wrapper.y_size, yac_wrapper.x_size),
                            buffer=data, dtype="f8")

        # Handles the data assignment to the NetCDF variable
        # If the variable is a 2D field NOT time-dependent
        if len(self.nc_variables[field_name].shape) == 2:
            self.nc_variables[field_name][:, :] = values[0]

        # If the variable is a time-dependent 2D field
        elif len(self.nc_variables[field_name].shape) == 3:
            self.nc_variables[field_name][self.time_index, :, :] = values[0]

        # If the variable is a time-dependent 3D field
        else:
            assert collection_size > 1
    
            # Temporary array for reordering the data for the NetCDF variable
            tmp = np.ndarray(shape=(yac_wrapper.y_size, yac_wrapper.x_size, collection_size),
                                dtype=values.dtype)
    
            for c in range(collection_size):
                tmp[:, :, c] = values[c]
    
            self.nc_variables[field_name][self.time_index, :, :, :] = tmp

    # Receives a non-spatial field through the coupler intercommunicator and
    # writes its data to the corresponding NetCDF variable
    def receive_non_spatial_field(self, variable_name, yac_wrapper):
        variable_dims = self.variables_metadata[variable_name]["dimensions"]
        variable_metadata = self.variables_metadata[variable_name]
        status = MPI.Status() 

        # If the variables has any dimensions defined
        if len(variable_dims) > 0:
            # Find the size for the receive operation
            size = self.dimensions[variable_dims[0]] if self.dimensions[variable_dims[0]] > 0 else 1
            tmp_receival = np.empty(size, dtype = "f8")
            tag = int(variable_metadata["tag"])
            # If the variable is not a string, receive it as double and cast it to the actual variable type
            if(variable_metadata["dtype"] != "S1"):
                yac_wrapper.intercomm.Recv([tmp_receival, MPI.DOUBLE], source = yac_wrapper.remote_leader, tag = tag)
                self.variables_data[variable_name] = tmp_receival.astype(variable_metadata["dtype"])
            # If the variable is a string then receive the char array and use the status var
            # The status is required to know how many chars were received for the NetCDF var assignment later
            else:
                self.variables_data[variable_name] = np.empty(size, dtype = variable_metadata["dtype"])
                yac_wrapper.intercomm.Recv([self.variables_data[variable_name], MPI.CHAR],
                                            source = yac_wrapper.remote_leader,
                                            tag = tag, status = status)

        # Handles the assignment of the data to the NetCDF variable
        # If it is NOT time-dependent
        if "time" not in variable_dims:
            # If it is a text variable copy only the amount of received
            # characters to avoid garbage at the end of the string
            if variable_metadata["dtype"] == "S1":
                char_count = status.Get_count(MPI.CHAR)
                self.nc_variables[variable_name][:char_count] = self.variables_data[variable_name][:char_count]

            # Otherwise simply assign all the received data
            else:
                self.nc_variables[variable_name][:] = self.variables_data[variable_name]

        # If the variable is time-dependent, assign it to the current time index of the file.
        # So far there does not seem to be time-dependent text variables
        else:
            self.nc_variables[variable_name][self.time_index] = self.variables_data[variable_name]

    # Closes the underlying NetCDF dataset
    def close(self):
        self.nc_dataset.close()

# Receives the json string for the metadata of an action
def receive_action_metadata(yac_wrapper):
    metadata_array_length = np.empty(1, dtype='i')

    # First receive the length of the string array and then the array itself
    yac_wrapper.intercomm.Recv([metadata_array_length, MPI.INT], source = yac_wrapper.remote_leader, tag = 0)
    metadata_array = np.empty(metadata_array_length[0], dtype='S1')
    yac_wrapper.intercomm.Recv([metadata_array, MPI.CHAR], source = yac_wrapper.remote_leader, tag = 0)

    # Decode the data and construct a dictionary from the json
    return json.loads(metadata_array.tobytes().decode("utf-8"))

np.set_printoptions(threshold=np.inf)
yac_wrapper = YacWrapper()
server_action = np.empty(1, dtype='i')
files = {}

# Poll loop for listening for action requests from the client
while True:
    # Wait for an action
    yac_wrapper.intercomm.Recv([server_action, MPI.INT], 
                                source = yac_wrapper.remote_leader, 
                                tag = 0)
    
    # Handle each action based on the action id
    match server_action[0]:
        case ServerActions.FINISH.value:
            break

        case ServerActions.CREATE_FILE.value:
            file_metadata = receive_action_metadata(yac_wrapper)
            files[file_metadata["file_name"]] = OutputFile(file_metadata["file_name"])

        case ServerActions.SET_FILE_ATTRIBUTES.value:
            file_attributes = receive_action_metadata(yac_wrapper)
            files[file_attributes["file_name"]].set_attributes(file_attributes["attributes"])

        case ServerActions.SET_FILE_DIMENSION.value:
            file_dimension = receive_action_metadata(yac_wrapper)
            files[file_dimension["file_name"]].set_dimension(file_dimension)

        case ServerActions.START_YAC_INITIALIZATION.value:
           yac_wrapper.start_initialization() 

        case ServerActions.FINISH_YAC_INITIALIZATION.value:
           yac_wrapper.finish_initialization() 

        case ServerActions.DEFINE_NON_SPATIAL_VARIABLE.value:
            variable_metadata = receive_action_metadata(yac_wrapper)
            files[variable_metadata["file_name"]].define_variable(variable_metadata)

        case ServerActions.DEFINE_SPATIAL_VARIABLE.value:
            variable_metadata = receive_action_metadata(yac_wrapper)
            files[variable_metadata["file_name"]].define_variable(variable_metadata)

            if variable_metadata["variable_name"] not in yac_wrapper.fields:
                yac_wrapper.define_field(variable_metadata)

        case ServerActions.SEND_SPATIAL_VARIABLE.value:
            variable_info = receive_action_metadata(yac_wrapper)
            files[variable_info["file_name"]].receive_spatial_field(variable_info["variable_name"], yac_wrapper)

        case ServerActions.SEND_NON_SPATIAL_VARIABLE.value:
            variable_info = receive_action_metadata(yac_wrapper)
            files[variable_info["file_name"]].receive_non_spatial_field(variable_info["variable_name"], yac_wrapper)

        case ServerActions.UPDATE_TIME_LENGTH.value:
            file_info = receive_action_metadata(yac_wrapper)
            files[file_info["file_name"]].update_time_length(file_info["time_dimension_length"]);

# Close all the files
for file in files.values():
    file.close()
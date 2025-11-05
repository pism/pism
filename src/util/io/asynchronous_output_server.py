from netCDF4 import Dataset
from yac import *
import numpy as np
import json
from mpi4py import MPI
from enum import Enum

# Actions that the server can handle
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
time_reduction   = 0
       
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
    
        ignored_attributes = ['_FillValue', 'dimensions', 'tag', 'dtype', 
                              'collection_size', 'file_name', 'timestep', 'variable_name']
        for attr in variable_metadata:
            if attr not in ignored_attributes and variable_metadata[attr] != "":
                self.nc_variables[field_name].setncattr(attr, variable_metadata[attr])

    def receive_spatial_field(self, field_name, yac_wrapper):
        var_dims = self.variables_metadata[field_name]["dimensions"]
        collection_size = yac_wrapper.fields[field_name].collection_size
    
        values = []
        data = yac_wrapper.fields[field_name].get()[0]

        for level in range(collection_size):
            data[level] = data[level, np.argsort(yac_wrapper.global_vertex_indices)]

        values = np.ndarray(shape=(collection_size, yac_wrapper.y_size, yac_wrapper.x_size),
                            buffer=data, dtype="f8")

        if len(self.nc_variables[field_name].shape) == 2:
            self.nc_variables[field_name][:, :] = values[0]
        elif len(self.nc_variables[field_name].shape) == 3:
            self.nc_variables[field_name][self.time_index, :, :] = values[0]
        else:
            assert collection_size > 1
    
            tmp = np.ndarray(shape=(yac_wrapper.y_size, yac_wrapper.x_size, collection_size),
                                dtype=values.dtype)
    
            for c in range(collection_size):
                tmp[:, :, c] = values[c]
    
            self.nc_variables[field_name][self.time_index, :, :, :] = tmp

    def get_variable_metadata(self, variable_name):
        return self.variables_metadata[variable_name]

    def receive_non_spatial_field(self, variable_name, yac_wrapper):
        var_dims = self.variables_metadata[variable_name]["dimensions"]
        var_metadata = self.variables_metadata[variable_name]
        status = MPI.Status() 
        if len(var_dims) > 0:
            size = self.dimensions[var_dims[0]]
            if len(var_dims) == 1 and var_dims[0] == "time" and size == 0:
                size = 1

        if var_dims[0] == "time" or len(var_dims) == 1:
            tmp_receival = np.empty(size, dtype = "f8")
            tag = int(var_metadata["tag"])
            if(var_metadata["dtype"] != "S1"):
                yac_wrapper.intercomm.Recv([tmp_receival, MPI.DOUBLE], source = yac_wrapper.remote_leader, tag = tag)
                self.variables_data[variable_name] = tmp_receival.astype(var_metadata["dtype"])
            else:
                self.variables_data[variable_name] = np.empty(size, dtype = var_metadata["dtype"])
                yac_wrapper.intercomm.Recv([self.variables_data[variable_name], MPI.CHAR], source = yac_wrapper.remote_leader, tag = tag, status = status)

        if "time" not in var_dims:
            if var_metadata["dtype"] == "S1":
                char_count = status.Get_count(MPI.CHAR)
                self.nc_variables[variable_name][:char_count] = self.variables_data[variable_name][:char_count]
            else:
                self.nc_variables[variable_name][:] = self.variables_data[variable_name]
        else:
            self.nc_variables[variable_name][self.time_index] = self.variables_data[variable_name]

    def close(self):
        self.nc_dataset.close()

def receive_action_metadata_string(yac_wrapper):
    metadata_array_length = np.empty(1, dtype='i')

    yac_wrapper.intercomm.Recv([metadata_array_length, MPI.INT], source = yac_wrapper.remote_leader, tag = 0)
    metadata_array = np.empty(metadata_array_length[0], dtype='S1')
    yac_wrapper.intercomm.Recv([metadata_array, MPI.CHAR], source = yac_wrapper.remote_leader, tag = 0)

    return metadata_array.tobytes().decode("utf-8")

np.set_printoptions(threshold=np.inf)
yac_wrapper = YacWrapper()
server_action = np.empty(1, dtype='i')
yac_data = None
files = {}

while True:
    yac_wrapper.intercomm.Recv([server_action, MPI.INT], 
                                source = yac_wrapper.remote_leader, 
                                tag = 0)
    
    match server_action[0]:
        case ServerActions.FINISH.value:
            break

        case ServerActions.CREATE_FILE.value:
            file_name = receive_action_metadata_string(yac_wrapper)
            files[file_name] = OutputFile(file_name)

        case ServerActions.SET_FILE_ATTRIBUTES.value:
            file_attributes = json.loads(receive_action_metadata_string(yac_wrapper))
            files[file_attributes["file_name"]].set_attributes(file_attributes["attributes"])

        case ServerActions.SET_FILE_DIMENSION.value:
            file_dimension = json.loads(receive_action_metadata_string(yac_wrapper))
            files[file_dimension["file_name"]].set_dimension(file_dimension)

        case ServerActions.START_YAC_INITIALIZATION.value:
           yac_wrapper.start_initialization() 

        case ServerActions.FINISH_YAC_INITIALIZATION.value:
           yac_wrapper.finish_initialization() 

        case ServerActions.DEFINE_NON_SPATIAL_VARIABLE.value:
            variable_metadata = json.loads(receive_action_metadata_string(yac_wrapper))
            files[variable_metadata["file_name"]].define_variable(variable_metadata)

        case ServerActions.DEFINE_SPATIAL_VARIABLE.value:
            variable_metadata = json.loads(receive_action_metadata_string(yac_wrapper))
            files[variable_metadata["file_name"]].define_variable(variable_metadata)

            if variable_metadata["variable_name"] not in yac_wrapper.fields:
                yac_wrapper.define_field(variable_metadata)

        case ServerActions.SEND_SPATIAL_VARIABLE.value:
            variable_info = json.loads(receive_action_metadata_string(yac_wrapper))
            files[variable_info["file_name"]].receive_spatial_field(variable_info["variable_name"], yac_wrapper)

        case ServerActions.SEND_NON_SPATIAL_VARIABLE.value:
            variable_info = json.loads(receive_action_metadata_string(yac_wrapper))
            files[variable_info["file_name"]].receive_non_spatial_field(variable_info["variable_name"], yac_wrapper)

        case ServerActions.UPDATE_TIME_LENGTH.value:
            file_info = json.loads(receive_action_metadata_string(yac_wrapper))
            files[file_info["file_name"]].update_time_length(file_info["time_dimension_length"]);

for file in files.values():
    file.close()
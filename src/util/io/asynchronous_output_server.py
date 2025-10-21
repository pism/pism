from netCDF4 import Dataset
from yac import *
import numpy as np
import json
from mpi4py import MPI
from enum import Enum

class ServerActions(Enum):
    FINISH = 0
    CREATE_FILE = 1
    SET_FILE_DIMENSION = 2
    INIT_YAC_GRID = 3
    SET_FILE_ATTRIBUTES = 4
    DEFINE_NON_SPATIAL_VARIABLE = 5
    SEND_NON_SPATIAL_VARIABLE = 6
    FINISH_YAC_INITIALIZATION = 7
    DEFINE_SPATIAL_VARIABLE = 8
    SEND_SPATIAL_VARIABLE = 9
    UPDATE_TIME_LENGTH = 10

start_datetime   = "1850-01-01T00:00:00"
end_datetime     = "1850-01-01T00:06:00"
src_comp_name    = "pism"
src_grid_name    = "pism_grid"
target_grid_name = "pism_grid_output"
target_comp_name = "pism_output_server"
time_unit        = TimeUnit.ISO_FORMAT
time_reduction   = 0
fields           = {}

class OutputFile:
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

    def set_attributes(self, file_attributes):
        for attr in file_attributes:
            setattr(self.nc_dataset, attr, file_attributes[attr])

    def update_time_length(self, time_dimension_length):
        self.time_index = time_dimension_length
            
    def set_dimension(self, dimension):
        self.dimensions[dimension["dimension_name"]] = dimension["dimension_length"]
        self.nc_dataset.createDimension(dimension["dimension_name"],
                                        dimension["dimension_length"] if dimension["dimension_length"] > 0
                                        else None)
                                    
    def define_variable(self, variable_metadata):
        attributes = variable_metadata 
        field_name = variable_metadata["variable_name"]
        self.variables_metadata[field_name] = variable_metadata
    
        fill_value = attributes["_FillValue"] if '_FillValue' in attributes else None
    
        self.nc_variables[field_name] = self.nc_dataset.createVariable(field_name,
                                                                   attributes["dtype"],
                                                                   attributes["dimensions"],
                                                                   fill_value=fill_value)
    
        special = ['_FillValue', 'dimensions', 'tag', 'dtype', 'collection_size', 'file_name', 'timestep', 'variable_name']
        for attr in attributes:
            if attr not in special and attributes[attr] != "":
                self.nc_variables[field_name].setncattr(attr, attributes[attr])

    def define_yac_field(self, variable_metadata, grid, interpolation_stack, x_size, y_size):
        field_name = variable_metadata["variable_name"]
        timestep = variable_metadata["timestep"] 
        collection_size = variable_metadata["collection_size"] 
    
        fields[field_name] = Field.create(
            field_name, component, grid.corner_points, collection_size, timestep, TimeUnit.ISO_FORMAT
        )
    
        yac.def_couple(
        src_comp_name, src_grid_name, field_name,
        target_comp_name, target_grid_name, field_name,
        timestep, time_unit, time_reduction,
        interpolation_stack)

    def receive_spatial_field(self, field_name, global_vertex_indices):
        var_dims = self.variables_metadata[field_name]["dimensions"]
        collection_size = fields[field_name].collection_size
    
        values = []
        data = fields[field_name].get()[0]

        for level in range(collection_size):
            data[level] = data[level, np.argsort(global_vertex_indices)]

        values = np.ndarray(shape=(collection_size, y_size, x_size),
                            buffer=data, dtype="f8")

        if len(self.nc_variables[field_name].shape) == 2:
            self.nc_variables[field_name][:, :] = values[0]
        elif len(self.nc_variables[field_name].shape) == 3:
            self.nc_variables[field_name][self.time_index, :, :] = values[0]
        else:
            assert collection_size > 1
    
            tmp = np.ndarray(shape=(y_size, x_size, collection_size),
                                dtype=values.dtype)
    
            for c in range(collection_size):
                tmp[:, :, c] = values[c]
    
            self.nc_variables[field_name][self.time_index, :, :, :] = tmp

    def get_variable_metadata(self, variable_name):
        return self.variables_metadata[variable_name]

    def receive_non_spatial_field(self, intercomm, variable_name): 
        var_dims = self.variables_metadata[variable_name]["dimensions"]
        var_metadata = self.variables_metadata[variable_name]
        status = MPI.Status() 
        if len(var_dims) > 0:
            size = self.dimensions[var_dims[0]]
            if len(var_dims) == 1 and var_dims[0] == "time" and size == 0:
                size = 1

        if var_dims[0] == "time" or len(var_dims) == 1:
            self.variables_data[variable_name] = np.empty(size, dtype = var_metadata["dtype"])
            tag = int(var_metadata["tag"])
            if(var_metadata["dtype"] == "f8"):
                intercomm.Recv([self.variables_data[variable_name], MPI.DOUBLE], source = remote_leader, tag = tag)
            elif(var_metadata["dtype"] == "i4"):
                intercomm.Recv([self.variables_data[variable_name], MPI.INT], source = remote_leader, tag = tag)
            else:
                intercomm.Recv([self.variables_data[variable_name], MPI.CHAR], source = remote_leader, tag = tag, status = status)

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

def receive_action_metadata_string():
    metadata_array_length = np.empty(1, dtype='i')

    intercomm.Recv([metadata_array_length, MPI.INT], source = remote_leader, tag = 0)
    metadata_array = np.empty(metadata_array_length[0], dtype='S1')
    intercomm.Recv([metadata_array, MPI.CHAR], source = remote_leader, tag = 0)

    return metadata_array.tobytes().decode("utf-8")

np.set_printoptions(threshold=np.inf)

yac = YAC()

component = yac.def_comp(target_comp_name)

global_comm = yac.get_comps_comm(["pism", "pism_output_server"])
local_comm  = component.comp_comm

global_group = global_comm.Get_group()
local_group  = local_comm.Get_group()

local_leader = local_group.Translate_ranks([0], global_group)[0]

sendbuf = np.array([local_leader], dtype='i')
recvbuf = np.empty(global_comm.Get_size(), dtype='i')

global_comm.Allgather([sendbuf, MPI.INT], [recvbuf, MPI.INT])
remote_leader = recvbuf[0]
remote_size = global_comm.Get_size() - local_comm.Get_size()

intercomm = local_comm.Create_intercomm(
    0,
    global_comm,
    remote_leader,
    tag=0
)

server_action = np.empty(1, dtype='i')
interpolation_stack = None
grid = None
global_vertex_indices = None
x_size = 0
y_size = 0

def initialize_yac_grid():
    x_size = np.empty(1, dtype='i')
    y_size = np.empty(1, dtype='i')
    
    intercomm.Recv([x_size, MPI.INT], source = remote_leader, tag = 0)
    intercomm.Recv([y_size, MPI.INT], source = remote_leader, tag = 0)
    grid_points = x_size[0] * y_size[0]
    
    longitudes = np.empty(grid_points, dtype='d')
    latitudes = np.empty(grid_points, dtype='d')
    global_vertex_indices = np.empty(grid_points, dtype='i')
    displacements = np.empty(remote_size, dtype='i')
    gather_buf = np.empty(remote_size, dtype='i')
    
    intercomm.Gather(None, gather_buf, root = MPI.ROOT)
    
    displacements[0] = 0
    for i in range(1, remote_size):
        displacements[i] = displacements[i-1] + gather_buf[i-1]
    
    intercomm.Gatherv(None, (global_vertex_indices, (gather_buf, displacements)), root = MPI.ROOT)
    intercomm.Gatherv(None, (latitudes, (gather_buf, displacements)), root = MPI.ROOT)
    intercomm.Gatherv(None, (longitudes, (gather_buf, displacements)), root = MPI.ROOT)
    
    grid = CloudGrid(target_grid_name, longitudes, latitudes)
    grid.set_global_index(global_vertex_indices, Location.CORNER)
    
    vertex_points = grid.def_points(longitudes, latitudes)
    grid.corner_points = vertex_points
    
    interpolation_stack = InterpolationStack()
    interpolation_stack.add_nnn(NNNReductionType.AVG, 1, 1., 0.)

    yac.sync_def()

    return grid, interpolation_stack, global_vertex_indices, x_size[0], y_size[0]
    
def finish_yac_initialization():
    yac.enddef()

time_independent_var_values = {}
snapshot_counter = 0
fields_nc_var = {}
file_type_split = True
no_files_initialized = True
received_non_grid_time_independent = False

files = {}

while True:
    intercomm.Recv([server_action, MPI.INT], source = remote_leader, tag = 0)
    
    match server_action[0]:
        case ServerActions.FINISH.value:
            break

        case ServerActions.CREATE_FILE.value:
            file_name = receive_action_metadata_string()
            files[file_name] = OutputFile(file_name)

        case ServerActions.SET_FILE_ATTRIBUTES.value:
            file_attributes = json.loads(receive_action_metadata_string())
            files[file_attributes["file_name"]].set_attributes(file_attributes["attributes"])

        case ServerActions.SET_FILE_DIMENSION.value:
            file_dimension = json.loads(receive_action_metadata_string())
            files[file_dimension["file_name"]].set_dimension(file_dimension)

        case ServerActions.INIT_YAC_GRID.value:
            grid, interpolation_stack, global_vertex_indices, x_size, y_size = initialize_yac_grid()

        case ServerActions.FINISH_YAC_INITIALIZATION.value:
            finish_yac_initialization()

        case ServerActions.DEFINE_NON_SPATIAL_VARIABLE.value:
            variable_metadata = json.loads(receive_action_metadata_string())
            files[variable_metadata["file_name"]].define_variable(variable_metadata)

        case ServerActions.DEFINE_SPATIAL_VARIABLE.value:
            variable_metadata = json.loads(receive_action_metadata_string())
            files[variable_metadata["file_name"]].define_variable(variable_metadata)

            if variable_metadata["variable_name"] not in fields :
                files[variable_metadata["file_name"]].define_yac_field(variable_metadata, grid, interpolation_stack, x_size, y_size)

        case ServerActions.SEND_SPATIAL_VARIABLE.value:
            variable_info = json.loads(receive_action_metadata_string())
            files[variable_info["file_name"]].receive_spatial_field(variable_info["variable_name"], global_vertex_indices)

        case ServerActions.SEND_NON_SPATIAL_VARIABLE.value:
            variable_info = json.loads(receive_action_metadata_string())
            files[variable_info["file_name"]].receive_non_spatial_field(intercomm, variable_info["variable_name"])

        case ServerActions.UPDATE_TIME_LENGTH.value:
            file_info = json.loads(receive_action_metadata_string())
            files[file_info["file_name"]].update_time_length(file_info["time_dimension_length"]);

intercomm.Barrier()

for file in files.values():
    file.close()
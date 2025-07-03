from netCDF4 import Dataset
from yac import *
import numpy as np
import json
from mpi4py import MPI

start_datetime   = "1850-01-01T00:00:00"
end_datetime     = "1850-01-01T00:06:00"
src_comp_name    = "pism"
src_grid_name    = "pism_grid"
target_grid_name = "pism_grid_output"
target_comp_name = "pism_output_server"
time_unit        = TimeUnit.ISO_FORMAT
time_reduction   = 0
fields           = {}

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

continue_recieving = np.empty(1, dtype='i')

x_size = np.empty(1, dtype='i')
y_size = np.empty(1, dtype='i')

intercomm.Bcast(x_size, root = remote_leader)
intercomm.Bcast(y_size, root = remote_leader)
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
# print("latitudes: ", latitudes)
# print("longitudes: ", longitudes)

# latitudes = latitudes.reshape(x_size[0], y_size[0])
# longitudes = longitudes.reshape(x_size[0], y_size[0])

grid = CloudGrid(target_grid_name, longitudes, latitudes)
grid.set_global_index(global_vertex_indices, Location.CORNER)

vertex_points = grid.def_points(longitudes, latitudes)
grid.corner_points = vertex_points

interpolation_stack = InterpolationStack()
interpolation_stack.add_nnn(NNNReductionType.AVG, 1, 1., 0.)

yac.sync_def()

pism_field_names = yac.get_field_names("pism", "pism_grid")

for field_name in pism_field_names:
    timestep                    = yac.get_field_timestep("pism", "pism_grid", field_name)
    collection_size             = yac.get_field_collection_size("pism", "pism_grid", field_name)

    fields[field_name] = Field.create(
        field_name, component, grid.corner_points, collection_size, timestep, TimeUnit.ISO_FORMAT
    )

    yac.def_couple(
      src_comp_name, src_grid_name, field_name,
      target_comp_name, target_grid_name, field_name,
      timestep, time_unit, time_reduction,
      interpolation_stack)

yac.enddef()

fields_metadata     = {}
non_field_variables = json.loads(yac.get_component_metadata(src_comp_name))
dimensions          = json.loads(yac.get_grid_metadata(src_grid_name))

for variable in non_field_variables["non_field_variables"]:
    fields_metadata[variable] = non_field_variables["non_field_variables"][variable]

for field_name in pism_field_names:
    json_string                 = yac.get_field_metadata(src_comp_name, src_grid_name, field_name)
    fields_metadata[field_name] = json.loads(json_string)

#
# ###### NC FILE CREATION LOOP ########
#
time_independent_var_values = {}
snapshot_counter = 0
while True:
    intercomm.Bcast(continue_recieving, root = remote_leader)
    if continue_recieving[0] == False:
        break

    fields_nc_var = {}
    output_dataset = Dataset("snapshot_" + str(snapshot_counter) + ".nc", 'w')

    for dimension in dimensions:
        output_dataset.createDimension(dimension,
                                       dimensions[dimension] if dimensions[dimension] > 0
                                       else None)

    for attr, val in non_field_variables.get("global", {}).items():
        setattr(output_dataset, attr, val)

    for field_name in fields_metadata:
        if('_FillValue' in fields_metadata[field_name]):
            fields_nc_var[field_name] = output_dataset.createVariable(field_name,
                                                                      fields_metadata[field_name]["dtype"],
                                                                      fields_metadata[field_name]["dimensions"],
                                                                      fill_value = fields_metadata[field_name]["_FillValue"])
        else:
            fields_nc_var[field_name] = output_dataset.createVariable(field_name,
                                                                      fields_metadata[field_name]["dtype"],
                                                                      fields_metadata[field_name]["dimensions"])
#
        for attribute_name in fields_metadata[field_name]:
            if (attribute_name != '_FillValue' and
                attribute_name != 'dimensions' and
                attribute_name != 'output_units' and
                attribute_name != 'tag' and
                attribute_name != 'dtype'):
                if fields_metadata[field_name][attribute_name] != "":
                    fields_nc_var[field_name].setncattr(attribute_name,
                                                        fields_metadata[field_name][attribute_name])

#
# ###### DATA RECEIVAL - LOOP ########
#
    comm_reqs = []
    values_vars = {}
    for variable in non_field_variables["non_field_variables"]:
        var_dims = fields_metadata[variable]["dimensions"]
        if len(var_dims) > 0:
            size = dimensions[var_dims[0]]
            if len(var_dims) == 1 and var_dims[0] == "time":
                size = size + 1
            values_vars[variable] = np.empty(size, dtype = fields_metadata[variable]["dtype"])
            tag = int(fields_metadata[variable]["tag"])
            if(fields_metadata[variable]["dtype"] == "f8"):
                comm_reqs.append(intercomm.Irecv([values_vars[variable], MPI.DOUBLE], source = remote_leader, tag = tag))
            else:
                comm_reqs.append(intercomm.Irecv([values_vars[variable], MPI.INT], source = remote_leader, tag = tag))
#
    for field_name in pism_field_names:
        var_dims = fields_metadata[field_name]["dimensions"]

        values = []
        if field_name not in time_independent_var_values:
            data = fields[field_name].get()[0]

            # if(field_name == "enthalpy"):
            #     print("enthalpy server: ", data[0])

            if (fields[field_name].collection_size > 1):
                values = np.ndarray(shape = (x_size[0], y_size[0], fields[field_name].collection_size),
                                    buffer = data,
                                    dtype = "f8")
            else:
                values = np.ndarray(shape = (x_size[0], y_size[0]),
                                    buffer = data[0, :],
                                    dtype = "f8")

            if fields_metadata[field_name]["dtype"] != "f8":
                values = values.astype(fields_metadata[field_name]["dtype"])

            if "time" not in var_dims:
                time_independent_var_values[field_name] = values
        else:
            values = time_independent_var_values[field_name]

        if len(fields_nc_var[field_name].shape) == 2:
            fields_nc_var[field_name][:, :] = values
        elif len(fields_nc_var[field_name].shape) == 3:
            fields_nc_var[field_name][0, :, :] = values
        else:
            fields_nc_var[field_name][0, :, :, :] = values

    MPI.Request.Waitall(comm_reqs)
    for variable in values_vars:
        fields_nc_var[variable][:] = values_vars[variable]

    output_dataset.close()
    snapshot_counter = snapshot_counter + 1

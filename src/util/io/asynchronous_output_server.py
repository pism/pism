# Copyright (C) 2025, 2026 PISM Authors
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

import logging

import netCDF4
import yac
import numpy as np
import json
from mpi4py import MPI
from enum import Enum

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("pism_output")


class ServerActions(Enum):
    """Actions that the server can handle.

    These have to match the actions which are defined on the client
    """

    OPEN_FILE = 1
    CLOSE_FILE = 2
    DEFINE_DIMENSION = 3
    SET_FILE_ATTRIBUTES = 4
    APPEND_HISTORY = 5
    DEFINE_YAC_GRID = 6
    DEFINE_YAC_FIELD = 7
    FINISH_YAC_INITIALIZATION = 8
    DEFINE_VARIABLE = 9
    SEND_VARIABLE = 10
    SEND_GRIDDED_VARIABLE = 11
    APPEND_TIME = 12
    FINISH = 13


# YAC general component, grid and configuration variables
source_comp_name = "pism"
target_comp_name = "pism_output"


class YacWrapper:
    """Class holding all YAC variables and definitions."""

    # In the constructor the YAC initialization happens, followed
    # by the creation of the intercommunicator, which is necessary for
    # direct communication with the client and receiving its information
    def __init__(self):
        self.interpolation_stack = None

        # YAC grids, references by name
        self.grids = {}

        # vertex indexes for each grid (referenced by the name of the grid)
        #
        # These indexes are used to re-order data received from PISM
        self.vertex_indexes = {}

        # number of grid points (an array [x_size, y_size]) for a given grid (refenced by name)
        self.grid_size = {}

        # YAC fields (by name)
        self.fields = {}

        self.yac = yac.YAC()

        self.component = self.yac.def_comp(target_comp_name)

        # Get the local server communicator and group
        local_comm = self.component.comp_comm
        local_group = local_comm.Get_group()

        # Get the global communicator and group (contatining all processes from client and server)
        global_comm = self.yac.get_comps_comm([source_comp_name, target_comp_name])
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

    def _output_grid_name(self, pism_grid_name):
        """Return the name of the grid corresponding to PISM's 'pism_grid_name'."""
        return pism_grid_name + "_output"

    def define_yac_grid(self, metadata):
        """Define the YAC grid for a given variable."""
        client_local_domain_sizes = np.empty(self.remote_size, dtype='i')
        displacements = np.empty(self.remote_size, dtype='i')

        # the grid name in pism_output has to be different
        grid_name = metadata['grid_name']

        # Receive the global x and y sizes from the leader in the client
        grid_size = np.empty(2, dtype='i')
        self.intercomm.Recv([grid_size, MPI.INT], source=self.remote_leader, tag=0)
        self.grid_size[grid_name] = grid_size

        # Total number of grid points
        grid_points = grid_size[0] * grid_size[1]

        # Gather the sizes of the local domains of all client processes
        self.intercomm.Gather(None, client_local_domain_sizes, root=MPI.ROOT)

        # Calculate the array displacements for receiving geometrical data from client processes
        displacements[0] = 0
        for i in range(1, self.remote_size):
            displacements[i] = displacements[i-1] + client_local_domain_sizes[i-1]

        # Create arrays for geometrical data and global data ordering
        longitudes = np.empty(grid_points, dtype='d')
        latitudes = np.empty(grid_points, dtype='d')
        self.vertex_indexes[grid_name] = np.empty(grid_points, dtype='i')

        # Receive latitudes, longitudes and the global ordering of points
        self.intercomm.Gatherv(None, (self.vertex_indexes[grid_name],
                                      (client_local_domain_sizes, displacements)),
                               root=MPI.ROOT)
        self.intercomm.Gatherv(None, (latitudes, (client_local_domain_sizes, displacements)),
                               root=MPI.ROOT)
        self.intercomm.Gatherv(None, (longitudes, (client_local_domain_sizes, displacements)),
                               root=MPI.ROOT)

        # Create grid and define point locations for corners (vertices)
        grid = yac.CloudGrid(self._output_grid_name(grid_name), longitudes, latitudes)
        grid.corner_points = grid.def_points(longitudes, latitudes)

        self.grids[grid_name] = grid

        # Create the interpolation stack and add an nearest-neighbor inteporlation
        # For the current purposes of the output server we just want the data to be transferred
        # from the client to the server without any interpolation. Therefore, we use the NNN with a
        # single neighbor, and since the client and server grids match, the data should arrive the same.
        self.interpolation_stack = yac.InterpolationStack()
        self.interpolation_stack.add_nnn(yac.NNNReductionType.AVG, 1, 1., 0.)

    def finish_initialization(self):
        """Finalize the YAC definitions phase."""
        self.yac.enddef()

    # In the future, when multiple files share gridded variables, either a
    # single field definition could be reused for multiple files or
    # individual dedicated fields could be created for each variable/file combination
    def define_field(self, metadata):
        """Define a YAC field for a spatial variable."""
        field_name = metadata["variable_name"]
        timestep = metadata["timestep"]
        collection_size = metadata["collection_size"]

        grid_name = metadata["grid_name"]

        time_unit = yac.TimeUnit.ISO_FORMAT
        time_reduction = yac.Reduction.TIME_NONE

        corner_points = self.grids[grid_name].corner_points
        self.fields[field_name] = yac.Field.create(field_name, self.component,
                                                   corner_points, collection_size,
                                                   timestep, yac.TimeUnit.ISO_FORMAT)

        self.yac.def_couple(source_comp_name, grid_name, field_name,
                            target_comp_name, self._output_grid_name(grid_name), field_name,
                            timestep, time_unit, time_reduction,
                            self.interpolation_stack)


class OutputFile:
    """Class wrapping NetCDF functionalities for writing files."""

    # The constructor does empty initialization of most members and
    # creates the dataset for the NetCDF file
    def __init__(self, file_name, mode="w"):
        self.nc_dataset = netCDF4.Dataset(file_name, mode)

        if "time" in self.nc_dataset.dimensions:
            time_length = len(self.nc_dataset.dimensions["time"])
        else:
            time_length = 0

        self.time_index = time_length - 1 if time_length > 0 else 0

    def set_attributes(self, file_attributes):
        """Set attributes in the NetCDF dataset."""
        for attr, value in file_attributes.items():
            self.nc_dataset.setncattr(attr, value)

    def append_time(self, time_seconds):
        """Create one more time record by adding a value to the time dimension."""
        self.time_index = len(self.nc_dataset.dimensions["time"])
        self.nc_dataset.variables["time"][self.time_index] = time_seconds

    def append_history(self, string):
        """Append a `string` to the global attribute 'history'."""
        try:
            old_history = self.nc_dataset["history"]
        except IndexError:
            old_history = ""

        self.nc_dataset.setncattr("history", old_history + string)

    def define_dimension(self, metadata):
        """Create a dimension in the NetCDF dataset."""
        name = metadata["name"]
        length = metadata["length"]
        self.nc_dataset.createDimension(name, length if length > 0 else None)

    def define_variable(self, metadata):
        """Define a new variable in the NetCDF dataset."""
        name = metadata["variable_name"]
        attributes = metadata["attributes"]

        fill_value = attributes["_FillValue"] if '_FillValue' in attributes else None
        nc_variable = self.nc_dataset.createVariable(name,
                                                     metadata["dtype"],
                                                     metadata["dimensions"],
                                                     fill_value=fill_value)

        ignored_attributes = ['_FillValue']
        for attr, value in attributes.items():
            if attr not in ignored_attributes and value != "":
                nc_variable.setncattr(attr, value)

    def receive_gridded_variable(self, metadata, yac_wrapper):
        """Receive a spatial field and write to the corresponding NetCDF variable."""

        name = metadata["variable_name"]
        ndims = metadata["ndims"]
        time_dependent = metadata["time_dependent"]

        nc_variable = self.nc_dataset.variables[name]

        # Number of vertical levels for that field in YAC
        collection_size = yac_wrapper.fields[name].collection_size

        # Grid size
        grid_name = metadata["grid_name"]
        x_size, y_size = yac_wrapper.grid_size[grid_name]

        # Variable for holding the raw data provided by the YAC get, which is directly issued here
        data = yac_wrapper.fields[name].get()[0]

        # For each vertical level, reorders the received data according to the horizontal
        # global vertex indices
        for level in range(collection_size):
            data[level] = data[level, np.argsort(yac_wrapper.vertex_indexes[grid_name])]

        # Creates a 3D array for reshaping the raw data received from YAC
        values = np.ndarray(shape=(collection_size, y_size, x_size),
                            buffer=data, dtype="f8")

        if ndims == 2:
            # 2D variable
            if time_dependent:
                nc_variable[self.time_index, :, :] = values[0]
            else:
                nc_variable[:, :] = values[0]
        else:
            # 3D variable
            assert collection_size > 1

            # Temporary array for reordering the data for the NetCDF variable
            tmp = np.ndarray(shape=(y_size, x_size, collection_size),
                             dtype=values.dtype)

            for c in range(collection_size):
                tmp[:, :, c] = values[c]

            if time_dependent:
                nc_variable[self.time_index, :, :, :] = tmp
            else:
                nc_variable[:, :, :] = tmp

    def receive_variable(self, metadata, yac_wrapper):
        """Receive and write non-gridded variable data."""
        name = metadata["variable_name"]
        start = metadata["start"]
        count = metadata["count"]
        tag = int(metadata["tag"])

        dtype = "f8"
        mpi_dtype = MPI.DOUBLE
        if "text" in metadata:
            dtype = "S1"
            mpi_dtype = MPI.CHAR

        data_size = 1
        for c in count:
            data_size *= c

        tmp = np.empty(data_size, dtype=dtype)
        yac_wrapper.intercomm.Recv([tmp, mpi_dtype],
                                   source=yac_wrapper.remote_leader, tag=tag)

        # Convert start and count arrays into a list of slices (one per dimension).
        ndim = len(count)
        hyperslab = [slice(start[i], start[i] + count[i]) for i in range(ndim)]

        # Unfortunately we have to use __setitem__ instead of [] indexing notation because
        # the "*hyperslab" syntax (Unpacking Argument Lists) does not work in square
        # brackets (invalid syntax).
        #
        # Ugh. Oh well.
        self.nc_dataset.variables[name].__setitem__(*hyperslab, tmp)

    def close(self):
        """Close the underlying NetCDF dataset."""
        self.nc_dataset.close()


def receive_action(yac_wrapper):
    """Receive the JSON string encoding an action and its metadata."""
    comm = yac_wrapper.intercomm

    status = MPI.Status()
    comm.Probe(source=yac_wrapper.remote_leader, tag=0, status=status)

    length = status.Get_count(MPI.CHAR)

    string = np.empty(length, dtype='S1')
    comm.Recv([string, MPI.CHAR], source=yac_wrapper.remote_leader, tag=0)

    # Decode the data and construct a dictionary from the json string
    string = string.tobytes().decode("utf-8")
    try:
        return json.loads(string)
    except BaseException:
        logger.critical(f"failed to parse action string '{string}'")
        raise


np.set_printoptions(threshold=np.inf)
yac_wrapper = YacWrapper()
files = {}


def get_file(name):
    """Return the file object, opening it if necessary."""
    if name not in files:
        files[name] = OutputFile(name, mode="w")
    return files[name]


# Poll loop for listening for action requests from the client
while True:
    # Wait for an action
    message = receive_action(yac_wrapper)

    action_id = message['action']
    metadata = message['info']

    try:
        file_name = metadata['file_name']
    except BaseException:
        file_name = None

    # Handle each action based on the action id
    match action_id:
        case ServerActions.FINISH.value:
            logger.debug("DONE")
            break

        case ServerActions.OPEN_FILE.value:
            # This can be the used as the stub of the "append" action.
            logger.debug(f"OPEN_FILE {file_name}")
            F = OutputFile(file_name, mode="rw")
            files[file_name] = F

            time_length = np.empty(1, dtype="i")
            if "time" in F.dimensions:
                time_length[0] = len(F.dimensions["time"])
            else:
                time_length[0] = 0

            last_time = np.empty(1, dtype="f8")
            last_time[0] = F.variables["time"][F.time_index]

            # Send time length and last time to PISM:
            yac_wrapper.intercomm.Send([time_length, 1, MPI.INT], dest=0, tag=0)
            yac_wrapper.intercomm.Send([last_time, 1, MPI.DOUBLE], dest=0, tag=0)

        case ServerActions.CLOSE_FILE.value:
            files[file_name].close()
            del files[file_name]

        case ServerActions.SET_FILE_ATTRIBUTES.value:
            logger.debug(f"SET_FILE_ATTRIBUTES {file_name}")
            get_file(file_name).set_attributes(metadata["attributes"])

        case ServerActions.APPEND_HISTORY.value:
            logger.debug(f"APPEND_HISTORY {file_name}")
            get_file(file_name).append_history(metadata["history"])

        case ServerActions.DEFINE_DIMENSION.value:
            logger.debug(f"DEFINE_DIMENSION {metadata['name']} in {file_name}")
            get_file(file_name).define_dimension(metadata)

        case ServerActions.DEFINE_YAC_GRID.value:
            logger.debug(f"DEFINE_YAC_GRID {metadata['grid_name']}")
            yac_wrapper.define_yac_grid(metadata)

        case ServerActions.DEFINE_YAC_FIELD.value:
            logger.debug(f"DEFINE_YAC_FIELD {metadata['variable_name']}")

            if metadata["variable_name"] not in yac_wrapper.fields:
                yac_wrapper.define_field(metadata)

        case ServerActions.FINISH_YAC_INITIALIZATION.value:
            logger.debug("FINISH_YAC_INITIALIZATION")
            yac_wrapper.finish_initialization()

        case ServerActions.DEFINE_VARIABLE.value:
            logger.debug(f"DEFINE_VARIABLE {metadata['variable_name']} in {file_name}")
            get_file(file_name).define_variable(metadata)

        case ServerActions.SEND_GRIDDED_VARIABLE.value:
            logger.debug(f"SEND_GRIDDED_VARIABLE {metadata['variable_name']} in {file_name}")
            get_file(file_name).receive_gridded_variable(metadata, yac_wrapper)

        case ServerActions.SEND_VARIABLE.value:
            logger.debug(f"SEND_VARIABLE {metadata['variable_name']} in {file_name}")
            get_file(file_name).receive_variable(metadata, yac_wrapper)

        case ServerActions.APPEND_TIME.value:
            logger.debug(f"APPEND_TIME in {file_name}")
            get_file(file_name).append_time(metadata["time"])

# Close all the files
for file in files.values():
    file.close()

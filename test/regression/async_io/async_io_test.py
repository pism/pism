#!/usr/bin/env python3

import PISM
import numpy as np

ctx = PISM.Context()
config = ctx.config

def create_grid():
    # these values are not important, but they should all be different
    config.set_number("grid.Mx", 200)
    config.set_number("grid.My", 100)
    config.set_number("grid.Mz", 41)
    config.set_number("grid.Lx", 2e3)
    config.set_number("grid.Ly", 1e3)

    return PISM.Grid.FromOptions(ctx.ctx)

def global_attributes():
    strings = PISM.StringMap()
    strings["Conventions"] = "CF-1.12"
    numbers = PISM.DoubleVectorMap()
    numbers["test"] = [1, 2]

    return strings, numbers

def create_array2d(grid):

    array2d = PISM.Scalar(grid, "thk")
    array2d.metadata().units("meter")

    with PISM.vec.Access(array2d):
        for i, j in grid.points():
            x = grid.x(i) * 1e-6
            y = grid.y(j) * 1e-6

            array2d[i, j] = x**2 + y**2

    return array2d

def create_array3d(grid):
    z = [0, 1, 2, 3, 4, 5]
    array3d = PISM.Array3D(grid, "enthalpy", PISM.WITHOUT_GHOSTS, z)
    array3d.metadata().units("J kg^-1")

    z = np.array(z)
    with PISM.vec.Access(array3d):
        for i, j in grid.points():
            x = grid.x(i) * 1e-6
            y = grid.y(j) * 1e-6

            array3d.set_column(i, j, x**2 + y**2 + z**2)

    return array3d


def test_writer(grid, writer, output_filename):
    """Test an output writer.

    Note: this can be called once per writer since the YAC-based writer can be initialized
    only once.
    """
    array2d = create_array2d(grid)
    array3d = create_array3d(grid)

    array2d_no_time = create_array2d(grid)
    array2d_no_time.metadata().set_name("topg")
    array2d_no_time.metadata().set_time_dependent(False)

    # ensures that we test type conversion in the output server
    array2d_int = create_array2d(grid)
    array2d_int.metadata().set_name("mask")
    array2d_int.metadata().units("")
    array2d_int.metadata().set_output_type(PISM.PISM_INT)

    array3d_no_time = create_array3d(grid)
    array3d_no_time.metadata().set_name("temp")
    array3d_no_time.metadata().units("kelvin")
    array3d_no_time.metadata().set_time_dependent(False)

    variables = [v.metadata() for v in [array2d, array2d_no_time, array2d_int,
                                        array3d, array3d_no_time]]
    variables.append(PISM.config_metadata(ctx.config))

    writer.initialize(PISM.VariableSet(variables))

    output = PISM.OutputFile(writer, output_filename)
    output.set_global_attributes(*global_attributes())
    output.append_history("some history")

    # the time variable has to be defined *first*
    first_time = 10.0
    output.define_variable(ctx.time.metadata())

    PISM.define_variables(output, variables, grid.get_mapping_info(), False)

    output.append_time(first_time)
    array2d.write(output)
    array2d_no_time.write(output)
    array2d_int.write(output)

    array3d.write(output)
    array3d_no_time.write(output)

    PISM.write_config(ctx.config, "pism_config", output)

    output.sync()

    output.close()

    return first_time


def test_appending(grid, writer, filename, old_time, old_time_length):
    """Test appending to an existing file."""
    output = PISM.OutputFile(writer, filename)
    output.append()

    assert output.time_dimension_length() == old_time_length
    assert output.last_time_value() == old_time

    time = old_time + 1.0
    output.append_time(time)

    assert output.time_dimension_length() == 2
    assert output.last_time_value() == time

    # test the ability to ingnore definition of an existing variable
    array2d = create_array2d(grid)
    output.define_variable(array2d.metadata())

    # test writing to a new time record while appending
    array2d.write(output)

    output.close()

    return time


if __name__ == "__main__":

    grid = create_grid()

    async_writer = PISM.YacOutputWriter(ctx.com, ctx.config)
    sync_writer = PISM.SynchronousOutputWriter(ctx.com, ctx.config)

    test_writer(grid, async_writer, "output_async.nc")
    test_writer(grid, sync_writer, "output_sync.nc")

    # test appending using the async writer
    #
    # Note: here we create a file using the sync writer so that the async writer has to
    # get time dimension info from the file (cannot re-use it).
    time1 = test_writer(grid, sync_writer, "file1.nc")
    test_appending(grid, async_writer, "file1.nc", time1, 1)

    # test appending using the sync writer
    time1 = test_writer(grid, sync_writer, "file2.nc")
    test_appending(grid, sync_writer, "file2.nc", time1, 1)

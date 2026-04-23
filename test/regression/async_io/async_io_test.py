#!/usr/bin/env python3

"""Test PISM's asynchronous output code.

This script creates two output files: one using old synchronous code and the other using
asynchronous code, then compares results.

"""

import PISM
import numpy as np
import netCDF4

ctx = PISM.Context()
config = ctx.config


def create_grid(Mx=200, My=100):
    """Create a PISM grid that will be used to allocate 2D and 3D arrays."""
    # these values are not important, but the grid we use should not be square
    config.set_number("grid.Mx", Mx)
    config.set_number("grid.My", My)
    config.set_number("grid.Mz", 41)
    config.set_number("grid.Lx", 2e3)
    config.set_number("grid.Ly", 1e3)

    return PISM.Grid.FromOptions(ctx.ctx)


def global_attributes():
    """Create global attributes."""
    strings = PISM.StringMap()
    strings["Conventions"] = "CF-1.12"
    numbers = PISM.DoubleVectorMap()
    numbers["test"] = [1, 2]

    return strings, numbers


def create_array2d(grid):
    """Create a 2D array for testing."""
    array2d = PISM.Scalar(grid, "thk")
    array2d.metadata().units("meter")

    with PISM.vec.Access(array2d):
        for i, j in grid.points():
            x = grid.x(i) * 1e-6
            y = grid.y(j) * 1e-6

            array2d[i, j] = x**2 + y**2

    return array2d


def create_array3d(grid):
    """Create a 3D array for testing."""
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


def create_arrays(grid):
    """Create 2D and 3D arrays used for testing."""
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

    return [array2d, array2d_int, array3d, array3d_no_time]


def test_writer(writer, arrays, mapping_info, output_filename):
    """Test an output writer.

    Note: this can be called once per writer since the YAC-based writer can be initialized
    only once.
    """
    variables = [v.metadata() for v in arrays]
    variables.append(PISM.config_metadata(ctx.config))

    writer.initialize(PISM.VariableSet(variables))

    output = PISM.OutputFile(writer, output_filename)
    output.set_global_attributes(*global_attributes())
    output.append_history("some history")

    # the time variable has to be defined *first*
    first_time = 10.0
    output.define_variable(ctx.time.metadata())

    PISM.define_variables(output, variables, mapping_info, False)

    output.append_time(first_time)
    for a in arrays:
        a.write(output)

    PISM.write_config(ctx.config, "pism_config", output)

    output.sync()

    output.close()

    return first_time


def test_appending(writer, arrays, filename, old_time, old_time_length):
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
    array = arrays[0]
    output.define_variable(array.metadata())

    # test writing to a new time record while appending
    array.write(output)

    output.close()

    return time


def compare(file_name_1, file_name_2):
    """Compare contents of two NetCDF files."""
    with netCDF4.Dataset(file_name_1) as f1:
        with netCDF4.Dataset(file_name_2) as f2:

            # compare global attributes:
            assert f1.ncattrs() == f2.ncattrs()
            for attr in f1.ncattrs():
                print(f"comparing global attribute {attr}")
                a1 = f1.getncattr(attr)
                a2 = f2.getncattr(attr)
                np.testing.assert_array_equal(a1, a2)

            # compare dimensions:
            assert f1.dimensions.keys() == f2.dimensions.keys()
            for dimension in f1.dimensions.keys():
                print(f"comparing dimension {dimension}")
                assert len(f1.dimensions[dimension]) == len(f2.dimensions[dimension])

            # compare variables
            assert f1.variables.keys() == f2.variables.keys()
            for variable in f1.variables.keys():
                print(f"comparing variable {variable}")
                v1 = f1.variables[variable]
                v2 = f1.variables[variable]

                # compare variable attributes:
                assert v1.ncattrs() == v2.ncattrs()
                for attr in v1.ncattrs():
                    print(f"comparing attribute {variable}:{attr}")
                    a1 = v1.getncattr(attr)
                    a2 = v2.getncattr(attr)
                    np.testing.assert_array_equal(a1, a2)

                # compare variable values
                print(f"comparing values stored in variable {variable}")
                np.testing.assert_array_equal(v1[:], v2[:])


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Test for the asynchronous output writer")
    parser.add_argument(
        "-k",
        "--keep_files",
        dest="keep",
        action="store_true",
        help="keep output files",
    )

    options = parser.parse_args()

    files = []

    def cleanup():
        """Remove files created by this test."""
        if ctx.rank != 0:
            return

        import os
        for f in files:
            print(f"Removing {f}")
            os.remove(f)

    if not options.keep:
        import atexit
        atexit.register(cleanup)

    grid1 = create_grid()
    grid2 = create_grid(Mx=400, My=200)

    array2d_grid2 = create_array2d(grid2)
    array2d_grid2.metadata().set_name("thk_large")
    array2d_grid2.metadata().dimension("x").set_name("x2")
    array2d_grid2.metadata().dimension("y").set_name("y2")

    array3d_grid2 = create_array3d(grid2)
    array3d_grid2.metadata().set_name("enthalpy_large")
    array3d_grid2.metadata().dimension("x").set_name("x2")
    array3d_grid2.metadata().dimension("y").set_name("y2")

    arrays = create_arrays(grid1)
    arrays.append(array2d_grid2)
    arrays.append(array3d_grid2)

    async_writer = PISM.YacOutputWriter(ctx.com, ctx.config)
    sync_writer = PISM.SynchronousOutputWriter(ctx.com, ctx.config)

    test_writer(async_writer, arrays, grid1.get_mapping_info(), "output_async.nc")
    test_writer(sync_writer, arrays, grid1.get_mapping_info(), "output_sync.nc")
    files += ["output_async.nc", "output_sync.nc"]
    # test appending using the async writer
    #
    # Note: here we create a file using the sync writer so that the async writer has to
    # get time dimension info from the file (cannot re-use it).
    time1 = test_writer(sync_writer, arrays, grid1.get_mapping_info(), "file1.nc")
    test_appending(async_writer, arrays, "file1.nc", time1, 1)
    files.append("file1.nc")

    # test appending using the sync writer
    time1 = test_writer(sync_writer, arrays, grid1.get_mapping_info(), "file2.nc")
    test_appending(sync_writer, arrays, "file2.nc", time1, 1)
    files.append("file2.nc")

    # compare output files created using asynchronous and "regular" output code:
    if ctx.rank == 0:
        compare("output_sync.nc", "output_async.nc")

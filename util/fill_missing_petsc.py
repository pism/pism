#!/usr/bin/env python

# @package fill_missing
# @brief This script solves the Laplace equation as a method of filling holes in map-plane data.
#
# Uses an approximation to Laplace's equation
# 	@f[ \nabla^2 u = 0 @f]
# to smoothly replace missing values in two-dimensional NetCDF
# variables with the average of the ``nearby'' non-missing values.
#
# Here is hypothetical example, filling the missing values in the variables
# `topg` and `usurf` in `data.nc` :
# \code
# fill_missing.py -v topg,usurf data.nc data_smoothed.nc
# \endcode
# Generally variables should be filled one at a time.
#
# Each of the requested variables must have missing values specified
# using the _FillValue attribute.

import petsc4py
import sys
petsc4py.init(sys.argv)

from petsc4py import PETSc
import numpy as np


def assemble_matrix(mask):
    """Assemble the matrix corresponding to the standard 5-point stencil
    approximation of the Laplace operator on the domain defined by
    mask == True, where mask is a 2D NumPy array. The stencil wraps
    around the grid, i.e. this is an approximation of the Laplacian
    on a torus.

    The grid spacing is ignored, which is equivalent to assuming equal
    spacing in x and y directions.
    """
    PETSc.Sys.Print("Assembling the matrix...")
    # grid size
    nrow, ncol = mask.shape

    # create sparse matrix
    A = PETSc.Mat()
    A.create(PETSc.COMM_WORLD)
    A.setSizes([nrow * ncol, nrow * ncol])
    A.setType('aij')  # sparse
    A.setPreallocationNNZ(5)

    # precompute values for setting
    # diagonal and non-diagonal entries
    diagonal = 4.0
    offdx = - 1.0
    offdy = - 1.0

    def R(i, j):
        "Map from the (row,column) pair to the linear row number."
        return i * ncol + j

    # loop over owned block of rows on this
    # processor and insert entry values
    row_start, row_end = A.getOwnershipRange()
    for row in range(row_start, row_end):
        A[row, row] = diagonal

        i = row // ncol    # map row number to
        j = row - i * ncol  # grid coordinates

        if mask[i, j] == False:
            continue

        # i
        if i == 0:              # top row
            col = R(nrow - 1, j)
            A[row, col] = offdx

        if i > 0:               # interior
            col = R(i - 1, j)
            A[row, col] = offdx

        if i < nrow - 1:        # interior
            col = R(i + 1, j)
            A[row, col] = offdx

        if i == nrow - 1:       # bottom row
            col = R(0, j)
            A[row, col] = offdx

        # j
        if j == 0:              # left-most column
            col = R(i, ncol - 1)
            A[row, col] = offdy

        if j > 0:               # interior
            col = R(i, j - 1)
            A[row, col] = offdy

        if j < ncol - 1:        # interior
            col = R(i, j + 1)
            A[row, col] = offdy

        if j == ncol - 1:       # right-most column
            col = R(i, 0)
            A[row, col] = offdy

    # communicate off-processor values
    # and setup internal data structures
    # for performing parallel operations
    A.assemblyBegin()
    A.assemblyEnd()

    PETSc.Sys.Print("done.")
    return A


def assemble_rhs(rhs, X):
    """Assemble the right-hand side of the system approximating the
    Laplace equation.

    Modifies rhs in place; sets Dirichlet BC using X where X.mask ==
    False.
    """
    # PETSc.Sys.Print("Setting Dirichlet BC...")
    nrow, ncol = X.shape
    row_start, row_end = rhs.getOwnershipRange()

    # The right-hand side is zero everywhere except for Dirichlet
    # nodes.
    rhs.set(0.0)

    for row in range(row_start, row_end):
        i = row // ncol    # map row number to
        j = row - i * ncol  # grid coordinates

        if X.mask[i, j] == False:
            rhs[row] = 4.0 * X[i, j]

    rhs.assemble()
    # PETSc.Sys.Print("done.")


def create_solver():
    "Create the KSP solver"
    # create linear solver
    ksp = PETSc.KSP()
    ksp.create(PETSc.COMM_WORLD)

    # Use algebraic multigrid:
    pc = ksp.getPC()
    pc.setType(PETSc.PC.Type.GAMG)
    ksp.setFromOptions()

    ksp.setInitialGuessNonzero(True)

    return ksp


def fill_missing(field, matrix=None):
    """Fill missing values in a NumPy array 'field' using the matrix
    'matrix' approximating the Laplace operator."""

    ksp = create_solver()

    if matrix is None:
        A = assemble_matrix(field.mask)
    else:
        # PETSc.Sys.Print("Reusing the matrix...")
        A = matrix

    # obtain solution & RHS vectors
    x, b = A.getVecs()

    assemble_rhs(b, field)

    initial_guess = np.mean(field)

    # set the initial guess
    x.set(initial_guess)

    ksp.setOperators(A)

    # Solve Ax = b
    # PETSc.Sys.Print("Solving...")
    ksp.solve(b, x)
    # PETSc.Sys.Print("done.")

    # transfer solution to processor 0
    vec0, scatter = create_scatter(x)
    scatter_to_0(x, vec0, scatter)

    return vec0, A


def create_scatter(vector):
    "Create the scatter to processor 0."
    comm = vector.getComm()
    rank = comm.getRank()
    scatter, V0 = PETSc.Scatter.toZero(vector)
    scatter.scatter(vector, V0, False, PETSc.Scatter.Mode.FORWARD)
    comm.barrier()

    return V0, scatter


def scatter_to_0(vector, vector_0, scatter):
    "Scatter a distributed 'vector' to 'vector_0' on processor 0 using 'scatter'."
    comm = vector.getComm()
    scatter.scatter(vector, vector_0, False, PETSc.Scatter.Mode.FORWARD)
    comm.barrier()


def scatter_from_0(vector_0, vector, scatter):
    "Scatter 'vector_0' on processor 0 to a distributed 'vector' using 'scatter'."
    comm = vector.getComm()
    scatter.scatter(vector, vector_0, False, PETSc.Scatter.Mode.REVERSE)
    comm.barrier()


def fill_2d_record(data, matrix=None):
    "Fill missing values in a 2D record."

    if getattr(data, "mask", None) is None:
        return data, None
    filled_data, A = fill_missing(data, matrix)
    if PETSc.COMM_WORLD.getRank() == 0:
        filled_data = filled_data[:].reshape(data.shape)

    return filled_data, A


def test():
    "Test fill_missing() using synthetic data."
    N = 201
    M = N * 1.5
    x = np.linspace(-1, 1, N)
    y = np.linspace(-1, 1, M)
    xx, yy = np.meshgrid(x, y)
    zz = np.sin(2.5 * np.pi * xx) * np.cos(2.0 * np.pi * yy)

    K = 10
    mask = np.random.randint(0, K, zz.size).reshape(zz.shape) / float(K)

    mask[(xx - 1.0) ** 2 + (yy - 1.0) ** 2 < 1.0] = 1
    mask[(xx + 1.0) ** 2 + (yy + 1.0) ** 2 < 1.0] = 1

    field = np.ma.array(zz, mask=mask)

    zzz, _ = fill_missing(field)

    rank = PETSc.COMM_WORLD.getRank()
    if rank == 0:
        zzz0_np = zzz[:].reshape(field.shape)

        import pylab as plt

        plt.figure(1)
        plt.imshow(zz, interpolation='nearest')

        plt.figure(2)
        plt.imshow(field, interpolation='nearest')

        plt.figure(3)
        plt.imshow(zzz0_np, interpolation='nearest')

        plt.show()


def fill_variable(nc, name):
    "Fill missing values in one variable."
    PETSc.Sys.Print("Processing %s..." % name)
    t0 = time()

    var = nc.variables[name]

    comm = PETSc.COMM_WORLD
    rank = comm.getRank()

    if var.ndim == 3:
        A = None
        n_records = var.shape[0]
        for t in range(n_records):
            PETSc.Sys.Print("Processing record %d/%d..." % (t + 1, n_records))
            data = var[t, :, :]

            filled_data, A = fill_2d_record(data, A)
            if rank == 0:
                var[t, :, :] = filled_data

            comm.barrier()
        PETSc.Sys.Print("Time elapsed: %5f seconds." % (time() - t0))
    elif var.ndim == 2:
        data = var[:, :]

        filled_data, _ = fill_2d_record(data)
        if rank == 0:
            var[:, :] = filled_data

        comm.barrier()
        PETSc.Sys.Print("Time elapsed: %5f seconds." % (time() - t0))
    else:
        PETSc.Sys.Print("Skipping the %dD variable %s." % (var.ndim, name))
        return

    # Remove the _FillValue attribute:
    try:
        delattr(var, '_FillValue')
    except:
        pass

    # Remove the missing_value attribute:
    try:
        delattr(var, 'missing_value')
    except:
        pass


def add_history(nc):
    "Update the history attribute in a NetCDF file nc."
    comm = PETSc.COMM_WORLD
    rank = comm.getRank()

    if rank != 0:
        return

    # add history global attribute (after checking if present)
    historysep = ' '
    historystr = asctime() + ': ' + historysep.join(sys.argv) + '\n'
    if 'history' in nc.ncattrs():
        nc.history = historystr + nc.history  # prepend to history string
    else:
        nc.history = historystr


if __name__ == "__main__":
    from argparse import ArgumentParser
    import os
    import os.path
    import tempfile
    import shutil
    from time import time, asctime

    try:
        from netCDF4 import Dataset as NC
    except:
        PETSc.Sys.Print("netCDF4 is not installed!")
        sys.exit(1)

    parser = ArgumentParser()
    parser.description = "Fill missing values by solving the Laplace equation in on the missing values and using present values as Dirichlet B.C."

    parser.add_argument("INPUT", nargs=1, help="Input file name.")
    parser.add_argument("OUTPUT", nargs=1, help="Output file name.")
    parser.add_argument("-a", "--all", dest="all", action="store_true",
                        help="Process all variables.")
    parser.add_argument("-v", "--vars", dest="variables",
                        help="comma-separated list of variables to process")

    options, _ = parser.parse_known_args()

    input_filename = options.INPUT[0]
    output_filename = options.OUTPUT[0]

    if options.all:
        nc = NC(input_filename)
        variables = list(nc.variables.keys())
        nc.close()
    else:
        try:
            variables = (options.variables).split(',')
        except:
            PETSc.Sys.Print("Please specify variables using the -v option.")
            sys.exit(-1)

    # Done processing command-line options.

    comm = PETSc.COMM_WORLD
    rank = comm.getRank()

    t0 = time()

    PETSc.Sys.Print("Filling missing values in %s and saving results to %s..." % (input_filename,
                                                                                  output_filename))
    if rank == 0:
        try:
            PETSc.Sys.Print("Creating a temporary file...")
            # find the name of the directory with the output file:
            dirname = os.path.dirname(os.path.abspath(output_filename))
            (handle, tmp_filename) = tempfile.mkstemp(prefix="fill_missing_",
                                                      suffix=".nc",
                                                      dir=dirname)

            os.close(handle)  # mkstemp returns a file handle (which we don't need)
        except IOError:
            PETSc.Sys.Print("ERROR: Can't create %s, Exiting..." % tmp_filename)

        try:
            PETSc.Sys.Print("Copying input file %s to %s..." % (input_filename,
                                                                tmp_filename))
            shutil.copy(input_filename, tmp_filename)
        except IOError:
            PETSc.Sys.Print("ERROR: Can't copy %s, Exiting..." % input_filename)

    try:
        if rank == 0:
            nc = NC(tmp_filename, 'a')
        else:
            nc = NC(input_filename, 'r')
    except Exception as message:
        PETSc.Sys.Print(message)
        PETSc.Sys.Print("Note: %s was not modified." % output_filename)
        sys.exit(-1)

    add_history(nc)

    for name in variables:
        try:
            fill_variable(nc, name)
        except Exception as message:
            PETSc.Sys.Print("ERROR:", message)
            PETSc.Sys.Print("Note: %s was not modified." % output_filename)
            sys.exit(-1)
    nc.close()

    try:
        if rank == 0:
            shutil.move(tmp_filename, output_filename)
    except:
        PETSc.Sys.Print("Error moving %s to %s. Exiting..." % (tmp_filename,
                                                               output_filename))
        sys.exit(-1)

    PETSc.Sys.Print("Total time elapsed: %5f seconds." % (time() - t0))

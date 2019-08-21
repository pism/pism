#!/usr/bin/env python

# @package fill_missing
# \brief This script solves the Laplace equation as a method of filling holes in map-plane data.
#
# \details The script is an implementation of the SOR method with Chebyshev
# acceleration for the Laplace equation, as described in 'Numerical Recipes in
# Fortran: the art of scientific computing' by William H. Press et al -- 2nd
# edition.
#
# Note also that this script can be used both from the command line and as a
# Python module -- by adding 'from fill_missing import laplace' to your
# program.
# Uses an approximation to Laplace's equation
# 	\f[ \nabla^2 u = 0 \f]
# to smoothly replace missing values in two-dimensional NetCDF variables with the average of the ``nearby'' non-missing values.
# Here is hypothetical example, filling the missing values in the variables \c topg and \c usurf, using a convergence tolerance of \f$10^{-4}\f$ and the initial guess of \f$100\f$, on data in the NetCDF file \c data.nc :
# \code
# fill_missing.py -f data.nc -v topg,usurf --eps=1.0e-4 \
#                 -i 100.0 -o data_smoothed.nc
# \endcode
# Options \c -i and \c -e specify the initial guess and the convergence tolerance for \e all the specified variables, so using these options only makes sense if all the variables have the same units.  Moreover, making a good initial guess can noticeably reduce the time needed to fill in the holes.  Generally variables should be filled one at a time.
#
# Each of the requested variables must have missing values specified
# according to CF Metadata conventions, namely one of the following:
# \c valid_range or both of \c valid_min and
# \c valid_max (if the values are in a specific range); one of
# \c valid_min (\c valid_max) if values are greater (less)
# than some value, or \c _FillValue.  Also \c _FillValue is
# interpreted as \c valid_max if it is positive, and as
# \c valid_min otherwise, and the \c missing_value  attribute is deprecated
# by the NetCDF User's Guide, but is supported for backward compatibility.  For more information see
# <a href="http://www.unidata.ucar.edu/software/netcdf/guide_10.html#SEC76">NetCDF User's Guide: Attributes</a>.
# Run \verbatim fill_missing.py --help \endverbatim for the list of available
# command-line options.


# CK, 08/12/2008

from numpy import *

# Computes \f$\rho_{Jacobi}\f$, see formula (19.5.24), page 858.


def rho_jacobi(dimensions):
    (J, L) = dimensions
    return (cos(pi / J) + cos(pi / L)) / 2

# This makes the stencil wrap around the grid. It is unclear if this should be
#  done, but it allows using a 4-point stencil for all points, even if they
#  are on the edge of the grid (otherwise we need to use three points on the
#  sides and two in the corners).
#
#  Is and Js are arrays with row- and column-indices, M and N are the grid
#  dimensions.


def fix_indices(Is, Js, dimensions):

    (M, N) = dimensions
    Is[Is == M] = 0
    Is[Is == -1] = M - 1
    Js[Js == N] = 0
    Js[Js == -1] = N - 1
    return (Is, Js)

# \brief laplace solves the Laplace equation
# \details laplace solves the Laplace equation using the SOR method with Chebyshev
#    acceleration as described in 'Numerical Recipes in Fortran: the art of
#    scientific computing' by William H. Press et al -- 2nd edition, section
#    19.5.
#
#    data is a 2-d array (computation grid)
#
#    mask is a boolean array; setting mask to 'data == 0', for example, results
#         in only modifying points where 'data' is zero, all the other points
#         are left as is. Intended use: if in an array the value of -9999.0
#         signifies a missing value, then setting mask to 'data == -9999.0'
#         fills in all the missing values.
#
#    eps1 is the first stopping criterion: the iterations stop if the norm of
#         residual becomes less than eps1*initial_norm, where 'initial_norm' is
#         the initial norm of residual. Setting eps1 to zero or a negative
#         number disables this stopping criterion.
#
#    eps2 is the second stopping criterion: the iterations stop if the absolute
#         value of the maximal change in value between successive iterations is
#         less than eps2. Setting eps2 to zero or a negative number disables
#         this stopping criterion.
#
#    initial_guess is the initial guess used for all the values in the domain;
#         the default is 'mean', i.e. use the mean of all the present values as
#         the initial guess for missing values. initial_guess has to be 'mean'
#         or a number.
#
#    max_iter is the maximum number of iterations allowed. The default is 10000.


def laplace(data, mask, eps1, eps2, initial_guess='mean', max_iter=10000):

    dimensions = data.shape
    rjac = rho_jacobi(dimensions)
    i, j = indices(dimensions)
    # This splits the grid into 'odd' and 'even' parts, according to the
    # checkerboard pattern:
    odd = (i % 2 == 1) ^ (j % 2 == 0)
    even = (i % 2 == 0) ^ (j % 2 == 0)
    # odd and even parts _in_ the domain:
    odd_part = list(zip(i[mask & odd], j[mask & odd]))
    even_part = list(zip(i[mask & even], j[mask & even]))
    # relative indices of the stencil points:
    k = array([0, 1, 0, -1])
    l = array([-1, 0, 1, 0])
    parts = [odd_part, even_part]

    try:
        initial_guess = float(initial_guess)
    except:
        if initial_guess == 'mean':
            present = array(ones_like(mask) - mask, dtype=bool)
            initial_guess = mean(data[present])
        else:
            print("""ERROR: initial_guess of '%s' is not supported (it should be a number or 'mean').
    Note: your data was not modified.""" % initial_guess)
            return

    data[mask] = initial_guess
    print("Using the initial guess of %10f." % initial_guess)

    # compute the initial norm of residual
    initial_norm = 0.0
    for m in [0, 1]:
        for i, j in parts[m]:
            Is, Js = fix_indices(i + k, j + l, dimensions)
            xi = sum(data[Is, Js]) - 4 * data[i, j]
            initial_norm += abs(xi)
    print("Initial norm of residual =", initial_norm)
    print("Criterion is (change < %f) OR (res norm < %f (initial norm))." % (eps2, eps1))

    omega = 1.0
    # The main loop:
    for n in arange(max_iter):
        anorm = 0.0
        change = 0.0
        for m in [0, 1]:
            for i, j in parts[m]:
                # stencil points:
                Is, Js = fix_indices(i + k, j + l, dimensions)
                residual = sum(data[Is, Js]) - 4 * data[i, j]
                delta = omega * 0.25 * residual
                data[i, j] += delta

                # record the maximal change and the residual norm:
                anorm += abs(residual)
                if abs(delta) > change:
                    change = abs(delta)
                # Chebyshev acceleration (see formula 19.5.30):
                if n == 1 and m == 1:
                    omega = 1.0 / (1.0 - 0.5 * rjac ** 2)
                else:
                    omega = 1.0 / (1.0 - 0.25 * rjac ** 2 * omega)
        print("max change = %10f, residual norm = %10f" % (change, anorm))
        if (anorm < eps1 * initial_norm) or (change < eps2):
            print("Exiting with change=%f, anorm=%f after %d iteration(s)." % (change,
                                                                               anorm, n + 1))
            return
    print("Exceeded the maximum number of iterations.")
    return


if __name__ == "__main__":
    from optparse import OptionParser
    from sys import argv, exit
    from shutil import copy, move
    from tempfile import mkstemp
    from os import close
    from time import time, asctime
    try:
        from netCDF4 import Dataset as NC
    except:
        print("netCDF4 is not installed!")
        sys.exit(1)

    parser = OptionParser()

    parser.usage = "%prog [options]"
    parser.description = "Fills missing values in variables selected using -v in the file given by -f."
    parser.add_option("-f", "--file", dest="input_filename",
                      help="input file")
    parser.add_option("-v", "--vars", dest="variables",
                      help="comma-separated list of variables to process")
    parser.add_option("-o", "--out_file", dest="output_filename",
                      help="output file")
    parser.add_option("-e", "--eps", dest="eps",
                      help="convergence tolerance",
                      default="1.0")
    parser.add_option("-i", "--initial_guess", dest="initial_guess",
                      help="initial guess to use; applies to all selected variables",
                      default="mean")

    (options, args) = parser.parse_args()

    if options.input_filename == "":
        print("""Please specify the input file name
(using the -f or --file command line option).""")
        exit(-1)
    input_filename = options.input_filename

    if options.variables == "":
        print("""Please specify the list of variables to process
(using the -v or --variables command line option).""")
        exit(-1)
    variables = (options.variables).split(',')

    if options.output_filename == "":
        print("""Please specify the output file name
(using the -o or --out_file command line option).""")
        exit(-1)
    output_filename = options.output_filename

    eps = float(options.eps)

    # Done processing command-line options.

    print("Creating the temporary file...")
    try:
        (handle, tmp_filename) = mkstemp()
        close(handle)  # mkstemp returns a file handle (which we don't need)
        copy(input_filename, tmp_filename)
    except IOError:
        print("ERROR: Can't create %s, Exiting..." % tmp_filename)

    try:
        nc = NC(tmp_filename, 'a')
    except Exception as message:
        print(message)
        print("Note: %s was not modified." % output_filename)
        exit(-1)

    # add history global attribute (after checking if present)
    historysep = ' '
    historystr = asctime() + ': ' + historysep.join(argv) + '\n'
    if 'history' in nc.ncattrs():
        nc.history = historystr + nc.history  # prepend to history string
    else:
        nc.history = historystr

    t_zero = time()
    for name in variables:
        print("Processing %s..." % name)
        try:
            var = nc.variables[name]

            attributes = ["valid_range", "valid_min", "valid_max",
                          "_FillValue", "missing_value"]
            adict = {}
            print("Reading attributes...")
            for attribute in attributes:
                print("* %15s -- " % attribute, end=' ')
                if attribute in var.ncattrs():
                    adict[attribute] = getattr(var, attribute)
                    print("found")
                else:
                    print("not found")

            if (var.ndim == 3):
                nt = var.shape[0]
                for t in range(0, nt):
                    print("\nInterpolating time step %i of %i\n" % (t, nt))

                    data = asarray(squeeze(var[t, :, :].data))

                    if "valid_range" in adict:
                        range = adict["valid_range"]
                        mask = ((data >= range[0]) & (data <= range[1]))
                        print("Using the valid_range attribute; range = ", range)

                    elif "valid_min" in adict and "valid_max" in adict:
                        valid_min = adict["valid_min"]
                        valid_max = adict["valid_max"]
                        mask = ((data < valid_min) | (data > valid_max))
                        print("""Using valid_min and valid_max attributes.
        valid_min = %10f, valid_max = %10f.""" % (valid_min, valid_max))

                    elif "valid_min" in adict:
                        valid_min = adict["valid_min"]
                        mask = data < valid_min
                        print("Using the valid_min attribute; valid_min = %10f" % valid_min)

                    elif "valid_max" in adict:
                        valid_max = adict["valid_max"]
                        mask = data > valid_max
                        print("Using the valid_max attribute; valid_max = %10f" % valid_max)

                    elif "_FillValue" in adict:
                        fill_value = adict["_FillValue"]
                        if fill_value <= 0:
                            mask = data <= fill_value + 2 * finfo(float).eps
                        else:
                            mask = data >= fill_value - 2 * finfo(float).eps
                        print("Using the _FillValue attribute; _FillValue = %10f" % fill_value)

                    elif "missing_value" in adict:
                        missing = adict["missing_value"]
                        mask = abs(data - missing) < 2 * finfo(float).eps
                        print("""Using the missing_value attribute; missing_value = %10f
        Warning: this attribute is deprecated by the NUG.""" % missing)

                    else:
                        print("No missing values found. Skipping this variable...")
                        continue

                    count = int(sum(mask))
                    if count == 0:
                        print("No missing values found. Skipping this variable...")
                        continue
                    print("Filling in %5d missing values..." % count)
                    t0 = time()
                    laplace(data, mask, -1, eps, initial_guess=options.initial_guess)
                    var[t, :, :] = data

                    # now REMOVE missing_value and _FillValue attributes
                    try:
                        delattr(var, '_FillValue')
                    except:
                        pass
                    try:
                        delattr(var, 'missing_value')
                    except:
                        pass
                    print("This took %5f seconds." % (time() - t0))

            elif (var.ndim == 2):

                data = asarray(squeeze(var[:]))

                if "valid_range" in adict:
                    range = adict["valid_range"]
                    mask = ((data >= range[0]) & (data <= range[1]))
                    print("Using the valid_range attribute; range = ", range)

                elif "valid_min" in adict and "valid_max" in adict:
                    valid_min = adict["valid_min"]
                    valid_max = adict["valid_max"]
                    mask = ((data < valid_min) | (data > valid_max))
                    print("""Using valid_min and valid_max attributes.
    valid_min = %10f, valid_max = %10f.""" % (valid_min, valid_max))

                elif "valid_min" in adict:
                    valid_min = adict["valid_min"]
                    mask = data < valid_min
                    print("Using the valid_min attribute; valid_min = %10f" % valid_min)

                elif "valid_max" in adict:
                    valid_max = adict["valid_max"]
                    mask = data > valid_max
                    print("Using the valid_max attribute; valid_max = %10f" % valid_max)

                elif "_FillValue" in adict:
                    fill_value = adict["_FillValue"]
                    if fill_value <= 0:
                        mask = data <= fill_value + 2 * finfo(float).eps
                    else:
                        mask = data >= fill_value - 2 * finfo(float).eps
                    print("Using the _FillValue attribute; _FillValue = %10f" % fill_value)

                elif "missing_value" in adict:
                    missing = adict["missing_value"]
                    mask = abs(data - missing) < 2 * finfo(float).eps
                    print("""Using the missing_value attribute; missing_value = %10f
    Warning: this attribute is deprecated by the NUG.""" % missing)

                else:
                    print("No missing values found. Skipping this variable...")
                    continue

                count = int(sum(mask))
                if count == 0:
                    print("No missing values found. Skipping this variable...")
                    continue
                print("Filling in %5d missing values..." % count)
                t0 = time()
                laplace(data, mask, -1, eps, initial_guess=options.initial_guess)
                var[:] = data

                # now REMOVE missing_value and _FillValue attributes
                try:
                    delattr(var, '_FillValue')
                except:
                    pass
                try:
                    delattr(var, 'missing_value')
                except:
                    pass
                print("This took %5f seconds." % (time() - t0))
            else:
                print('wrong shape')

        except Exception as message:
            print("ERROR:", message)
            print("Note: %s was not modified." % output_filename)
            exit(-1)

    print("Processing all the variables took %5f seconds." % (time() - t_zero))
    nc.close()
    try:
        move(tmp_filename, output_filename)
    except:
        print("Error moving %s to %s. Exiting..." % (tmp_filename,
                                                     output_filename))
        exit(-1)

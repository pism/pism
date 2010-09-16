#!/usr/bin/env python
# update the history attribute

from optparse import OptionParser

try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF

import numpy as np

from sys import argv,exit
from time import asctime

def copy_dim(nc1, nc2, name, direction):
    """Copy a dimension from nc1 to nc2."""
    if (name == direction):
        return
    dim1 = nc1.dimensions[name]
    dim2 = nc2.createDimension(name, len(dim1))

def copy_attributes(var1, var2):
    """Copy attributes of var1 to var2."""
    for each in var1.ncattrs():
        var2.setncattr(each, var1.getncattr(each))


def process(input, output, direction, collapse):
    """Process the file 'input', expanding or collapsing data according to
    'action' and 'direction'. Saves result in 'output'."""
    try:
        nc = CDF(input)
    except:
        print "ERROR: Can't open %s" % input
        exit(1)

    try:
        out = CDF(output, 'w', format="NETCDF3_CLASSIC")
    except:
        print "ERROR: Can't open %s" % output
        exit(1)
        
    copy_attributes(nc, out)

    for name in nc.dimensions.keys():
        copy_dim(nc, out, name, direction)

    if collapse:
        for name in nc.variables.keys():
            if name == direction:
                continue
            collapse_var(nc, out, name, direction)
        message = "Collapsed using flowline.py"
    else:
        out.createDimension(direction, 3)

        if direction == 'x':
            dim = 'y'
        else:
            dim = 'x'

        var1 = nc.variables[dim]
        delta = np.diff(var1[:])[0]

        var2 = out.createVariable(direction, 'f8', (direction,))
        var2.axis          = "%s" % direction.upper()
        var2.long_name     = "%s-coordinate in Cartesian system" % direction.upper()
        var2.standard_name = "projection_%s_coordinate" % direction
        var2.units = var1.units
        var2[:] = [-delta, 0, delta]
        
        for name in nc.variables.keys():
            expand_var(nc, out, name, direction)

    message = asctime() + ': ' + ' '.join(argv) + '\n'
    if 'history' in out.ncattrs():
      out.history = message + out.history  # prepend to history string
    else:
      out.history = message
    out.close()

def collapse_var(nc, out, name, direction):
    """Saves a collapsed (according to 'direction')
    copy of a variable 'name' in 'nc' to 'out'."""
    var1 = nc.variables[name]
    N = (len(nc.dimensions[direction]) - 1) / 2

    print "Processing %s..." % name
    dims = var1.dimensions
    if len(dims) > 1:                   # only collapse spatial fields
        dims = filter(lambda(x): x != direction, dims)
        
    var2 = out.createVariable(name, var1.dtype, dimensions=dims)
    copy_attributes(var1, var2)

    # Note: the following depends on the current PISM variable storage order: (t, z, y, x)
    if direction == 'x':
        if var1.ndim == 4:
            var2[:] = var1[:, :, :, N]  # (t, z, y, x)
        elif var1.ndim == 3:
            var2[:] = var1[:, :, N]     # (t, y, x)
        elif var1.ndim == 2:
            var2[:] = var1[:, N]        # (y, x)
        elif var1.ndim == 1:
            var2[:] = var1[:]           # (t), (y) or (x)
    elif direction == 'y':
        if var1.ndim == 4:
            var2[:] = var1[:, :, N, :]
        elif var1.ndim == 3:
            var2[:] = var1[:, N, :]
        elif var1.ndim == 2:
            var2[:] = var1[N, :]
        elif var1.ndim == 1:
            var2[:] = var1[:]

def expand_var(nc, out, name, direction):
    """Saves an expanded (according to 'direction')
    copy of a variable 'name' in 'nc' to 'out'."""
    if name == direction:
        return
    
    var1 = nc.variables[name]

    print "Processing %s..." % name

    # Copy coordinate variables and stop:
    if name in ['t', 'z', 'y', 'x', 'zb']:
        var2 = out.createVariable(name, var1.dtype, (name,))
        var2[:] = var1[:]
        copy_attributes(var1, var2)
        return

    dims = var1.dimensions
    if   len(dims) == 1:
        dims = ('y', 'x')
    elif len(dims) == 2:
        dims = ('t', 'y', 'x')
    elif len(dims) == 3:
        if name == "litho_temp":        # litho_temp is the only variable depending on 'zb'.
            dims = ('t', 'zb', 'y', 'x')
        else:
            dims = ('t', 'z', 'y', 'x')

    var2 = out.createVariable(name, var1.dtype, dims)
    copy_attributes(var1, var2)

    # Note: the following depends on the current PISM variable storage order: (t, z, y, x)
    if direction == 'x':
        if var1.ndim == 3:
            for j in range(3):
                var2[:,:,:,j] = var1[:,:,:]
        elif var1.ndim == 2:
            for j in range(3):
                var2[:,:,j] = var1[:,:]
        elif var1.ndim == 1:
            for j in range(3):
                var2[:,j] = var1[:]
    elif direction == 'y':
        if var1.ndim == 3:
            for j in range(3):
                var2[:,:,j,:] = var1[:,:,:]
        elif var1.ndim == 2:
            for j in range(3):
                var2[:,j,:] = var1[:,:]
        elif var1.ndim == 1:
            for j in range(3):
                var2[j,:] = var1[:]

parser = OptionParser()
parser.usage = "usage: %prog -o foo.nc -d {x,y} {--collapse,--expand} file.nc"
parser.description = "Collapses or expands a NetCDF file in a specified direction."

parser.add_option("-d", "--direction", dest="direction",
                  help="direction: one of x,y")
parser.add_option("-o", "--output", dest="output_filename",
                  help="output file")
parser.add_option("-e", "--expand", dest="expand", action="store_true",
                  help="expand")
parser.add_option("-c", "--collapse", dest="collapse", action="store_true",
                  help="collapse")

(opts, args) = parser.parse_args()

if len(args) != 1:
    print "ERROR: File argument is missing."
    exit(1)


if (opts.expand and opts.collapse) or ((not opts.expand) and (not opts.collapse)):
    print "ERROR: exactly one of -e and -c is required."
    exit(1)

if not opts.direction:
    if opts.collapse or (not opts.expand):
        opts.direction = 'y'
    else:
        nc = CDF(args[0])
        try:
            x = nc.variables['x']
            opts.direction = 'y'
        except:
            opts.direction = 'x'
        nc.close()
elif opts.direction not in ['x', 'y']:
    print "ERROR: Please specify direction using the -d option. (Choose one of x,y.)"
    exit(1)

if (not opts.output_filename):
    print "ERROR: Please specify the output file name using the -o option."
    exit(1)

if opts.collapse:
    print "Collapsing %s in the %s direction, writing to %s..." % (args[0], opts.direction, opts.output_filename)
else:
    print "Expanding %s in the %s direction, writing to %s..." % (args[0], opts.direction, opts.output_filename)
    
process(args[0], opts.output_filename, opts.direction, opts.collapse or (not opts.expand))

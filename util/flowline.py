#!/usr/bin/env python3

# @package flowline
# \brief This script aids in setting up and visualizing flow-line modeling runs.
#
# \details This script expands and collapses NetCDF datasets along the 'x' or
# 'y' dimension. Run flowline.py --help to see the command-line options.
#
# To set up a data-set for flow-line modeling, create a NetCDF file with an 'x'
# dimension and appropriate variables (usually x(x), acab(x), artm(x),
# bheatflx(x), thk(x) and topg(x)), then run
# \code
# flowline.py -e -o foo_pism.nc foo.nc
# \endcode
# Then bootstrap from foo_pism.nc, adding -periodicity y to PISM command-line options.
#
# To collapse a PISM flow-line model output for plotting, etc, run
# \code
# flowline.py -c -o bar.nc pism_output.nc
# \endcode

from optparse import OptionParser

import xarray as xr
import numpy as np

from sys import argv, exit
from time import asctime

def collapse(infile: str, outfile: str, direction: str):
    """
    Collapse dataset.
    """

    with xr.open_dataset(infile) as ds:

        message = asctime() + ': ' + ' '.join(argv) + '\n'
        if 'history' in ds.attrs:
            ds.attrs["history"] = message + ds.history  # prepend to history string
        else:
            ds.attrs["history"] = message
            
        ds.mean(dim=direction).to_netcdf(outfile, mode="w")

def expand(infile: str, outfile: str, direction: str):
    """
    Expand dataset.
    """

    with xr.open_dataset(infile) as ds:
        
        if direction == 'x':
            dim = 'y'
        else:
            dim = 'x'

        var1 = ds.variables[dim]
        delta = np.diff(var1[:])[0]

        var2 = [-delta, 0, delta]
        ds.expand_dims({direction: var2})
        
        ds.assign({direction: ([direction],
                               var2,
                               {"units": "m",
                                "axis": f"{direction.upper()}",
                                "standard_name": f"projection_{direction}_coordinate",
                                "long_name": f"{direction}-coordinate in projected coordinate system",
                                "_FillValue": False,
                                },
                               ),
                          }
                         )

        dims = [direction] + list(ds.sizes.keys())
        vals = [3] + list(ds.sizes.values())
        expand_vars = []
        for name in list(ds.variables.keys()):
             if name not in ['time', 'z', 'y', 'x', 'zb']:
                 expand_vars.append(name)
                 
        expand_vars_dict = {name: ds.variables[name].set_dims(dims, vals) for name in expand_vars}

        ds.update(expand_vars_dict)
        
        message = asctime() + ': ' + ' '.join(argv) + '\n'
        if 'history' in ds.attrs:
            ds.attrs["history"] = message + ds.history  # prepend to history string
        else:
            ds.attrs["history"] = message
            
        ds.to_netcdf(outfile, mode="w")
    

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
    print("ERROR: File argument is missing.")
    exit(1)


if (opts.expand and opts.collapse) or ((not opts.expand) and (not opts.collapse)):
    print("ERROR: exactly one of -e and -c is required.")
    exit(1)

if not opts.direction:
    if opts.collapse or (not opts.expand):
        opts.direction = 'y'
    else:
        ds = xr.open_dataset(args[0])
        try:
            x = ds.variables['x']
            opts.direction = 'y'
        except:
            opts.direction = 'x'
elif opts.direction not in ['x', 'y']:
    print("ERROR: Please specify direction using the -d option. (Choose one of x,y.)")
    exit(1)

if (not opts.output_filename):
    print("ERROR: Please specify the output file name using the -o option.")
    exit(1)

if opts.collapse:
    print("Collapsing %s in the %s direction, writing to %s..." % (args[0], opts.direction, opts.output_filename))
    collapse(args[0], opts.output_filename, opts.direction)
else:
    print("Expanding %s in the %s direction, writing to %s..." % (args[0], opts.direction, opts.output_filename))
    expand(args[0], opts.output_filename, opts.direction)


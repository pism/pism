#!/usr/bin/env python

from sys import argv, exit
from getopt import getopt, GetoptError
from numpy import squeeze

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

usage="""
THKMASK.PY Uses ice thickness 'thk' in one NetCDF file to mask-out ice-free areas
in identified variables in another NetCDF file.  Modifies the latter file.

Example.  To mask-out variables acab,smelt in foo.nc using 'thk' in start.nc:
  thkmask.py -H start.nc -v acab,smelt foo.nc    # 

Note:  start.nc and foo.nc must share a common grid for this to work!

Run with --help or --usage to get a this usage message.
"""

def usagefailure(message):
    print message
    print
    print usage
    exit(2)

def maskout(nc, thk, name, thktol=1.0):

    fill_value = -9.99999e+35
    valid_min = -9.0e+35

    try:
        nc.variables[name]._FillValue = fill_value
        nc.variables[name].valid_min = valid_min
    except:
        usagefailure("ERROR: VARIABLE '%s' NOT FOUND IN FILE" % name)

    dims = nc.variables[name].dimensions
    var = nc.variables[name][:]
    
    if 't' in dims:
        if var.shape[0] > 1:
            var = squeeze(var)
        for j in range(var.shape[0]):
            tmp = var[j]
            tmp[thk<thktol] = fill_value
            var[j] = tmp
    else:
        var = squeeze(var)
        var[thk<thktol] = fill_value
        
    return var


if __name__ == "__main__":
    try:
      opts, args = getopt(argv[1:], "H:v:", ["help","usage"])
    except GetoptError:
      usagefailure('ERROR: INCORRECT COMMAND LINE ARGUMENTS FOR climmask.py')
    dovars = []
    for (opt, optarg) in opts:
        if opt == "-H":
            thkfile = optarg
        if opt == "-v":
            dovars = optarg.split(",")
        if opt in ("--help", "--usage"):
            print usage
            exit(0)

    if len(args) != 1:
      usagefailure('ERROR: WRONG NUMBER OF ARGUMENTS FOR climmask.py')

    try:
      nc_thk = NC(thkfile, 'r')
    except:
      usagefailure("ERROR: FILE '%s' CANNOT BE OPENED FOR READING" % thkfile)

    try:
      nc = NC(args[0], 'a')
    except:
      usagefailure("ERROR: FILE '%s' CANNOT BE OPENED FOR MODIFYING" % args[0])

    try:
      thk = squeeze(nc_thk.variables["thk"][:])
    except:
      usagefailure("ERROR: ICE THICKNESS VARIABLE 'thk' NOT FOUND IN FILE %s" % args[0])

    for varname in dovars:
      print "  masking out variable '%s' ..." % varname
      var = nc.variables[varname]
      var[:] = maskout(nc, thk, varname)

    nc.close()
    nc_thk.close()


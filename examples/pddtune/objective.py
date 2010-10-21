#!/usr/bin/env python
from sys import argv, exit
from getopt import getopt, GetoptError

## @package objective
## Compares variables by computing an objective function.
##
## FIXME:  need help file
## 
## Run with --help to get a "usage" message.

usage="""OBJECTIVE.PY  Computes weighted L^2 norm of difference of variables.

usage: to evaluate closeness of 'acab' in foo.nc compared to 'smb' in bar.nc:
  objective.py -v acab,smb foo.nc bar.nc
"""

def usagefailure(message):
    print message
    print
    print usage
    exit(2)

if __name__ == "__main__":
    try:
      opts, args = getopt(argv[1:], "v:", ["help","usage"])
    except GetoptError:
      usagefailure('ERROR: INCORRECT COMMAND LINE ARGUMENTS FOR objective.py')
    dovars = []
    for (opt, optarg) in opts:
        if opt == "-v":
            dovars = optarg.split(",")
        if opt in ("--help", "--usage"):
            print usage
            exit(0)

    if len(args) != 3:
      usagefailure('ERROR: objective.py requires exactly three file names')

    if len(dovars) != 2:
      usagefailure('ERROR: -v option for objective.py requires exactly two variable names')

    try:
      from netCDF4 import Dataset as NC
    except:
      from netCDF3 import Dataset as NC

    try:
      nc_pism = NC(args[0], 'r')
    except:
      usagefailure("ERROR: FILE '%s' CANNOT BE OPENED FOR READING" % args[0])
    try:
      nc_reference = NC(args[1], 'r')
    except:
      usagefailure("ERROR: FILE '%s' CANNOT BE OPENED FOR READING" % args[1])

    from numpy import squeeze, shape, double

    print "  comparing variable '%s' in %s to '%s' in %s ..." % \
        (dovars[0], args[0], dovars[1], args[1])

    pism_var = squeeze(nc_pism.variables[dovars[0]])
    ref_var = squeeze(nc_reference.variables[dovars[1]])
    Mx, My = shape(pism_var)
    valid_min = nc_pism.variables[dovars[0]].valid_min
    
    sumdiff = 0.0
    sumL2 = 0.0
    countvalid = 0
    for i in range(Mx):
      for j in range(My):
        if pism_var[i,j] >= valid_min:
          countvalid += 1
          diff_vars = pism_var[i,j] - ref_var[i,j]
          sumdiff += diff_vars
          sumL2 += diff_vars * diff_vars
    avdiff = sumdiff / double(countvalid)
    avL2 = sumL2 / double(countvalid)
    
    nc_pism.close()
    nc_reference.close()

    print "  %d locations for valid comparison:" % countvalid
    print "    average of signed differences is  %14.5f" % avdiff
    print "    average of squared differences is %14.5f" % avL2

    print "  writing these values to text file %s ..." % args[2]
    diffs = file(args[2],'a')
    diffs.write("%s %14.5f %14.5f\n" % (args[0], avdiff, avL2))
    diffs.close()



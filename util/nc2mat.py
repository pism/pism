#!/usr/bin/env python

## @package nc2mat
## Reads specified variables from a NetCDF file and writes them to an output file in the MATLAB binary data file format \c .mat, supported by MATLAB version 5 and later). It depends on \c netcdf4-python and \c SciPy.
##
## The options go before the input file name.  Here are some examples:
## \verbatim nc2mat.py ross.nc \endverbatim will read all the variables from \c ross.nc and write them to \c ross.mat, \verbatim nc2mat.py -v thk,acab -o ross_matlab.mat ross.nc \endverbatim will read variables \c thk and \c acab from \c ross.nc and write them to \c ross_matlab.nc, and
## \verbatim nc2mat.py -e x,y,z ross.nc \endverbatim
## will read all the variables except \c x, \c y and \c z from \c ross.nc and write them to \c ross.mat.
## Because PISM saves all diagnostic variables using the single precision (type \c float), while most MATLAB plotting functions expect data to have double precision (type \c double), so one might need to convert them explicitly.  For example, assuming that \c ross.nc contains data produced by a PISM run, running
## \verbatim nc2mat.py -v cbar ross.nc \endverbatim and then typing
## \code
## >> load ross.mat;
## >> cbar = double(reshape(cbar, 147, 147));
## >> pcolor(cbar); colorbar();
## \endcode
## in the MATLAB shell will produce a picture of the Ross Ice Shelf shaded by modeled ice speed.

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
    
from getopt import getopt, GetoptError
from os.path import splitext
from sys import argv

def print_options_and_exit(exit_code=0):
    print """nc2mat.py converts a NetCDF file to a MATLAB .mat binary data file
It uses MATLAB v5 format.

        nc2mat.py [-o<output file name>] [-v<variable list>] [-e<list of variables to exclude>] <file>
Options:
        -o<file name> or --output=<file name>
            specifies the output file name
            (omit to use the original file name; "nc2mat.py foo.nc" will produce foo.mat) 
        -v<comma-separated list of variables> or --variables=<comma-separated list of variables>
            specifies the list of variables to write to the output file
            (omit to dump all the variables)
        -e<comma-separated list of variables> or --exclude=<comma-separated list of variables>
            specifies the list of variables to exclude

        Note: options should go before the input file name.
Examples:
        nc2mat.py -e x,y,z ross.nc
            will write all the variables from ross.nc except x,y,z to ross.mat
        nc2mat.py -v thk,acab -o ross_matlab.mat ross.nc 
            will write thk and acab to ross_matlab.mat
        nc2mat.py ross.nc
            will write all the variables from ross.nc to ross.mat
"""
    exit(exit_code)

# process command-line options:
try:
    opts, args = getopt(argv[1:], "o:v:e:", ["output=", "variables=", "exclude="])
    # defaults:
    variables = []    # all the variables
    exclude_list = [] # do not exclude anything
    output_name = ""  # same as input_name, but with .mat instead of .nc
    for opt, arg in opts:
        if opt in ("-o", "--output"):
            output_name = arg
        if opt in ("-v", "--variables"):
            variables = arg.split(',')
        if opt in ("-e", "--exclude"):
            exclude = arg.split(',')
    if len(args) != 1:
        print "ERROR! Please specify one input file. Exiting..."
        exit(-1)
    input_name = args[0]
    if output_name == "":
        output_name = splitext(input_name)[0] + '.mat'
except GetoptError:
    print_options_and_exit(-1)

try:
    input_file = NC(input_name, "r")
except Exception, message:
    print "ERROR: ", message
    print "Can't open %s. Exiting..." %(input_name)
    exit(-1)

print "Reading data from %s..." % (input_name)

# if no variables were specified, assume that we need to dump them all
if variables == []:
    variables = input_file.variables.keys()

print "Loading variables...",
mdict = {}
missing = []
for each in variables:
    try:
        if each in exclude_list:
            continue
        mdict[each] = input_file.variables[each][:]
        print each,
    except Exception:
        missing.append(each)
        continue
print ""

if len(missing) > 0:
    print "Warning! the following variable(s) was(were) not found in '%s':" % (input_name), missing

print "Excluded variables: ", exclude_list

if mdict.keys() == []:
    print "WARNING:  No data to write. Exiting ..."
    print_options_and_exit(-2)

print "Writing data to %s..." % (output_name)
try:
    from scipy.io.matlab.mio import savemat
except:
    print "ERROR! Can't import 'savemat' from scipy.io.matlab.mio.  Exiting ..."
    exit(-3)
from scipy.io.matlab.mio import savemat
try:
    savemat(output_name, mdict, appendmat=False, format='5')
except:
    print "ERROR! Can't write to %s. Exiting ..." % (output_name)
    exit(-4)

#!/usr/bin/env python

from scipy.io.matlab.mio import savemat
from sys import argv, exit
from netCDF3 import Dataset as NC
from getopt import getopt, GetoptError
from os.path import splitext

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
    print "No data to write. Exiting..."
    exit(-1)

print "Writing data to %s..." % (output_name)
try:
    savemat(output_name, mdict, appendmat=False, format='5')
except:
    print "ERROR! Can't write to %s. Exiting..." % (output_name)
    exit(-1)

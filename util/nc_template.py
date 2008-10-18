#!/usr/bin/env python
from netCDF3 import *
from getopt import getopt, GetoptError
from subprocess import Popen, PIPE
from sys import argv, exit, stdout
import os

# Reads pism_state.cdl, replaces magic numbers with actual values and runs
# ncgen on the result to produce a .nc file having the right structure.

def make_template(x, y, vars, pism_state_cdl, output_filename, remove=False):
    """Reads 'pism_state_cdl', replaces magic numbers with given grid dimensions,
    filters out all the variables except ones specified in 'vars' and runs
    ncgen to create 'output_filename', a NetCDF file having the structure of a
    PISM input file.

    If remove == True, then 'output_filename' is removed first.
    
    If vars == [], then includes x,y,z,zb,t plus all the variables
       having the 'pism_intent' attribute of 'mapping', 'pism_state',
       'climate_steady' or 'climate_snapshot'.

    All variables listed in 'vars' but not defined in pism_state_cdl are
    ignored silently.
    """
    def replace_magic_numbers(line):
        pairs = [("91", str(x)),
                 ("92", str(y)),
                 ("93", "1"),
                 ("94", "1")]
        for each in pairs:
            line = line.replace(each[0], each[1])
        return line

    def read_variable(j, printp = False):
        """Reads until the next variable definition is found or the end of the
        'variables' section is reached. If printp = True, then also prints the
        lines read."""
        from numpy import any
        for j in range(j, N):
            positions = map(lines[j].strip().find,
                            ["int ", "double ",
                             "float ", "byte ",
                             "// global attributes"])
            condition = map(lambda(x): x == 0, positions)
            if any(condition):
                break
            if printp:
                ncprint(lines[j])
        return j

    def check_variables(line):
        """Checks if 'line' is the beginning of a variable definition for all variables
        in 'vars'."""
        for v in vars:
            s = " " + v + "("
            if line.strip().find(s) >= 0:
                return True
            s = " " + v + ";"
            if line.strip().find(s) >= 0:
                return True
        return False
        
    try:
        if remove:
            os.remove(output_filename)
    except:
        print "Warning! Could not remove %s." % output_filename

    try:
        pism_cdl = open(pism_state_cdl, 'r')
    except:
        print "Error! Could not open %s. Exiting..." % pism_state_cdl
        exit(-1)

    try:
        ncgen = Popen(["ncgen", "-o%s" % (output_filename), "-"], stdin=PIPE)
        ncprint = ncgen.stdin.write
#        ncprint = stdout.write
    except:
        print "Error! Could not run ncgen. Exiting..."
        exit(-1)

    lines = pism_cdl.readlines()
    pism_cdl.close()

    # if vars == [], then we need to find all the mapping, model_state and
    # climate-related variables:
    if vars == []:
        for line in lines:
            if line.find("//") == 0:
                continue
            if line.find("pism_intent") >= 0:
                parts = line.split(":")
                if (parts[1].find("mapping") >= 0 or
                    parts[1].find("model_state") >= 0 or
                    parts[1].find("climate") >= 0):
                    name = parts[0].strip()
                    if name not in ("age", "temp", "litho_temp"):
                        vars.append(name)
    vars += ["x", "y", "z", "zb", "t", "polar_stereographic"]

    N = len(lines); j = 0
    # Skip all the comments, then replace the header:
    fname = os.path.basename(output_filename)
    netcdf_name = os.path.splitext(fname)[0]
    for j in range(0, N):
        if lines[j].find("netcdf pism_state") >= 0:
            ncprint("netcdf %s {\n" % netcdf_name)
            j += 1
            break

    # Replace magic numbers:
    for j in range(j, N):
        if lines[j].find("variables:") >= 0:
            ncprint(lines[j])
            j += 1
            break
        ncprint(replace_magic_numbers(lines[j]))

    # Process variables:
    while j < N:
        if lines[j].find("global attributes") >= 0:
            break
        if check_variables(lines[j]):
            ncprint(lines[j])
            j = read_variable(j + 1, printp = True)
        else:
            j = read_variable(j + 1, printp = False)

    # Feed the rest (i.e. global attributes) to ncgen
    for j in range(j,N):
        ncprint(lines[j])

    # Waiting for ncgen...
    ncgen.stdin.close()
    ncgen.wait()

    # Do some basic initialization:
    try:
        nc = Dataset(output_filename, 'a')
    except:
        print "Could not open %s. Exiting..." % output_filename
        exit(-1)

    dummy_variables = ["t", "z", "zb"]
    for name in dummy_variables:
        var = nc.variables[name]
        var[0] = 0.0
    nc.close()

if __name__ == "__main__":
    try:
        opts, args = getopt(argv[1:], "x:y:v:p:o:r", ["x=", "y=", "variables="
                                                      "pism-state=","output=",
                                                      "remove"])
        x = 0;y = 0
        vars = []
        pism_state_cdl = "pism_state.cdl"
        output_filename = ""
        for opt, arg in opts:
            if opt in ("-p", "--pism-state"):
                pism_state_cdl = arg
            if opt in ("-o", "--output"):
                output_filename = arg
            if opt in ("-x", "--x"):
                x = int(arg)
            if opt in ("-y", "--y"):
                y = int(arg)
            if opt in ("-v", "--variables"):
                vars = arg.split(",")
        if (x <= 0) or (y <= 0):
            print "Both x and y should be greater than zero. Exiting..."
            exit(-1)
        if output_filename == "":
            print "Please specify the output file name. Exiting..."
            exit(-1)
    except GetoptError:
        print """Options:
  -p [--pism-state=] location of pism_state.cdl
                     (nc_template.py loo
  -o [--output=]     output file name
  -x [--x=]          grid width (points)
  -y [--y=]          grid height (points)
  -v [--variables=]  comma-separated list of variable names to include
                     if this list is omitted, then all the mapping, climate and 2D
                     model_state variables are included"""
        exit(-1)
    make_template(x, y, vars, pism_state_cdl, output_filename, False)

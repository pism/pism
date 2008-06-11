#!/usr/bin/env python
from getopt import getopt, GetoptError
from subprocess import Popen, PIPE
from sys import argv, exit
import os

# Read pism_state.cdl, replace magic numbers with actual values and run ncgen
# on the result to produce a .nc file having the right structure.

def make_template(x, y, pism_state_cdl, output_filename, remove=False):
    """Reads 'pism_state_cdl', replaces magic numbers with given grid dimensions
    and runs ncgen to create 'output_filename', a NetCDF file having the
    structure of a PISM input file.

    If remove == True, then 'output_filename' is removed first.
    """
    try:
        if remove: os.remove(output_filename)
    except:
        print "Warning! Could not remove %s." % output_filename

    pism_cdl = open(pism_state_cdl, 'r')
    def replace_magic_numbers(line):
        pairs = [("pism_state", os.path.splitext(output_filename)[0]),
                 ("91", str(x)),
                 ("92", str(y)),
                 ("93", "1"),
                 ("94", "1")]
        for each in pairs:
            line = line.replace(each[0], each[1])
        return line

    ncgen = Popen(["ncgen", "-o%s" % (output_filename), "-"], stdin=PIPE)

    for line in pism_cdl:
        ncgen.stdin.write(replace_magic_numbers(line))
    pism_cdl.close()
    ncgen.stdin.close()

    # waiting for ncgen...
    ncgen.wait()

if __name__ == "__main__":
    opts, args = getopt(argv[1:], "x:y:p:o:r", ["x=", "y=", 
                                                "pism-state=","output=",
                                                "remove="])
    x = 0;y = 0
    pism_state_cdl = "pism_state.cdl"
    output_filename = ""
    remove = False
    for opt, arg in opts:
        if opt in ("-p", "--pism-state"):
            pism_state_cdl = arg
        if opt in ("-o", "--output"):
            output_filename = arg
        if opt in ("-x", "--x"):
            x = int(arg)
        if opt in ("-y", "--y"):
            y = int(arg)
        if opt in ("-r", "--remove"):
            remove = True
    if (x <= 0) or (y <= 0):
        print "Both x and y should be greater than zero. Exiting..."
        exit(-1)
    if output_filename == "":
        print "Please specify the output file name. Exiting..."
        exit(-1)

    make_template(x, y, pism_state_cdl, output_filename, remove)

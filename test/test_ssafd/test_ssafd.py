#!/usr/bin/env python

from PISMNC import PISMDataset as NC
from optparse import OptionParser
import numpy as np
import subprocess
import shlex
import sys

parser = OptionParser()

parser.usage = "%prog [options]"
parser.description = "Test the SSAFD solver using various geometric configurations."
parser.add_option("-o", dest="output_file_prefix", default="ssafd_test",
                  help="output file prefix")
parser.add_option("-L", dest="L", type=float, default=10.0,
                  help="horizontal domain dimensions, km")
parser.add_option("-H", dest="H", type=float, default=500.0,
                  help="ice thickness in icy areas")

(options, args) = parser.parse_args()

M = 3


def generate_input(N):
    """Generate a PISM bootstrapping file for a neighborhood number N."""
    output_filename = options.output_file_prefix + "_%08d.nc" % N
    pism_output_filename = options.output_file_prefix + "_o_%08d.nc" % N
    nc = NC(output_filename, 'w')

    format_string = "{0:0%db}" % (M * M)
    string = format_string.format(N)

    x = np.linspace(-options.L * 1000, options.L * 1000, M)

    zeros = np.zeros((M, M))
    thk = np.zeros(M * M)

    topg = zeros.copy() - 0.1 * options.H
    topg[1, 1] = -options.H

    for j in range(M * M):
        if string[j] == '1':
            thk[j] = options.H

    nc.create_dimensions(x, x)
    nc.write("topg", topg, attrs={"units": "m", "long_name": "bed_topography"})
    nc.write("climatic_mass_balance", zeros, attrs={"units": "kg m-2 year-1"})
    nc.write("ice_surface_temp", zeros, attrs={"units": "Celsius"})
    nc.write("thk", thk.reshape(M, M),
             attrs={"units": "m", "long_name": "land_ice_thickness"})

    nc.close()

    return output_filename, pism_output_filename


def run_pismr(input_filename, output_filename):
    command = "pismr -i %s -bootstrap -o %s -Mx %d -My %d -Lz 1000 -Mz 5 -stress_balance ssa+sia -cfbc -y 0.001 -verbose 1" % (
        input_filename, output_filename, M, M)
    print("Running %s" % command)
    subprocess.call(shlex.split(command))


for k in range(2 ** (M * M)):
    input_file, output_file = generate_input(k)
    run_pismr(input_file, output_file)

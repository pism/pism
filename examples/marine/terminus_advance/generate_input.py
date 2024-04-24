#!/usr/bin/env python3

from PISMNC import PISMDataset as NC
from optparse import OptionParser
import numpy as np

parser = OptionParser()

parser.usage = "%prog [options]"
parser.description = "Fills missing values in variables selected using -v in the file given by -f."
parser.add_option("-o", dest="output_file_name", default="advance_test.nc",
                  help="output file name")
parser.add_option("-L", dest="length", default=100e3, type=float,
                  help="domain length, meters")
parser.add_option("-M", dest="Mx", default=101, type=int,
                  help="number of grid points in the flow direction")

(options, args) = parser.parse_args()

My = 3
bed_slope = 0.02
ice_thickness = 1000.0
U = 100.0

nc = NC(options.output_file_name, 'w')

x = np.linspace(0, options.length, options.Mx)
y = x[0:3] - x[1]

xx, yy = np.meshgrid(x, y)

zeros = np.zeros((My, options.Mx))
topg = (0.5 * options.length - xx) * bed_slope

thk = zeros.copy()
thk[topg > 200] = ice_thickness

v = zeros.copy()
u = zeros.copy() + U

nc.create_dimensions(x, y)
nc.write("topg", topg, attrs={"units": "m", "long_name": "bed_topography"})
nc.write("climatic_mass_balance", zeros, attrs={"units": "kg m-2 year-1"})
nc.write("ice_surface_temp", zeros, attrs={"units": "degree_Celsius"})
nc.write("thk", thk,
         attrs={"units": "m", "standard_name": "land_ice_thickness"})
nc.write("ubar", u,
         attrs={"units": "m/year", "long_name": "x-component of velocity"})
nc.write("vbar", v,
         attrs={"units": "m/year", "long_name": "y-component of velocity"})

nc.close()

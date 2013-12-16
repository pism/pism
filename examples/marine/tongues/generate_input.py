#!/usr/bin/env python

from PISMNC import PISMDataset as NC
import numpy as np

import argparse

parser = argparse.ArgumentParser()
parser.description = "Generates an input file for the 'tongues' experiment."

parser.add_argument("-M", dest="M", type=int,
                    help="grid size", default=101)
parser.add_argument("-L", dest="L", type=float,
                    help="domain size, meters", default=1e5)
parser.add_argument("-o", dest="output", help="output file name",
                    default="tongues.nc")
options = parser.parse_args()
L = options.L
M = options.M

x = np.linspace(-L, L, M)
y = np.linspace(-L, L, M)

dx = x[1] - x[0]

xx,yy = np.meshgrid(x, y)

z = np.zeros_like(xx) - 1000.0

def tongue(xx, x0, width):
    result = np.zeros_like(xx)
    result[:,x0:x0+width] = 100.0
    return result

thk = np.zeros_like(z)

x0 = 3
width = 1
spacing = 5
while x0 + width < M - 1:
    thk += tongue(xx, x0, width)
    x0 += width + spacing
    width += 1

thk[5:,:] = 0

bcflag = np.zeros_like(thk)
bcflag[thk > 0] = 1

z[thk > 0] = -(910.0 / 1028.0) * 100.0 + 1

ubar = np.zeros_like(thk)

vbar = np.zeros_like(thk)
vbar[bcflag == 1] = 100.0

try:
    nc = NC(options.output, 'w')
except:
    nc = NC(options.output, 'a')

try:
    nc.create_dimensions(x, y, time_dependent=False)

    nc.define_2d_field("topg", attrs={"units" : "m",
                                      "long_name" : "bedrock topography"})
    nc.define_2d_field("thk", attrs={"units" : "m",
                                     "long_name" : "ice thickness"})

    nc.define_2d_field("climatic_mass_balance", attrs={"units" : "m/year"})
    nc.define_2d_field("ice_surface_temp", attrs={"units" : "Kelvin"})

    nc.define_2d_field("u_ssa_bc", attrs={"units" : "m/year"})
    nc.define_2d_field("v_ssa_bc", attrs={"units" : "m/year"})
except:
    pass

nc.write("topg", z)
nc.write("thk", thk)
nc.write("climatic_mass_balance", np.zeros_like(xx))
nc.write("ice_surface_temp", np.zeros_like(xx) + 273.15 - 30.0)
nc.write("u_ssa_bc", ubar)
nc.write("v_ssa_bc", vbar)
nc.write("bcflag", bcflag)

nc.close()

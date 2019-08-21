#!/usr/bin/env python
from PISMNC import PISMDataset as NC
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.description = "Generates an input file for the 'tongues' experiment."

parser.add_argument("-M", dest="M", type=int, help="grid size", default=101)
parser.add_argument("-L", dest="L", type=float, help="domain size, meters", default=1e5)
parser.add_argument("-o", dest="output", help="output file name", default="tongues.nc")
options = parser.parse_args()

L = options.L
M = options.M

# grid
x = np.linspace(-L, L, M)
y = np.linspace(-L, L, M)
xx, yy = np.meshgrid(x, y)


def tongue(xx, x0, width):
    "create one ice tongue"
    result = np.zeros_like(xx)
    result[:, x0:x0 + width] = 100.0
    return result


thk = np.zeros_like(xx)

x0 = 3
width = 1
spacing = 5
while x0 + width < M - 1:
    thk += tongue(xx, x0, width)
    x0 += width + spacing
    width += 1

# make tongues shorter
thk[5:, :] = 0

bc_mask = np.zeros_like(thk)
bc_mask[thk > 0] = 1

# make the bed deep everywhere except in icy areas, where it is barely
# grounded
z = np.zeros_like(xx) - 1000.0
z[thk > 0] = -(910.0 / 1028.0) * 100.0 + 1

# Velocity Dirichlet B.C.:
ubar = np.zeros_like(thk)
vbar = np.zeros_like(thk)
vbar[bc_mask == 1] = 100.0

try:
    nc = NC(options.output, 'w')

    nc.create_dimensions(x, y, time_dependent=False)

    nc.define_2d_field("topg", attrs={"units": "m",
                                      "long_name": "bedrock topography"})
    nc.define_2d_field("thk", attrs={"units": "m",
                                     "long_name": "ice thickness"})

    nc.define_2d_field("climatic_mass_balance", attrs={"units": "kg m-2 year-1"})
    nc.define_2d_field("ice_surface_temp", attrs={"units": "Celsius"})

    nc.define_2d_field("u_ssa_bc", attrs={"units": "m/year"})
    nc.define_2d_field("v_ssa_bc", attrs={"units": "m/year"})
except:
    nc = NC(options.output, 'a')

nc.write("topg", z)
nc.write("thk", thk)
nc.write("climatic_mass_balance", np.zeros_like(xx))
nc.write("ice_surface_temp", np.zeros_like(xx) - 30.0)  # irrelevant
nc.write("u_ssa_bc", ubar)
nc.write("v_ssa_bc", vbar)
nc.write("bc_mask", bc_mask)

nc.close()

#!/usr/bin/env python3

import sys
import os
import subprocess
import shlex
import numpy as np
import netCDF4

def run(cmd):
    sys.stderr.write(cmd + "\n")
    if subprocess.call(shlex.split(cmd)) != 0:
        sys.exit(1)

pism_path = sys.argv[1]
mpiexec = sys.argv[2]

cmd = "{path}/pismv -test F -y 10 -verbose 1 -o_size small -no_report -o in-temp-continuity.nc".format(path=pism_path)

run(cmd)

deltas = []
dts = [100, 50]
for dt in dts:
    try:
        cmd = "{path}/pismr -eisII B -y 2400 -Mx 16 -My 16 -Mz 21 -Lbz 1000 -Mbz 11 -energy enthalpy -regrid_file in-temp-continuity.nc -regrid_vars thk -verbose 1 -max_dt {dt} -o out-temp-continuity.nc -output.sizes.medium temp -gradient mahaffy".format(path=pism_path, dt=dt)

        run(cmd)

        with netCDF4.Dataset("out-temp-continuity.nc") as f:
            # note: this file stores 3D variables in the time,y,x,(z|zb) order
            temp = f.variables["temp"][0, :, :, 0]              # pick the bottom layer
            litho_temp = f.variables["litho_temp"][0, :, :, -1] # pick the top layer
            deltas.append(np.abs(temp - litho_temp).max())
    finally:
        os.remove("out-temp-continuity.nc")

# these deltas are observed to decrease O(dt^1) approximately, which is expected from theory
for (dt, delta) in zip(dts, deltas):
    sys.stderr.write("dt = %f, delta = %f\n" % (dt, delta))

# the only test is whether they decrease; no rate measured
if any(np.diff(deltas) >= 0):
    print(np.diff(deltas))
    sys.exit(1)

os.remove("in-temp-continuity.nc")

sys.exit(0)

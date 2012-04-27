#!/usr/bin/env python

from sys import exit, argv, stderr
from os import system
from numpy import squeeze, abs, diff

try:
    from netCDF3 import Dataset as NC
except:
    from netCDF4 import Dataset as NC

pism_path=argv[1]
mpiexec=argv[2]

stderr.write("Testing: temperature continuity at ice-bed interface (polythermal case).\n")

cmd = "%s %s/pismv -test F -y 10 -verbose 1 -o bar.nc" % (mpiexec, pism_path)
stderr.write(cmd + '\n')

e = system(cmd)
if e != 0:
  exit(1)

deltas = []
dts = [200, 100]
for dt in dts:
    cmd = "%s %s/pisms -eisII B -y 5000 -Mx 16 -My 16 -Mz 21 -Lbz 1000 -Mbz 11 -no_cold -regrid_file bar.nc -regrid_vars thk -verbose 1 -max_dt %f -o foo.nc -o_size big" % (mpiexec, pism_path, dt)
    stderr.write(cmd + '\n')

    e = system(cmd)
    if e != 0:
        exit(1)

    e = system("ncks -O -v temp -d z,0 foo.nc temp.nc")
    if e != 0:
        exit(1)

    e = system("ncks -O -v litho_temp -d zb,10 foo.nc litho_temp.nc")
    if e != 0:
        exit(1)

    nc1 = NC("temp.nc")
    nc2 = NC("litho_temp.nc")

    temp = squeeze(nc1.variables['temp'][:])
    litho_temp = squeeze(nc2.variables['litho_temp'][:])

    deltas.append(abs(temp - litho_temp).max())

# these deltas are observed to decrease O(dt^1) approximately, which is expected from theory
for (dt, delta) in zip(dts, deltas):
    stderr.write("dt = %f, delta = %f\n" % (dt, delta))

# the only test is whether they decrease; no rate measured
if any(diff(deltas) > 0):
    exit(1)

system("rm foo.nc foo.nc~ bar.nc temp.nc litho_temp.nc")
exit(0)


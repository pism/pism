#!/usr/bin/env python

from sys import exit, argv, stderr
from os import system
from numpy import squeeze, abs, diff

try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

pism_path = argv[1]
mpiexec = argv[2]

stderr.write("Testing: temperature continuity at ice-bed interface (polythermal case).\n")

cmd = "%s %s/pismv -test F -y 10 -verbose 1 -o bar-temp-continuity.nc" % (mpiexec, pism_path)
stderr.write(cmd + '\n')

e = system(cmd)
if e != 0:
    exit(1)

deltas = []
dts = [200, 100]                # FIXME: this is fragile and the test fails if I add smaller dt like 50 here
for dt in dts:
    cmd = "%s %s/pisms -eisII B -y 2400 -Mx 16 -My 16 -Mz 21 -Lbz 1000 -Mbz 11 -energy enthalpy -regrid_file bar-temp-continuity.nc -regrid_vars thk -verbose 1 -max_dt %f -o foo-temp-continuity.nc -o_size big" % (
        mpiexec, pism_path, dt)
    stderr.write(cmd + '\n')

    e = system(cmd)
    if e != 0:
        exit(1)

    e = system("ncks -O -v temp -d z,0 foo-temp-continuity.nc temp-temp-continuity.nc")
    if e != 0:
        exit(1)

    e = system("ncks -O -v litho_temp -d zb,10 foo-temp-continuity.nc litho_temp-temp-continuity.nc")
    if e != 0:
        exit(1)

    nc1 = NC("temp-temp-continuity.nc")
    nc2 = NC("litho_temp-temp-continuity.nc")

    temp = squeeze(nc1.variables['temp'][:])
    litho_temp = squeeze(nc2.variables['litho_temp'][:])

    deltas.append(abs(temp - litho_temp).max())

# these deltas are observed to decrease O(dt^1) approximately, which is expected from theory
for (dt, delta) in zip(dts, deltas):
    stderr.write("dt = %f, delta = %f\n" % (dt, delta))

# the only test is whether they decrease; no rate measured
if any(diff(deltas) > 0):
    print(diff(deltas))
    exit(1)

system("rm foo-temp-continuity.nc foo-temp-continuity.nc~ bar-temp-continuity.nc temp-temp-continuity.nc litho_temp-temp-continuity.nc")
exit(0)

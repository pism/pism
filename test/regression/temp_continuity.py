#!/usr/bin/env python

from sys import exit, argv
from os import system
from numpy import squeeze, abs, diff

try:
    from netCDF3 import Dataset as NC
except:
    from netCDF4 import Dataset as NC

pism_path=argv[1]
mpiexec=argv[2]

print "Testing: temperature continuity at ice-bed interface (cold case)."

deltas = []
dts = [2, 1, 0.5]
for dt in dts:
    e = system("%s %s/pisms -eisII A -Mx 31 -My 31 -Mz 31 -Lbz 1000 -Mbz 11 -y 100 -o foo.nc -max_dt %f" %
               (mpiexec, pism_path, dt))
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

for (dt, delta) in zip(dts, deltas):
    print "dt = %f, delta = %f" % (dt, delta)

diff = diff(deltas)

if any(diff > 0):
    exit(1)

system("rm foo.nc foo.nc~ temp.nc litho_temp.nc")
exit(0)

#!/bin/bash

#set -v

MPIDO="mpiexec -n 4"
PATHTO=""

if [ "" ]; then      # zero length string is false
#if [ "pre" ]; then  # <-- do this to create "start.n" for comparison runs
  $MPIDO ${PATHTO}pisms -Mx 121 -My 121 -Mz 101 -y 10 -o pre.nc
  $MPIDO ${PATHTO}pismr -i pre.nc -ye 8000 -o start.nc
fi

# only goal of this test script is to *LOOK* at the results; verification is another matter

$MPIDO ${PATHTO}pismr -i start.nc -hydrology tillcan -ys 0 -y 20 -o tillcan.nc
$MPIDO ${PATHTO}pismr -boot_file start.nc -hydrology tillcan -Mx 181 -My 181 -Lz 4100 -Mz 101 -ys 0 -y 2 -regrid_file tillcan.nc -regrid_vars enthalpy,lithotemp,bwat,bmelt -o bootandregridtillcan.nc

$MPIDO ${PATHTO}pismr -i start.nc -hydrology diffuseonly -report_mass_accounting -ys 0 -y 20 -o diffuse.nc

#$MPIDO ${PATHTO}pismr -i start.nc -hydrology lakes       -ys 0 -y 20 -o lakes.nc
$MPIDO ${PATHTO}pismr -i start.nc -hydrology lakes -report_mass_accounting -hydrology_use_const_bmelt -hydrology_const_bmelt 3.1689e-10 -ys 0 -y 20 -extra_file extras_lakes.nc -extra_times 0:0.5:20 -extra_vars bwat,bwp,bwatvel,thk,enwat -o lakes.nc

# see for basic tests of -hydrology distributed:
#    test/regression/test_29.py
#    examples/nbreen/run.sh

#set +v


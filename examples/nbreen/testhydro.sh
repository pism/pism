#!/bin/bash

#set -v

MPIDO="mpiexec -n 4"
PATHTO=""

if [ "" ]; then      # zero length string is false
#if [ "pre" ]; then  # <-- do this to create "start.nc" for comparison runs
  $MPIDO ${PATHTO}pisms -Mx 121 -My 121 -Mz 101 -y 10 -o pre.nc
  $MPIDO ${PATHTO}pismr -i pre.nc -ye 8000 -o start.nc
fi

# only goal of this test script is to *LOOK* at the results; verification is another matter

$MPIDO ${PATHTO}pismr -i start.nc -hydrology null -ys 0 -y 20 -o null.nc
$MPIDO ${PATHTO}pismr -boot_file start.nc -hydrology null -Mx 181 -My 181 -Lz 4100 -Mz 101 -ys 0 -y 2 -regrid_file null.nc -regrid_vars enthalpy,lithotemp,bwat,bmelt -o bootandregridnull.nc

cp start.nc start_withbwat.nc
ncrename -v tillwat,bwat start_withbwat.nc
ncks -A -v tillwat start.nc start_withbwat.nc

$MPIDO ${PATHTO}pismr -i start_withbwat.nc -hydrology routing -report_mass_accounting -hydrology_use_const_bmelt -hydrology_const_bmelt 3.1689e-10 -ys 0 -y 20 -extra_file extras_routing.nc -extra_times 0:0.5:20 -extra_vars bwat,bwp,bwatvel,thk,enwat -o routing.nc

# see for basic tests of -hydrology distributed:
#    test/regression/test_29.py
#    examples/nbreen/run.sh

#set +v


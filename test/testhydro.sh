#!/bin/bash

#set -v

MPIDO="mpiexec -n 4"
PATHTO="./"

if [ "" ]; then      # zero length string is false
#if [ "pre" ]; then  # <-- do this to create "start.n" for comparison runs
  $MPIDO ${PATHTO}pisms -Mx 121 -My 121 -Mz 101 -y 10 -o pre.nc
  $MPIDO ${PATHTO}pismr -i pre.nc -ye 8000 -o start.nc
fi

# only goal of this test script is to *LOOK* at the results; verification is another matter

#$MPIDO ${PATHTO}pismr -i start.nc                        -y 20 -o tillcandefault.nc

$MPIDO ${PATHTO}pismr -i start.nc -hydrology tillcan     -y 20 -o tillcan.nc
#$MPIDO ${PATHTO}pismr -i tillcan.nc -hydrology tillcan   -y 0.001 -o continuetillcan.nc
#$MPIDO ${PATHTO}pismr -boot_file tillcan.nc -hydrology tillcan -Mx 181 -My 181 -Lz 4100 -Mz 101 -y 0.001 -o boottillcan.nc
#$MPIDO ${PATHTO}pismr -boot_file tillcan.nc -hydrology tillcan -Mx 181 -My 181 -Lz 4100 -Mz 101 -y 0.001 -regrid_file continuetillcan.nc -regrid_vars enthalpy,lithotemp,bwat,bmelt -o bootandregridtillcan.nc

$MPIDO ${PATHTO}pismr -i start.nc -hydrology diffuseonly -y 20 -o diffuse.nc

$MPIDO ${PATHTO}pismr -i start.nc -hydrology lakes       -y 20 -o lakes.nc
#$MPIDO ${PATHTO}pismr -i lakes.nc -hydrology lakes       -y 0.001 -o continuelakes.nc
$MPIDO ${PATHTO}pismr -i start.nc -hydrology lakes -hydrology_use_const_bmelt -hydrology_const_bmelt 3.1689e-10 -y 20 -o lakes_const.nc

#$MPIDO ${PATHTO}pismr -i start.nc -hydrology distributed -y 20 -o distributed.nc

#set +v


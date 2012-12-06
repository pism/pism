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
$MPIDO ${PATHTO}pismr -i start.nc -hydrology diffuseonly -y 20 -o diffuse.nc
$MPIDO ${PATHTO}pismr -i start.nc -hydrology lakes       -y 20 -o lakes.nc
#$MPIDO ${PATHTO}pismr -i start.nc -hydrology distributed -y 20 -o distributed.nc

#set +v


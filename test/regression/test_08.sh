#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test # 8: regridding: coarse -> fine -> coarse (vertical direction)."
# The list of files to delete when done.
files="coarse1.nc coarse2.nc fine1.nc fine2.nc"

OPTS="-y 0"

set -e -x

# Create a file to regrid from:
$PISM_PATH/pisms -no_cold -Mx 11 -My 11 -Mz 11 -y 6000 -max_dt 300.0 -o coarse1.nc 
# Create another file with a finer grid:
$PISM_PATH/pisms -no_cold -Mx 11 -My 11 -Mz 21 -y 6000 -max_dt 300.0 -o fine1.nc 

# Coarse -> fine:
$PISM_PATH/pismr -i fine1.nc -regrid_file coarse1.nc -regrid_vars enthalpy $OPTS -o fine2.nc
# Fine -> coarse:
$PISM_PATH/pismr -i coarse1.nc -regrid_file fine2.nc -regrid_vars enthalpy $OPTS -o coarse2.nc

set +e

# Compare:
$PISM_PATH/nccmp.py -t 1e-16 -v enthalpy coarse1.nc coarse2.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

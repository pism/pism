#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 9: 3D regridding from files with different variable orders."
files="foo-09.nc bar-09.nc baz-09.nc"

OPTS="-Mx 61 -My 61 -Mz 21 -Lz 4000 -regrid_file foo-09.nc -regrid_vars topg,litho_temp,thk,bwat,enthalpy -y 0"

rm -f $files

set -e -x

# Create a file to bootstrap from (with a non-trivial bed topography):
$MPIEXEC -n 1 $PISM_PATH/pismr -eisII A -Mx 51 -My 60 -Mz 21 -Mbz 21 -Lbz 1000 -y 0 -o foo-09.nc

# Bootstrap from this file and run for 0 years:
$MPIEXEC -n 2 $PISM_PATH/pismr -i foo-09.nc -bootstrap $OPTS -o bar-09.nc

# Change the variable order in foo-09.nc to z,y,x:
ncpdq -O -a z,y,x foo-09.nc foo-09.nc

# Bootstrap from this file and run for 0 years:
$MPIEXEC -n 2 $PISM_PATH/pismr -i foo-09.nc -bootstrap $OPTS -o baz-09.nc

set +e
set +x

# Compare bar-09.nc and baz-09.nc:
$PISM_PATH/nccmp.py -x -v timestamp bar-09.nc baz-09.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

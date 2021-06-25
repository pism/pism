#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test #32: testing special enthalpy regridding code."
# The list of files to delete when done:
files="foo-32.nc bar-32.nc baz-32.nc"

OPTS="-ys 0 -y 0 -i foo-32.nc -bootstrap -regrid_file foo-32.nc -Lz 4000 -Mx 31 -My 31 -Mz 51"

set -e -x

# Create the file to regrid and bootstrap from:
$PISM_PATH/pismr -eisII A -energy enthalpy -y 500 -o foo-32.nc -o_size big

# Bootstrap from this file:
$PISM_PATH/pismr $OPTS -o bar-32.nc

# Remove enthalpy from foo-32.nc
ncks -x -v enthalpy -O foo-32.nc foo-32.nc

# Bootstrap again
$PISM_PATH/pismr $OPTS -o baz-32.nc

set +e

# Compare:
$PISM_PATH/nccmp.py -r -t 1e-6 -v enthalpy bar-32.nc baz-32.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

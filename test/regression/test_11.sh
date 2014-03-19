#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #11: automatic vertical grid extension."
# The list of files to delete when done.
files="foo-11.nc bar-11.nc baz-11.nc"

rm -f $files

set -x
set -e

OPTS="-My 121 -Mx 61 -eisII A -y 1000 -Mmax 0.925 -z_spacing equal"

echo "run with Lz set too low:"
$MPIEXEC -n 2 $PISM_PATH/pisms -Lz 900 -o foo-11.nc $OPTS

echo "run with Lz set just right:"
$MPIEXEC -n 2 $PISM_PATH/pisms -Mz 33 -Lz 960 -o bar-11.nc $OPTS

echo "regrid from the extended grid onto the one in bar-11.nc:"
$MPIEXEC -n 2 $PISM_PATH/pismr -i bar-11.nc -regrid_file foo-11.nc -regrid_vars enthalpy -y 0 -o baz-11.nc

# compare results
$PISM_PATH/nccmp.py -v enthalpy bar-11.nc baz-11.nc

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0


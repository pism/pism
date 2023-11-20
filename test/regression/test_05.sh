#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 5: bootstrapping from files with different variable orders."
files="foo-05.nc bar-05.nc baz-05.nc"

OPTS="-i foo-05.nc -bootstrap -Mx 41 -My 61 -Mz 21 -Lz 5000 -y 0 -o_size small"

set -e -x

# Create a file to bootstrap from (with a non-trivial bed topography):
$PISM_PATH/pismr -eisII I -Mx 121 -My 61 -Mz 21 -y 0 -o foo-05.nc

# Bootstrap from this file and run for 0 years:
$MPIEXEC -n 2 $PISM_PATH/pismr $OPTS -o bar-05.nc

# Change the variable order in foo-05.nc to z,y,x:
ncpdq -O -a z,y,x foo-05.nc foo-05.nc

# Bootstrap from this file and run for 0 years:
$MPIEXEC -n 2 $PISM_PATH/pismr $OPTS -o baz-05.nc

set +e

# Compare bar-05.nc and baz-05.nc:
$PISM_PATH/nccmp.py -x -v timestamp bar-05.nc baz-05.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

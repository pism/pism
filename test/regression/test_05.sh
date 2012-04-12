#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 5: bootstrapping from files with different variable orders."
files="foo.nc bar.nc baz.nc"

OPTS="-boot_file foo.nc -Mx 41 -My 61 -Mz 21 -Lz 5000 -y 0 -o_size small"

set -e -x

# Create a file to bootstrap from (with a non-trivial bed topography):
$PISM_PATH/pisms -eisII I -Mx 121 -My 61 -Mz 21 -y 0 -o foo.nc 

# Bootstrap from this file and run for 0 years:
$MPIEXEC -n 2 $PISM_PATH/pismr $OPTS -o bar.nc 

# Change the variable order in foo.nc to z,y,x:
ncpdq -O -a z,y,x foo.nc foo.nc

# Bootstrap from this file and run for 0 years:
$MPIEXEC -n 2 $PISM_PATH/pismr $OPTS -o baz.nc 

set +e

# Compare bar.nc and baz.nc:
$PISM_PATH/nccmp.py bar.nc baz.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

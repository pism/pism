#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 4: regridding during bootstrapping: coarse -> fine -> coarse."
files="foo-04.nc bar-04.nc baz-04.nc"

OPTS="-Lz 4000 -y 0 -o_size small "
COARSE="-Mx 11 -My 21 -Mz 31"
FINE=" -Mx 21 -My 31 -Mz 41"

set -e -x

# Create a file to bootstrap from:
$PISM_PATH/pismr -test G $COARSE -y 0 -o foo-04.nc

# Coarse -> fine:
$PISM_PATH/pismr -i foo-04.nc -bootstrap $FINE   $OPTS -o bar-04.nc

# Fine -> coarse:
$PISM_PATH/pismr -i bar-04.nc -bootstrap $COARSE $OPTS -o baz-04.nc

# Compare:
$PISM_PATH/nccmp.py -v topg foo-04.nc baz-04.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

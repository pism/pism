#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 3: no information loss on -y 0 runs (ignoring diagnostics)."
files="foo.nc bar.nc"

OPTS="-o_size small "

set -e -x

# Create a file to start from:
$MPIEXEC -n 2 $PISM_PATH/pisms -no_cold -y 1000 $OPTS -o foo.nc

# Run for 0 years:
$MPIEXEC -n 2 $PISM_PATH/pismr -i foo.nc -y 0 $OPTS -o bar.nc

# Compare, excluding irrelevant diagnostic variables:
$PISM_PATH/nccmp.py foo.nc bar.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

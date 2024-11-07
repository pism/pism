#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 3: no information loss on -y 0 runs (ignoring diagnostics)."
files="foo-03.nc bar-03.nc"

OPTS="-o_size small "

set -e -x

# Create a file to start from:
$MPIEXEC -n 2 $PISM_PATH/pism -eisII A -energy enthalpy -y 1000 $OPTS -o foo-03.nc -Mx 31 -My 41

# Run for 0 years:
$MPIEXEC -n 2 $PISM_PATH/pism -i foo-03.nc -y 0 $OPTS -o bar-03.nc

# Compare, excluding irrelevant diagnostic variables:
$PISM_PATH/pism_nccmp -x -v timestamp foo-03.nc bar-03.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

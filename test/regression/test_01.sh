#!/bin/bash

echo "Test # 1: pismr exact restartability (SIA only)."
PISM_PATH=$1
MPIEXEC=$2

files="foo.nc bar.nc baz.nc joe.nc"
OPTS="-max_dt 1 -o_size small"

set -e -x

# generate an interesting file
$PISM_PATH/pisms -no_cold -Mx 6 -My 6 -Mz 5 -y 5000 -max_dt 500.0 -o baz.nc

# run for ten years, fixed time step
$PISM_PATH/pismr -i baz.nc $OPTS -y 10 -o foo.nc

# chain two five year runs, fixed time step
$PISM_PATH/pismr -i baz.nc $OPTS -y 5 -o joe.nc
$PISM_PATH/pismr -i joe.nc $OPTS -y 5 -o bar.nc

set +e

# Compare output files at year 10:
$PISM_PATH/nccmp.py foo.nc bar.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

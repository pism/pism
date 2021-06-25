#!/bin/bash

echo "Test # 1: pismr exact restartability (SIA only)."
PISM_PATH=$1
MPIEXEC=$2

files="foo-01.nc bar-01.nc baz-01.nc joe-01.nc"
OPTS="-max_dt 1 -o_size small -energy enthalpy"

set -e -x

# generate an interesting file
$PISM_PATH/pismr -eisII A -energy enthalpy -Mx 6 -My 6 -Mz 5 -y 5000 -max_dt 500.0 -o baz-01.nc

# run for ten years, fixed time step
$PISM_PATH/pismr -i baz-01.nc $OPTS -ys 0 -y 10 -o foo-01.nc

# chain two five year runs, fixed time step
$PISM_PATH/pismr -i baz-01.nc $OPTS -ys 0 -y 5 -o joe-01.nc
$PISM_PATH/pismr -i joe-01.nc $OPTS -y 5 -o bar-01.nc

set +e

# Compare output files at year 10:
$PISM_PATH/nccmp.py -v basal_melt_rate_grounded,enthalpy,thk,topg,tillwat foo-01.nc bar-01.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

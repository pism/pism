#!/bin/bash

# Test that sia_forward.py runs.

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3
PYTHON_EXECUTABLE=$5

# List of files to remove when done:
files="in.nc out.nc"

rm -f $files

set -e
set -x

# prepare input files with temp and liqfrac and with temp only
$PISM_PATH/pism -eisII A -Mx 61 -My 61 -y 100 -o in.nc

$PYTHON_EXECUTABLE $PISM_SOURCE_DIR/examples/python/sia_forward.py -i in.nc -o out.nc

if [[ -f in.nc && -f out.nc ]];
then
    rm -f $files; exit 0
else
    exit 1
fi

set +e

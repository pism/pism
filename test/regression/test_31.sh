#!/bin/bash

# Test initialization using "temperature" and "liqfrac" to compute
# enthalpy.

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="temp.nc temp-liqfrac.nc output-1.nc output-2.nc"

rm -f $files

set -e
set -x

# prepare input files with temp and liqfrac and with temp only
$PISM_PATH/pism -eisII A -Mx 3 -My 3 -energy enthalpy -y 10 -o temp-liqfrac.nc -o_size big
ncks -x -v enthalpy -O temp-liqfrac.nc temp-liqfrac.nc
ncks -x -v liqfrac -O temp-liqfrac.nc temp.nc

# try initializing from these:
$PISM_PATH/pism -i temp-liqfrac.nc -o output-1.nc -y 10 -options_left
$PISM_PATH/pism -i temp.nc -o output-2.nc -y 10 -options_left

if [[ -f output-1.nc && -f output-2.nc ]];
then
    rm -f $files; exit 0
else
    exit 1
fi

set +e

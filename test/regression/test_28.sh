#!/bin/bash

# Tests bootstrapping from an incomplete dataset.

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="foo-28.nc bar-28.nc baz-28.nc"

rm -f $files

set -e
set -x

# create a (complete) dataset to bootstrap from:
$MPIEXEC -n 2 $PISM_PATH/pisms -y 100 -o foo-28.nc

OPTS="-i foo-28.nc -bootstrap -Mx 61 -My 61 -Mz 11 -y 10 -Lz 1000"
# bootstrap and run for 100 years:
$MPIEXEC -n 2 $PISM_PATH/pismr $OPTS -o bar-28.nc

# remove topg and tillwat (all contain zeros in the file and will default to zero):
ncks -x -v bmelt,tillwat,dbdt,lon,lat -O foo-28.nc foo-28.nc

# bootstrap and run for 100 years:
$MPIEXEC -n 2 $PISM_PATH/pismr $OPTS -o baz-28.nc

set +e
set +x

# Check results:
$PISM_PATH/nccmp.py bar-28.nc baz-28.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

#!/bin/bash

# Tests bootstrapping from an incomplete dataset.

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="foo.nc bar.nc baz.nc"

rm -f $files

set -e

# create a (complete) dataset to bootstrap from:
$MPIEXEC -n 2 $PISM_PATH/pisms -y 100 -o foo.nc

OPTS="-boot_file foo.nc -Mx 61 -My 61 -Mz 11 -y 10"
# bootstrap and run for 100 years:
$MPIEXEC -n 2 $PISM_PATH/pismr $OPTS -o bar.nc

# remove topg and bwat (all contain zeros in the file and will default to zero):
ncks -x -v bmelt,bwat,dbdt,lon,lat -O foo.nc foo.nc

# bootstrap and run for 100 years:
$MPIEXEC -n 2 $PISM_PATH/pismr $OPTS -o baz.nc

set +e

# Check results:
$PISM_PATH/nccmp.py bar.nc baz.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

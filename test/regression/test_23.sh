#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #23: comparing restart: \"-i\" vs \"-boot_file\" & \"-regrid_file\"."
# The list of files to delete when done.
files="foo.nc bar.nc"

bedrock="-Mbz 11 -Lbz 1000"
opts="-Mx 21 -My 21 -Mz 31 -Lz 4000 $bedrock -o_size small"

rm -f $files

set -e -x

# create foo.nc (at about 6500 years we get some basal melting...)
$MPIEXEC -n 2 $PISM_PATH/pisms -no_cold -y 6500 $opts -o foo.nc

# bootstrap from it, re-gridding all the variables we can
$PISM_PATH/pismr -boot_file foo.nc -regrid_file foo.nc -y 0 -no_temp $opts -o bar.nc

set +e

# compare results (foo.nc and bar.nc contain model_state variables only)
$PISM_PATH/nccmp.py foo.nc bar.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0


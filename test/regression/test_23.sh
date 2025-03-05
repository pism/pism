#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #23: comparing restart: '-i file.nc' vs '-i file.nc -bootstrap' & '-regrid_file'."
# The list of files to delete when done.
files="foo-23.nc bar-23.nc"

bedrock="-Mbz 11 -Lbz 1000"
opts="-Mx 21 -My 21 -Mz 31 -Lz 4000 $bedrock -o_size small"

rm -f $files

set -e -x

# create foo-23.nc (at about 6500 years we get some basal melting...)
$MPIEXEC -n 2 $PISM_PATH/pism -eisII A -energy enthalpy -y 6500 $opts -o foo-23.nc

# bootstrap from it, re-gridding all the variables we can
$PISM_PATH/pism -i foo-23.nc -bootstrap -regrid_file foo-23.nc -y 0 $opts -o bar-23.nc

set +e

# compare results (foo-23.nc and bar-23.nc contain model_state variables only)
$PISM_PATH/pism_nccmp -x -v timestamp foo-23.nc bar-23.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0


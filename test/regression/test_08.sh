#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test # 8: regridding: coarse -> fine -> coarse (vertical direction)."
# The list of files to delete when done.
files="coarse1-08.nc coarse2-08.nc fine1-08.nc fine2-08.nc"

# -grid.registration corner is needed to match default settings of pism -eisII A
OPTS="-y 0 -grid.registration corner"

set -e -x

# Create a file to regrid from:
$PISM_PATH/pism -eisII A -energy enthalpy -Mx 10 -My 10 -Mz 11 -y 6000 -max_dt 300.0 -o coarse1-08.nc
# Create another file with a finer grid:
$PISM_PATH/pism -eisII A -energy enthalpy -Mx 10 -My 10 -Mz 21 -y 6000 -max_dt 300.0 -o fine1-08.nc

# Coarse -> fine:
$PISM_PATH/pism -i fine1-08.nc -regrid_file coarse1-08.nc -regrid_vars enthalpy $OPTS -o fine2-08.nc
# Fine -> coarse:
$PISM_PATH/pism -i coarse1-08.nc -regrid_file fine2-08.nc -regrid_vars enthalpy $OPTS -o coarse2-08.nc

set +e

# Compare:
$PISM_PATH/pism_nccmp -t 1e-16 -v enthalpy coarse1-08.nc coarse2-08.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

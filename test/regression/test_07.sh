#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 7: comparing regridding via '-i file.nc -bootstrap' and -regrid_file."
# The list of files to delete when done:
files="foo-07.nc bar-07.nc baz-07.nc"

OPTS="-y 0"

set -e -x

# Create the file to regrid and bootstrap from:
$PISM_PATH/pism -test G -Lx 4000 -Ly 4000 -Lz 4000 -Mx 41 -My 51 -Mz 31 -y 0 -o foo-07.nc

# Bootstrap from this file:
$PISM_PATH/pism -i foo-07.nc -bootstrap -Lx 2000 -Ly 2000 -Lz 4000 -Mx 31 -My 41 -Mz 51 $OPTS -o bar-07.nc

# Overwrite topg using -regrig_file and save the result to baz-07.nc:
$PISM_PATH/pism -i bar-07.nc -regrid_file foo-07.nc -regrid_vars topg $OPTS -o baz-07.nc

set +e

# Compare:
$PISM_PATH/nccmp.py -v topg bar-07.nc baz-07.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

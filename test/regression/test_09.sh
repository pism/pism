#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 9: 3D regridding from files with different variable orders."
files="foo.nc bar.nc baz.nc"

OPTS="-Mx 61 -My 61 -Mz 21 -Lz 4000 -regrid_file foo.nc -regrid_vars topg,litho_temp,thk,bwat,enthalpy -y 0"

rm -f $files

set -e -x

# Create a file to bootstrap from (with a non-trivial bed topography):
$MPIEXEC -n 1 $PISM_PATH/pisms -eisII I -Mx 51 -My 60 -Mz 21 -Mbz 21 -Lbz 1000 -y 0 -o foo.nc

# Bootstrap from this file and run for 0 years:
$MPIEXEC -n 2 $PISM_PATH/pismr -boot_file foo.nc $OPTS -o bar.nc

# Change the variable order in foo.nc to z,y,x:
ncpdq -O -a z,y,x foo.nc foo.nc

# Bootstrap from this file and run for 0 years:
$MPIEXEC -n 2 $PISM_PATH/pismr -boot_file foo.nc $OPTS -o baz.nc

set +e
set +x

# Compare bar.nc and baz.nc:
$PISM_PATH/nccmp.py bar.nc baz.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

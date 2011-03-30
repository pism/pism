#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #31: testing whether runtime viewers break or not."
# The list of files to delete when done.
files="simp_exper.nc"

rm -f $files

set -e -x

$PISM_PATH/pisms -eisII A -y 1000 -view_sounding temp,litho_temp -view_map velsurf,thk -Mbz 11 -Lbz 1000 -o_size small

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0


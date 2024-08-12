#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #27: testing whether runtime viewers break or not."
# The list of files to delete when done.
files="simp_exper-27.nc"

rm -f $files

set -e -x

$PISM_PATH/pismr -eisII A -y 1000 -view velsurf,thk -Mx 61 -My 61 -Mbz 11 -Lbz 1000 -o_size small -o simp_exper-27.nc

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0


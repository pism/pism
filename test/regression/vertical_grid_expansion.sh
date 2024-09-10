#!/bin/bash

echo "Vertical grid expansion."
PISM_PATH=$1
MPIEXEC=$2

files="input-vertical-grid.nc fixed-vertical-grid.nc ok-vertical-grid.nc too-short-vertical-grid_max_thickness.nc"

short_grid="-Lz 1000 -Mz 11 -z_spacing equal"
tall_grid="-Lz 2000 -Mz 21 -z_spacing equal"

set -x
set -e

rm -f $files

# create an input file
${PISM_PATH}/pismr -eisII A -y 0 -Mx 5 -My 5 -o input-vertical-grid.nc ${short_grid}

# run in two steps: using the short grid and then re-starting with a better one
set +e
${PISM_PATH}/pismr -i input-vertical-grid.nc -bootstrap ${short_grid} \
        -y 3e3 \
        -o too-short-vertical-grid.nc
set -e
${PISM_PATH}/pismr -i too-short-vertical-grid_max_thickness.nc -bootstrap \
        -regrid_file too-short-vertical-grid_max_thickness.nc -allow_extrapolation \
        ${tall_grid} \
        -ye 3e3 \
        -o fixed-vertical-grid.nc

# run using the tall grid right away
${PISM_PATH}/pismr -i input-vertical-grid.nc -bootstrap \
        ${tall_grid} \
        -y 3e3 \
        -o ok-vertical-grid.nc

# compare results (use relative tolerance)
$PISM_PATH/nccmp.py -x -v timestamp -r -t 1e-7 fixed-vertical-grid.nc ok-vertical-grid.nc
if [ $? != 0 ];
then
    exit 1
fi

# compare results again (some variables should match exactly)
$PISM_PATH/nccmp.py -v thk,enthalpy fixed-vertical-grid.nc ok-vertical-grid.nc
if [ $? != 0 ];
then
    exit 1
fi

# clean up
rm -f $files; exit 0

#!/bin/bash

# Test exact restartability of the Lingle-Clark bed deformation model.

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="out0.nc out1.nc out2.nc out3.nc ex1.nc ex2.nc ex.nc ex-restart.nc"

rm -f $files

set -e

mpi="$MPIEXEC -n 3"
pismr="$PISM_PATH/pismr"

# time step length
dt=100
# use a non-square grid
Mx=41
My=61
options="-bed_def lc -extra_times 0,100,200,300 -extra_vars dbdt,topg,thk -stress_balance none -energy none -bed_deformation.lc.update_interval ${dt} -max_dt ${dt}"

grid="-Lz 5000 -Mz 3 -Mx ${Mx} -My ${My}"

# create the input file
${mpi} ${pismr} -eisII A ${grid} -y 1000 -o out0.nc -verbose 1

set -x

# Run with re-starting
#
# Note that this first run stops after 1.5 update intervals, so it runs the bed
# deformation model once at year 100 and saves that in the output file.
${mpi} ${pismr} ${options} -i out0.nc -o out1.nc -extra_file ex1.nc -ys 0 -ye 150 -bootstrap ${grid}
# This run reads the last bed deformation update time from its input file and updates the
# bed at years 200 and 300.
${mpi} ${pismr} ${options} -i out1.nc -o out2.nc -extra_file ex2.nc -ye 300

# run straight
#
# This run updates bed elevation at years 100, 200, and 300.
${mpi} ${pismr} ${options} -bootstrap ${grid} -i out0.nc -o out3.nc -extra_file ex.nc -ys 0 -ye 300

set +x

# combine output files
ncrcat -O ex1.nc ex2.nc ex-restart.nc

set +e

# Compare results:
$PISM_PATH/nccmp.py -v dbdt,topg ex.nc ex-restart.nc
if [ $? != 0 ];
then
    exit 1
fi

# rm -f $files; exit 0

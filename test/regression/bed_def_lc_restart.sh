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
pisms="$PISM_PATH/pisms"
pismr="$PISM_PATH/pismr"

# time step length
dt=100
# use a non-square grid
Mx=41
My=61
options="-bed_def lc -extra_times ${dt} -extra_vars dbdt,topg,thk -stress_balance none -energy none -calendar none -bed_deformation.update_interval 1 -max_dt ${dt}"

grid="-Lz 5000 -Mz 3 -Mx ${Mx} -My ${My}"

# create the input file
${mpi} ${pisms} ${grid} -y 1000 -o out0.nc -verbose 1

# run with re-starting
${mpi} ${pismr} ${options} -i out0.nc -o out1.nc -extra_file ex1.nc -ys 0 -ye 1000 -bootstrap ${grid}
${mpi} ${pismr} ${options} -i out1.nc -o out2.nc -extra_file ex2.nc -ye 2000

# combine output files
ncrcat -O ex1.nc ex2.nc ex-restart.nc

# run straight
${mpi} ${pismr} ${options} -bootstrap ${grid} -i out0.nc -o out3.nc -extra_file ex.nc -ys 0 -ye 2000

set +e

# Compare results:
$PISM_PATH/nccmp.py -v dbdt,topg ex.nc ex-restart.nc
if [ $? != 0 ];
then
    exit 1
fi

# rm -f $files; exit 0

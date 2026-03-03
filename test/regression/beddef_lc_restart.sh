#!/bin/bash

# Test exact restartability of the Lingle-Clark bed deformation model.

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="out0.nc out1.nc out2.nc out3.nc spatial1.nc spatial2.nc spatial.nc spatial-restart.nc"

rm -f $files

set -e

mpi="$MPIEXEC -n 3"
pism="$PISM_PATH/pism"

# Time step length
dt=100
# Use a non-square grid:
Mx=41
My=61

grid="-Lz 5000 -Mz 3 -Mx ${Mx} -My ${My}"

# Create the input file
${mpi} ${pism} -eisII A ${grid} -y 1000 -o out0.nc -verbose 1

set -x

# NB: The -spatial_times frequency (every 50 years) is chosen so that stopping the run at
# year 150 and then re-starting PISM and continuing to year 300 *does not* affect time
# stepping compared to the uninterrupted run from year 0 to year 300.
#
# The stress balance model is turned off, so all ice thickness changes are due to the SMB.
options="-bed_def lc -spatial_times 50 -spatial_vars dbdt,topg,thk -stress_balance none -energy none -bed_deformation.update_interval ${dt}"

# Run with re-starting
#
# Note that this first run stops after 1.5 bed deformation update intervals, so it runs
# the bed deformation model once at year 100 and saves that in the output file.
${mpi} ${pism} ${options} -i out0.nc -o out1.nc -spatial_file spatial1.nc -ys 0 -ye 150 -bootstrap ${grid}
# This run reads the last bed deformation update time from its input file and updates the
# bed at years 200 and 300.
${mpi} ${pism} ${options} -i out1.nc -o out2.nc -spatial_file spatial2.nc -ye 300

# Run without interruptions
#
# This run updates bed elevation at years 100, 200, and 300.
${mpi} ${pism} ${options} -bootstrap ${grid} -i out0.nc -o out3.nc -spatial_file spatial.nc -ys 0 -ye 300

set +x

# Combine output files from the stopped and re-started run:
ncrcat -O spatial1.nc spatial2.nc spatial-restart.nc

set +e

# Compare results:
$PISM_PATH/pism_nccmp -v dbdt,topg spatial.nc spatial-restart.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

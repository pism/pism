#!/bin/bash

# Test the ability to stop and re-start the isochrone tracking model.

set -u

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# create a temporary directory and set up automatic cleanup
temp_dir=$(mktemp -d --tmpdir pism-test-XXXX)
trap 'rm -rf "$temp_dir"' EXIT
cd $temp_dir

set -e

M=50

grid="
        -grid.Lz 5000
        -grid.Mx ${M}
        -grid.My ${M}
        -grid.Mz 21
"

bootstrap="
        -input.bootstrap
        -input.file input.nc
        -isochrones.bootstrapping.n_layers 0
        -time.start 0
"

common_options="
        -isochrones.deposition_times 10
        -stress_balance.model sia
        -stress_balance.sia.flow_law isothermal_glen
        -energy.model none
        -time_stepping.maximum_time_step 1
"

# number of MPI processes:
N=4

end_time=20
stop_time=10

set -x

# generate an input file
${MPIEXEC} -n ${N} ${PISM_PATH}/pism -eisII A \
        ${grid} \
        -output.file input.nc \
        -time.end 1s

# chain together two runs:

# from 0 to A:
${MPIEXEC} -n ${N} ${PISM_PATH}/pism \
           ${bootstrap} ${common_options} \
           -output.file o_part1.nc \
           -time.end ${stop_time}
# and from A to B:
${MPIEXEC} -n ${N} ${PISM_PATH}/pism \
           ${common_options} \
           -input.file o_part1.nc \
           -time.end ${end_time} \
           -output.file o_interrupted.nc

# Now run from 0 to B and compare:
${MPIEXEC} -n ${N} ${PISM_PATH}/pism \
           ${bootstrap} ${common_options} \
           -output.file o_uninterrupted.nc \
           -time.end ${end_time} \


set +e

ignored_vars=wall_clock_time,step_counter,pism_config,model_years_per_processor_hour

# Compare results:
$PISM_PATH/pism_nccmp -x -v ${ignored_vars} o_uninterrupted.nc o_interrupted.nc

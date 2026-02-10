#!/bin/bash

set -e
set -u
set -x

#################################################
# Test asynchronous output using YAC-based code #
#################################################

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3
PYTHON=$5

# Note: async_io_test.py should use "-n X" for some X > 1. This is needed to make sure
# that pism_async_writer uses the code assembling parts of the grid received from
# different ranks (on the PISM side).
${MPIEXEC} -n 3 ${PYTHON} ${PISM_SOURCE_DIR}/test/regression/async_io/async_io_test.py ${arg:-} :\
           -n 1 ${PYTHON} ${PISM_SOURCE_DIR}/util/pism_async_writer -d

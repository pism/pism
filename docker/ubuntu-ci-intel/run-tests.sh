#!/bin/bash

set +x
source /opt/intel/oneapi/setvars.sh
source $HOME/local/pism/bin/activate
set -x

set -e
set -u

# Test using this many jobs
N=${N:-4}

# Location of PISM's build directory
build_dir=${build_dir:-/tmp/pism-build}

# Tell Python where to look for the PISM module
export PYTHONPATH=${build_dir}/site-packages:$PYTHONPATH

# Run tests in parallel
export CTEST_PARALLEL_LEVEL=${N}

ctest --test-dir ${build_dir} --output-on-failure

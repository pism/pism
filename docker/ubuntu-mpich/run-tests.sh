#!/bin/bash

set -x
set -e
set -u

# Test using this many jobs
N=${N:-4}

# Location of PISM's build directory
build_dir=${build_dir:-/tmp/pism-build}

# Run tests in parallel
export CTEST_PARALLEL_LEVEL=${N}

ctest --test-dir ${build_dir} --output-on-failure

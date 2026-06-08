#!/bin/bash

# Activate the venv that build.sh populated (with mpi4py, xarray, ...) so
# tests invoked via "#!/usr/bin/env python3" pick up the venv's interpreter
# instead of falling back to the system Python (which lacks those packages).
set +x
source $HOME/local/pism/bin/activate
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

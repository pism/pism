#!/bin/bash

set -e
set -u
set -x

# The Make target to build
target=${target:-all}

# Compile using this many jobs
N=${N:-4}

# PISM's source directory in the container
source_dir=${source_dir:-$HOME/project}

# The build directory to use
build_dir=${build_dir:-/tmp/pism-build}

# Installation prefix for all prerequisites built from source
lib_prefix=${lib_prefix:-$HOME/local}

# Build Python bindings?
python=${python:-YES}

# Installation prefix
install_dir=${install_dir:-$HOME/local/pism}

git config --global --add safe.directory ${source_dir}

cmake -S ${source_dir} \
      -B ${build_dir} \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_FLAGS="-Werror" \
      -DCMAKE_C_FLAGS="-Werror" \
      -DCMAKE_INSTALL_PREFIX=${install_dir} \
      -DCMAKE_PREFIX_PATH="${lib_prefix}" \
      -DPism_BUILD_ICEBIN=YES \
      -DPism_BUILD_PYTHON_BINDINGS=${python} \
      -DPism_PEDANTIC_WARNINGS=YES \
      -DPism_USE_PARALLEL_NETCDF4=YES \
      -DPism_USE_PNETCDF=YES \
      -DPism_USE_PROJ=YES \
      -DPism_USE_YAC=YES

make --no-print-directory -C ${build_dir} -j ${N} ${target}

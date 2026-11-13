#!/bin/bash

set -x
set -e
set -u

# Uses environment variables CC and CXX to select the compiler

# Define environment variable "minimal" to build the minimal PISM configuration (no
# optional dependencies).
minimal=${minimal:-}

# The Make target to build
target=${target:-all}

# Compile using this many jobs
N=${N:-4}

# CMake installation prefix
CMAKE_PREFIX=${CMAKE_PREFIX:-/usr}

# PISM's source directory in the container
source_dir=${source_dir:-$HOME/project}

# The build directory to use
build_dir=${build_dir:-/tmp/pism-build}

# Clang's static analyzer:
run_analyzer=${run_analyzer:-}
if [ ! -z ${run_analyzer} ];
then
  scan_build=${CC/clang/scan-build}
  analyze="${scan_build} -o ${build_dir}/report --status-bugs -v --exclude ${source_dir}/src/external"
else
  analyze=""
fi

# Installation prefix for all prerequisites built from source
lib_prefix=${lib_prefix:-/usr/local/netcdf-parallel}

# PETSc directory to use
export PETSC_DIR=${PETSC_DIR:-/usr/local}

# Build Python bindings?
python=${python:-YES}

# Installation prefix
install_dir=${install_dir:-$HOME/local/pism}

git config --global --add safe.directory ${source_dir}

EXTRA_CMAKE_FLAGS=""
if [ -z ${minimal} ];
then
  EXTRA_CMAKE_FLAGS="
  -DPism_BUILD_ICEBIN=YES
  -DPism_BUILD_PYTHON_BINDINGS=${python}
  -DPism_USE_PARALLEL_NETCDF4=YES
  -DPism_USE_PNETCDF=YES
  -DPism_USE_PROJ=YES
  -DPism_USE_YAC=YES
"
fi

${analyze} ${CMAKE_PREFIX}/bin/cmake -S ${source_dir} \
               -B ${build_dir} \
               -DCMAKE_BUILD_TYPE=Debug \
               -DCMAKE_CXX_FLAGS="-Werror" \
               -DCMAKE_C_FLAGS="-Werror" \
               -DCMAKE_EXE_LINKER_FLAGS="-fuse-ld=lld" \
               -DCMAKE_INSTALL_PREFIX=${install_dir} \
               -DCMAKE_MODULE_LINKER_FLAGS="-fuse-ld=lld" \
               -DCMAKE_PREFIX_PATH="${lib_prefix}" \
               -DCMAKE_SHARED_LINKER_FLAGS="-fuse-ld=lld" \
               -DPism_PEDANTIC_WARNINGS=YES \
               ${EXTRA_CMAKE_FLAGS}

${analyze} make --no-print-directory -C ${build_dir} -j ${N} ${target}

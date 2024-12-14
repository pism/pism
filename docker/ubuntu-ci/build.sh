#!/bin/bash

set -x
set -e
set -u

# Uses environment variables CC and CXX to select the compiler

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

# Installation prefix for all prerequisites built from source
lib_prefix=${lib_prefix:-$HOME/local}

# PETSc directory to use
PETSC_DIR=${petsc_dir:-${lib_prefix}/petsc}
export PETSC_DIR

# Build Python bindings?
python=${python:-YES}

# Installation prefix
install_dir=${install_dir:-$HOME/local/pism}

hdf5_dir=${lib_prefix}/hdf5
netcdf_dir=${lib_prefix}/netcdf
pnetcdf_dir=${lib_prefix}/pnetcdf
yac_dir=${lib_prefix}/yac

git config --global --add safe.directory ${source_dir}

mkdir -p ${build_dir}

${CMAKE_PREFIX}/bin/cmake -S ${source_dir} \
               -B ${build_dir} \
               -DCMAKE_BUILD_TYPE=Debug \
               -DCMAKE_CXX_FLAGS="-Werror" \
               -DCMAKE_C_FLAGS="-Werror" \
               -DCMAKE_EXE_LINKER_FLAGS="-fuse-ld=lld" \
               -DCMAKE_MODULE_LINKER_FLAGS="-fuse-ld=lld" \
               -DCMAKE_SHARED_LINKER_FLAGS="-fuse-ld=lld" \
               -DCMAKE_PREFIX_PATH="${hdf5_dir};${netcdf_dir};${pnetcdf_dir};${yac_dir}" \
               -DPism_BUILD_PYTHON_BINDINGS=${python} \
               -DPism_BUILD_EXTRA_EXECS=YES \
               -DPism_PEDANTIC_WARNINGS=YES \
               -DPism_BUILD_ICEBIN=YES \
               -DPism_USE_PROJ=YES \
               -DPism_USE_PARALLEL_NETCDF4=YES \
               -DPism_USE_PNETCDF=YES \
               -DPism_USE_YAC_INTERPOLATION=YES \
               -DCMAKE_INSTALL_PREFIX=${install_dir}

make --no-print-directory -C ${build_dir} -j ${N} ${target}

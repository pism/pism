#!/bin/bash

set -x
set -e
set -u

# Compile using this many jobs
N=${N:-4}

# PISM's source directory in the container
source_dir=${source_dir:-/var/tmp/pism}

# The build directory to use
build_dir=${build_dir:-/tmp/pism-build}

# Installation prefix for all prerequisites built from source
lib_prefix=/opt
# Installation prefix
install_dir=/opt/pism

hdf5_dir=${lib_prefix}/hdf5
netcdf_dir=${lib_prefix}/netcdf
pnetcdf_dir=${lib_prefix}/pnetcdf
yac_dir=${lib_prefix}/yac

git config --global --add safe.directory ${source_dir}

cmake -S ${source_dir} \
      -B ${build_dir} \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_COMPILER=icpx \
      -DCMAKE_CXX_FLAGS="-Werror -Rno-debug-disables-optimization" \
      -DCMAKE_C_COMPILER=icx \
      -DCMAKE_C_FLAGS="-Werror -Rno-debug-disables-optimization" \
      -DCMAKE_INSTALL_PREFIX=${install_dir} \
      -DCMAKE_PREFIX_PATH="${hdf5_dir};${netcdf_dir};${pnetcdf_dir};${yac_dir}" \
      -DPism_BUILD_ICEBIN=YES \
      -DPism_BUILD_PYTHON_BINDINGS=YES \
      -DPism_PEDANTIC_WARNINGS=YES \
      -DPism_USE_PARALLEL_NETCDF4=YES \
      -DPism_USE_PNETCDF=YES \
      -DPism_USE_PROJ=YES \
      -DPism_USE_YAC_INTERPOLATION=YES

make --no-print-directory -C ${build_dir} -j ${N} all

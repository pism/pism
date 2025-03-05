#!/bin/bash

set -x
set -e
set -u

# PISM's source directory in the container
source_dir=${source_dir:-$HOME/project}

# The build directory to use
build_dir=${build_dir:-/tmp/pism-build}

# Installation prefix for all prerequisites built from source
lib_prefix=${lib_prefix:-$HOME/local}

# PETSc directory to use
export PETSC_DIR=${petsc_dir:-${lib_prefix}/petsc}

# Installation prefix
install_dir=${install_dir:-$HOME/local/pism}

hdf5_dir=${lib_prefix}/hdf5
netcdf_dir=${lib_prefix}/netcdf
pnetcdf_dir=${lib_prefix}/pnetcdf
yac_dir=${lib_prefix}/yac

git config --global --add safe.directory ${source_dir}

scan_build=${scan_build:-scan-build-18}
analyze="${scan_build} -o ${build_dir}/report --status-bugs -v --exclude ${source_dir}/src/external/cubature"

${analyze} cmake -S ${source_dir} \
               -B ${build_dir} \
               -DCMAKE_BUILD_TYPE=Debug \
               -DCMAKE_CXX_FLAGS="-Werror" \
               -DCMAKE_C_FLAGS="-Werror" \
               -DCMAKE_INSTALL_PREFIX=${install_dir} \
               -DCMAKE_PREFIX_PATH="${hdf5_dir};${netcdf_dir};${pnetcdf_dir};${yac_dir}" \
               -DPism_BUILD_EXTRA_EXECS=YES \
               -DPism_BUILD_ICEBIN=YES \
               -DPism_BUILD_PYTHON_BINDINGS=NO \
               -DPism_PEDANTIC_WARNINGS=YES \
               -DPism_USE_PARALLEL_NETCDF4=YES \
               -DPism_USE_PNETCDF=YES \
               -DPism_USE_PROJ=YES \
               -DPism_USE_YAC_INTERPOLATION=YES

${analyze} make --no-print-directory -C ${build_dir} -j ${N:-4} all

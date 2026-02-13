#!/bin/bash

. /opt/intel/oneapi/setvars.sh

set -e
set -u
set -x

N=${N:-8}
opt_flags="-O3 -axCORE-AVX512,CORE-AVX2,SSE4.2 -fp-model=precise"

# Prerequisites:
export PETSC_DIR=${PETSC_DIR:-/opt/petsc}
yaxt_prefix=${yaxt_prefix:-/opt/yaxt}
yac_prefix=${yac_prefix:-/opt/yac}
hdf5_prefix=${hdf5_prefix:-/opt/hdf5}
netcdf_prefix=${netcdf_prefix:-/opt/netcdf}

# Installation prefix and build location:
prefix=${prefix:-/opt/pism}
build_dir=${build_dir:-/var/tmp/build/pism}

mkdir -p ${build_dir}

export CC=mpiicx
export CXX=mpiicpx

cmake \
    -B ${build_dir} \
    -S ${src_dir} \
    -DCMAKE_PREFIX_PATH="${yaxt_prefix};${yac_prefix};${hdf5_prefix};${netcdf_prefix}" \
    -DCMAKE_BUILD_TYPE="Release" \
    -DCMAKE_CXX_FLAGS="${opt_flags} -diag-disable=cpu-dispatch,10006,2102" \
    -DCMAKE_C_FLAGS="${opt_flags} -diag-disable=cpu-dispatch,10006" \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DPism_BUILD_PYTHON_BINDINGS=YES \
    -DPism_USE_YAC=YES \
    -DPism_USE_PARALLEL_NETCDF4=YES \
    -DPism_USE_PROJ=YES  

make -j $N -C ${build_dir} install

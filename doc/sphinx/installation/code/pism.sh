#!/bin/bash

set -x

prefix=$HOME/local/
build_dir=$HOME/local/build/pism

mkdir -p ${build_dir}
pushd ${build_dir}

git clone https://github.com/pism/pism.git .

mkdir -p build
pushd build

export PETSC_DIR=/usr/lib/petsc

CC="mpicc" CXX="mpicxx" cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=${prefix}/pism \
    -DCMAKE_C_FLAGS="-march=native -mtune=native" \
    -DCMAKE_CXX_FLAGS="-march=native -mtune=native" \
    -DPism_BUILD_EXTRA_EXECS=YES \
    -DPython_ADDITIONAL_VERSIONS=3 \
    -DPism_BUILD_PYTHON_BINDINGS=YES \
    -DPism_USE_PNETCDF=YES \
    -DPism_USE_PROJ=YES \
    ..

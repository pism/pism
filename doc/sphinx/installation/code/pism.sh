#!/bin/bash

set -x
set -e

prefix=$HOME/local/
build_dir=$HOME/local/build/pism

mkdir -p ${build_dir}
pushd ${build_dir}

# FIXME: remove "-b dev"
git clone --depth=1 -b dev https://github.com/pism/pism.git .

mkdir -p build
pushd build

export PETSC_DIR=/usr/lib/petsc

CC=mpicc CXX=mpicxx cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=${prefix}/pism \
    -DCMAKE_C_FLAGS="-march=native -mtune=native" \
    -DCMAKE_CXX_FLAGS="-march=native -mtune=native" \
    -DPism_BUILD_PYTHON_BINDINGS=YES \
    -DPism_USE_PNETCDF=YES \
    -DPism_USE_PROJ=YES \
    ..

make -j8 all
# PYTHONPATH=$PWD/site-packages:$PYTHONPATH ctest --output-on-failure
make install
popd
popd

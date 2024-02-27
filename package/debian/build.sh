#!/bin/bash

# This script is used to create the Debian/Ubuntu package.

set -u
set -e
set -x

build_dir=${build_dir:-/var/tmp/pism-build-deb}
source_dir=${source_dir:?Please set source_dir}

rm -rf ${build_dir}
mkdir -p ${build_dir}

export PETSC_DIR=/usr/lib/petsc

CC=mpicc CXX=mpicxx cmake -S ${source_dir} -B ${build_dir} \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/usr \
    -DPism_BUILD_DOCS=YES \
    -DPism_BUILD_EXTRA_EXECS=NO \
    -DPism_BUILD_PYTHON_BINDINGS=NO \
    -DPism_USE_PROJ=YES

make -j4 -C ${build_dir} package

mv -v ${build_dir}/*.deb .

rm -rf ${build_dir}

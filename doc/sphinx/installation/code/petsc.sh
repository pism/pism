#!/bin/bash

set -e
set -u
set -x

# Install the latest PETSc in ~/local/petsc using ~/local/build/petsc as the build
# directory.

prefix=$HOME/local/petsc
build_dir=~/local/build/petsc

rm -rf ${build_dir}
mkdir -p ${build_dir}
pushd ${build_dir}

git clone -b release --depth=1 https://gitlab.com/petsc/petsc.git .

PETSC_DIR=$PWD
PETSC_ARCH="linux-opt"

./configure \
  --prefix=${prefix} \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpifort \
  --with-shared-libraries \
  --with-debugging=0 \
  --with-petsc4py \
  --download-f2cblaslapack

export PYTHONPATH=${prefix}/lib
make all
make install
make PETSC_DIR=${prefix} PETSC_ARCH="" check

popd

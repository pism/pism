#!/bin/bash

set -e
set -u
set -x

# Install the latest PETSc with 64 bit indices in ~/local/petsc using ~/local/build/petsc
# as the build directory.

build_dir=~/local/build/petsc

rm -rf ${build_dir}
mkdir -p ${build_dir}
pushd ${build_dir}

git clone -b release --depth=1 https://gitlab.com/petsc/petsc.git .

petsc_prefix=$HOME/local/petsc-64bit
PETSC_DIR=$PWD
PETSC_ARCH="linux-opt"

./configure \
  COPTFLAGS="-g -O3" \
  --prefix=${petsc_prefix} \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpifort \
  --with-shared-libraries \
  --with-debugging=0 \
  --with-petsc4py \
  --with-x=0 \
  --with-64-bit-indices \
  --download-f2cblaslapack

export PYTHONPATH=${petsc_prefix}/lib
make all
make install
make PETSC_DIR=${petsc_prefix} PETSC_ARCH="" check
popd

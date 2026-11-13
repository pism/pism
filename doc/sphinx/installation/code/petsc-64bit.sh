#!/bin/bash

set -e
set -u
set -x

# Install the latest PETSc with 64 bit indices in /usr/local/petsc-64bit using
# /tmp/build/petsc-64bit as the build directory.

build_dir=${build_dir:-/tmp/build/petsc-64bit}

rm -rf ${build_dir}
mkdir -p ${build_dir}
pushd ${build_dir}

git clone -b release --depth=1 https://gitlab.com/petsc/petsc.git .

prefix=/usr/local/petsc-64bit
PETSC_DIR=$PWD
PETSC_ARCH="linux-opt"

./configure \
  COPTFLAGS="-g -O3" \
  --prefix=${prefix} \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpifort \
  --with-shared-libraries \
  --with-debugging=0 \
  --with-petsc4py \
  --with-x=0 \
  --with-64-bit-indices \
  --download-f2cblaslapack

export PYTHONPATH=${prefix}/lib
make all
make install
make PETSC_DIR=${prefix} PETSC_ARCH="" check
popd

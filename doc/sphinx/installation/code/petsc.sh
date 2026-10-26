#!/bin/bash

set -e
set -u
set -x

# Install the latest PETSc in /usr/local/ using /tmp/build/petsc as the build directory.

build_dir=${build_dir:-/tmp/build/petsc}

rm -rf ${build_dir}
mkdir -p ${build_dir}
pushd ${build_dir}

git clone -b release --depth=1 https://gitlab.com/petsc/petsc.git .

# manual-begin
prefix=${prefix:-/usr/local/}
petsc4py=${petsc4py:-1}
PETSC_DIR=$PWD
PETSC_ARCH="linux-opt"

./configure \
  COPTFLAGS="-g -O3" \
  --prefix=${prefix} \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=0 \
  --with-shared-libraries \
  --with-debugging=0 \
  --with-petsc4py=${petsc4py} \
  --with-x=0 \
  --download-f2cblaslapack

make all
make install
make PETSC_DIR=${prefix} PETSC_ARCH="" check
# manual-end
popd

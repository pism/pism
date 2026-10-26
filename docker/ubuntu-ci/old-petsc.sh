#!/bin/bash

set -e
set -u
set -x

# Install an old PETSc version in /usr/local/petsc-${version} using
# /tmp/build/petsc-${version} as the build directory.

# The default below corresponds to the oldest version supported by PISM.
version=${version:-3.22.0}
build_dir=/tmp/build/petsc-${version}

rm -rf ${build_dir}
mkdir -p ${build_dir}
cd ${build_dir}

git clone -b v${version} --depth=1 https://gitlab.com/petsc/petsc.git .

prefix=/usr/local/petsc-${version}
PETSC_DIR=$PWD
PETSC_ARCH="linux-opt"

python3 ./configure \
  CFLAGS="-fPIC" \
  --prefix=${prefix} \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=0 \
  --with-shared-libraries \
  --with-debugging=0 \
  --with-x=0 || (cat configure.log && exit 1)

make all
make install
make PETSC_DIR=${prefix} PETSC_ARCH="" check

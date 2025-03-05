#!/bin/bash

set -e
set -u
set -x

# Install an old PETSc version in ~/local/petsc-${version} using
# ~/local/build/petsc-${version} as the build directory.

# The default below corresponds to the oldest version supported by PISM.
version=${version:-3.11.0}
build_dir=~/local/build/petsc-${version}

rm -rf ${build_dir}
mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${version}.tar.gz

tar xzf petsc-lite-${version}.tar.gz

pushd petsc-${version}

petsc_prefix=$HOME/local/petsc-${version}
PETSC_DIR=$PWD
PETSC_ARCH="linux-opt"

python3 ./configure \
  CFLAGS="-fPIC" \
  --prefix=${petsc_prefix} \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=0 \
  --with-shared-libraries \
  --with-debugging=0 \
  --with-x=0 \
  --download-f2cblaslapack

make all
make install
make PETSC_DIR=${petsc_prefix} PETSC_ARCH="" check

popd
popd

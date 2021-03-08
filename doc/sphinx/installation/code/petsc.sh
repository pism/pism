#!/bin/bash

set -e
set -u
set -x

# Install the latest PETSc release in ~/local/petsc
# using ~/local/build/petsc as the build directory.

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
  COPTFLAGS='-O3 -march=native -mtune=native' \
  CXXOPTFLAGS='-O3 -march=native -mtune=native' \
  FOPTFLAGS='-O3 -march=native -mtune=native' \
  --download-f2cblaslapack=1 \
  --download-petsc4py \
  --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch \
  --download-hypre

make all
make install
make PETSC_DIR=${prefix} PETSC_ARCH="" check

popd

#!/bin/bash

set -e
set -u
set -x

# Install the latest PETSc with MUMPS in ~/local/petsc-mumps using ~/local/build/petsc-mumps
# as the build directory.

prefix=$HOME/local/petsc-mumps
build_dir=~/local/build/petsc-mumps

rm -rf ${build_dir}
mkdir -p ${build_dir}
pushd ${build_dir}

git clone -b release --depth=1 https://gitlab.com/petsc/petsc.git .

PETSC_DIR=$PWD
PETSC_ARCH="linux-opt"

./configure \
  COPTFLAGS='-O3 -march=native -mtune=native' \
  CXXOPTFLAGS='-O3 -march=native -mtune=native' \
  FOPTFLAGS='-O3 -march=native -mtune=native' \
  --prefix=${prefix} \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpifort \
  --with-shared-libraries \
  --with-debugging=0 \
  --with-petsc4py \
  --download-scalapack \
  --download-mumps \
  --download-parmetis \
  --download-metis \
  --download-ptscotch \
  --download-f2cblaslapack --download-blis \
  --with-valgrind=0 \
  --with-x=0

export PYTHONPATH=${prefix}/lib
make all
make install
make PETSC_DIR=${prefix} PETSC_ARCH="" check

popd

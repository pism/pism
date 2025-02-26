#!/bin/bash

set -e
set -u
set -x

# Install PETSc in /opt/petsc using /var/tmp/build/petsc as the build
# directory.

MPICC=${MPICC:-mpiicx}
MPICXX=${MPICXX:-mpiicx}
COPTFLAGS=${COPTFLAGS:--O3 -march=native -mtune=native -fp-model=precise}

build_dir=${build_dir:-/var/tmp/build/petsc}
prefix=${prefix:-/opt/petsc}

mkdir -p ${build_dir}
cd ${build_dir}
version=${version:-3.22.2}

wget -nc \
     https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${version}.tar.gz
rm -rf petsc-${version}
tar xzf petsc-${version}.tar.gz

cd petsc-${version}

PETSC_DIR=$PWD
PETSC_ARCH="linux-opt"

python3 ./configure \
        COPTFLAGS="${COPTFLAGS}" \
        --prefix=${prefix} \
        --with-cc="${MPICC}" \
        --with-cxx="${MPICXX}" \
        --with-debugging=0 \
        --with-fc=0 \
        --with-petsc4py \
        --with-shared-libraries \
        --with-x=0 \
        --download-f2cblaslapack

make all
make install
make PETSC_DIR=${prefix} PETSC_ARCH="" check

#!/bin/bash

set -e
set -u
set -x

# Install YAXT 0.9.0 in ~/local/cdipio/,
# using ~/local/build/cdipio/yaxt as the build directory.

prefix=$HOME/local/cdipio/
version=0.9.0
build_dir=~/local/build/cdipio/yaxt
url=https://www.dkrz.de/redmine/attachments/download/498/yaxt-0.9.0.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar xzf yaxt-${version}.tar.gz

pushd yaxt-${version}

ignore_mpi_defect=""
#ignore_mpi_defect="--without-regard-for-quality"

export CC=mpicc
export CFLAGS="-std=gnu99 -O3"
./configure --prefix=${prefix} \
            --enable-shared \
            --enable-static \
            ${ignore_mpi_defect}
make
make install

popd
popd

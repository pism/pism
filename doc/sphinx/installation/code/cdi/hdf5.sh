#!/bin/bash

set -e
set -u
set -x

# Install HDF5 1.10.5 with parallel I/O in ~/local/hdf5,
# using ~/local/build/hdf5 as the build directory.

version=1.10.5
prefix=$HOME/local/cdipio/
build_dir=~/local/build/cdipio/hdf5
hdf5_site=https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10
url=${hdf5_site}/hdf5-${version}/src/hdf5-${version}.tar.gz

rm -rf ${build_dir}
mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar xzf hdf5-${version}.tar.gz

pushd hdf5-${version}

CC=mpicc CFLAGS=-w ./configure \
  --disable-static \
  --enable-parallel \
  --prefix=${prefix} 2>&1 | tee hdf5_configure.log

make all 2>&1 | tee hdf5_compile.log
make install 2>&1 | tee hdf5_install.log

popd
popd

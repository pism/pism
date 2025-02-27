#!/bin/bash

set -e
set -u
set -x

# Install PnetCDF 1.12.1 in ~/local/pnetcdf,
# using ~/local/build/pnetcdf as a build directory.

version=1.12.1
prefix=${prefix:-$HOME/local/pnetcdf}
build_dir=${build_dir:-$HOME/local/build/pnetcdf/}
url=https://parallel-netcdf.github.io/Release/pnetcdf-${version}.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar xzf pnetcdf-${version}.tar.gz

pushd pnetcdf-${version}

export CC=${CC:-mpicc}

./configure \
      --prefix=${prefix} \
      --enable-shared \
      --disable-static \
      --disable-cxx \
      --disable-fortran

make all
make install

popd
popd

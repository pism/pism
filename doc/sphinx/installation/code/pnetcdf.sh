#!/bin/bash

set -e
set -u
set -x

# Install PnetCDF 1.14.1 in /usr/local/,
# using /tmp/build/pnetcdf as a build directory.

version=1.14.1
prefix=${prefix:-/usr/local/}
build_dir=${build_dir:-/tmp/build/pnetcdf/}
url=https://parallel-netcdf.github.io/Release/pnetcdf-${version}.tar.gz

rm -rf ${build_dir}
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

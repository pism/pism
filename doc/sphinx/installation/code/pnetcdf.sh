#!/bin/bash

set -e
set -u
set -x

# Install PnetCDF 1.12.0 in ~/local/pnetcdf,
# using ~/local/build/pnetcdf as a build directory.

version=1.12.0
prefix=$HOME/local/pnetcdf
build_dir=~/local/build/pnetcdf/
url=https://parallel-netcdf.github.io/Release/pnetcdf-${version}.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar xzf pnetcdf-${version}.tar.gz

pushd pnetcdf-${version}

CFLAGS="-fPIC" ./configure \
      --prefix=${prefix} \
      --enable-shared \
      --disable-cxx \
      --disable-fortran 2>&1 | tee pnetcdf_configure.log

make all 2>&1 | tee pnetcdf_compile.log
make install 2>&1 | tee pnetcdf_install.log

popd
popd

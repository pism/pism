#!/bin/bash

set -u
set -e
set -x

# Install NetCDF Fortran using NetCDF in ~/local/netcdf.

netcdf=~/local/netcdf

version=4.4.3
prefix=$HOME/local/netcdff
build_dir=~/local/build/netcdff
url=https://github.com/Unidata/netcdf-fortran/archive/v4.4.3.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar zxf v${version}.tar.gz

pushd netcdf-fortran-${version}

LD_LIBRARY_PATH=${netcdf}/lib:${LD_LIBRARY_PATH}
CPPFLAGS=-I${netcdf}/include
LDFLAGS=-L${netcdf}/lib
./configure --prefix=${prefix} --disable-fortran-type-check
make check
make install

popd
popd
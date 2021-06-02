#!/bin/bash

set -u
set -e
set -x

# Install NetCDF with PNetCDF support using parallel HDF5 
# in ~/local/hdf5 and and PNetCDF in ~/local/pnetcdf and 
# ~/local/build/netcdf as a build directory.

hdf5=~/local/hdf5
pnetcdf=~/local/pnetcdf

version=4.4.0
prefix=$HOME/local/netcdf
build_dir=~/local/build/netcdf
url=https://github.com/Unidata/netcdf-c/archive/v4.4.0.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar zxf v${version}.tar.gz

pushd netcdf-c-${version}

CPPFLAGS="-I${hdf5}/include -I${pnetcdf}/include"
LDFLAGS="-L${hdf5}/lib -L${pnetcdf}/lib"

./configure --enable-pnetcdf --enable-parallel-tests --prefix=${prefix}
make check
make Install

popd
popd
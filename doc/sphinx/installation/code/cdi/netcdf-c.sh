#!/bin/bash

set -u
set -e
set -x

# Install NetCDF with PNetCDF support using parallel HDF5 
# in ~/local/cdipio/hdf5 and and PNetCDF in ~/local/pnetcdf and
# ~/local/build/cdipio/netcdf-c as a build directory.

hdf5=~/local/cdipio/
pnetcdf=~/local/pnetcdf

version=4.4.0
prefix=$HOME/local/cdipio/
build_dir=~/local/build/cdipio/netcdf-c
url=https://github.com/Unidata/netcdf-c/archive/v4.4.0.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar zxf v${version}.tar.gz

pushd netcdf-c-${version}

export CC=mpicc
export CPPFLAGS="-I${hdf5}/include -I${pnetcdf}/include"
export LDFLAGS="-L${hdf5}/lib -L${pnetcdf}/lib"

./configure --enable-pnetcdf \
            --enable-parallel-tests \
            --prefix=${prefix} | tee netcdf-configure.log

make check | tee netcdf-make-check.log
make install | tee netcdf-make-install.log

popd
popd

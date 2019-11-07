#!/bin/bash

set -u
set -e
set -x

# Install parallel NetCDF using parallel HDF5 in ~/local/hdf5 and
# ~/local/build/netcdf as a build directory.

hdf5=~/local/hdf5

# 4.7.1 does not work on macOS, 4.7.2 has a bug #1502 that breaks parallel I/O.
# Use 4.7.3 once it is released.
version=4.7.0
prefix=$HOME/local/netcdf
build_dir=~/local/build/netcdf
url=ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-c-${version}.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar zxf netcdf-c-${version}.tar.gz

pushd netcdf-c-${version}

CC=${hdf5}/bin/h5pcc ./configure \
  --enable-netcdf4 \
  --disable-dap \
  --prefix=${prefix} 2>&1 | tee netcdf_configure.log

make all 2>&1 | tee netcdf_compile.log
make install 2>&1 | tee netcdf_install.log

popd
popd

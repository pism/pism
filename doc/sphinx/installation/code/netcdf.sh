#!/bin/bash

set -u
set -e
set -x

# Install parallel NetCDF using parallel HDF5 and PnetCDF in ${lib_dir} and
# ~/local/build/netcdf as a build directory.

# directory containing HDF5 and PnetCDF
lib_dir=${lib_dir:-$HOME/local}

version=4.9.3
prefix=${prefix:-$HOME/local/netcdf}
build_dir=${build_dir:-$HOME/local/build/netcdf}
url=https://github.com/Unidata/netcdf-c/archive/refs/tags/v${version}.tar.gz

mkdir -p ${build_dir}
cd ${build_dir}

wget -nc ${url}
tar zxf v${version}.tar.gz

cd netcdf-c-${version}

export CPPFLAGS="-I${lib_dir}/include"
export LDFLAGS="-L${lib_dir}/lib"
export CC=${CC:-mpicc}

./configure \
  --enable-netcdf4 \
  --enable-parallel4 \
  --enable-pnetcdf \
  --disable-dap \
  --disable-libxml2 \
  --disable-byterange \
  --prefix=${prefix} 2>&1 | tee netcdf_configure.log

make all 2>&1 | tee netcdf_compile.log
make install 2>&1 | tee netcdf_install.log

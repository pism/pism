#!/bin/bash

. /opt/intel/oneapi/setvars.sh

set -u
set -e
set -x

# Install parallel NetCDF using parallel HDF5 in $hdf5_prefix and
# /var/tmp/build/netcdf as a build directory.

netcdf_version=${netcdf_version:-4.9.3}
prefix=${prefix:-/opt/netcdf}
build_dir=${build_dir:-/var/tmp/build/netcdf}
url=https://github.com/Unidata/netcdf-c/archive/refs/tags/v${netcdf_version}.tar.gz

mkdir -p ${build_dir}
cd ${build_dir}

wget -nc ${url}
tar zxf v${netcdf_version}.tar.gz

cd netcdf-c-${netcdf_version}

export CPPFLAGS="-I${hdf5_prefix}/include"
export LDFLAGS="-L${hdf5_prefix}/lib"
export CC=${CC:-mpiicx}

./configure \
  --enable-netcdf4 \
  --enable-parallel4 \
  --disable-dap \
  --disable-libxml2 \
  --disable-byterange \
  --prefix=${prefix} 2>&1 | tee netcdf_configure.log

make all 2>&1 | tee netcdf_compile.log
make install 2>&1 | tee netcdf_install.log

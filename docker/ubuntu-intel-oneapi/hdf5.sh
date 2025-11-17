#!/bin/bash

. /opt/intel/oneapi/setvars.sh

set -e
set -u
set -x

# Install HDF5 with parallel I/O in ${prefix} using ${build_dir} as
# the build directory.

# FIXME: HDF5 directory hierarchy is weird. This scripts will break
# when hdf5_version is not 1.12.x.
hdf5_version=${hdf5_version:-1.12.1}
prefix=${prefix:-/opt/hdf5/}
build_dir=${build_dir:-/var/tmp/build/hdf5}
hdf5_site=https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12
url=${hdf5_site}/hdf5-${hdf5_version}/src/hdf5-${hdf5_version}.tar.gz

mkdir -p ${build_dir}
cd ${build_dir}

wget -nc ${url}
tar xzf hdf5-${hdf5_version}.tar.gz

cd hdf5-${hdf5_version}

./configure \
  CC=${CC:-mpiicx} \
  CXX=${CC:-mpiicpx} \
  CFLAGS=-w \
  --disable-static \
  --enable-parallel \
  --prefix=${prefix} 2>&1 | tee hdf5_configure.log

make all 2>&1 | tee hdf5_compile.log
make install 2>&1 | tee hdf5_install.log

#!/bin/bash

# This script is used to build HDF5 and NetCDF with parallel I/O support in the
# Debian-based Docker image used to run PISM's regression tests.

set -x
set -e

prefix=${HOME}/local/
N=4

#!/bin/bash

build_hdf5() {
    # download and build HDF5
    build_dir=${prefix}/build/hdf5/
    mkdir -p ${build_dir}
    cd ${build_dir}
    version=1.8.21

    wget -nv -nc https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-${version}/src/hdf5-${version}.tar.gz

    tar xzf hdf5-${version}.tar.gz

    cd hdf5-${version}
    CC=mpicc CFLAGS=-w ./configure --enable-parallel --prefix=${prefix}/hdf5 2>&1 | tee hdf5_configure.log

    make all -j $N 2>&1 | tee hdf5_compile.log
    make install 2>&1 | tee hdf5_install.log

    cd ${prefix}
    rm -rf ${prefix}/build/hdf5
}

build_netcdf() {
    # download and build netcdf
    build_dir=${prefix}/build/netcdf/
    mkdir -p ${build_dir}
    cd ${build_dir}
    version=4.6.1

    wget -nv -nc ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-${version}.tar.gz
    tar xzf netcdf-${version}.tar.gz

    cd netcdf-${version}
    export CC=mpicc
    export CPPFLAGS="-I${prefix}/hdf5/include" LDFLAGS="-L${prefix}/hdf5/lib"

    ./configure \
        --enable-netcdf4 \
        --disable-dap \
        --prefix=${prefix}/netcdf 2>&1 | tee netcdf_configure.log

    make all -j $N 2>&1 | tee netcdf_compile.log
    make install 2>&1 | tee netcdf_install.log

    cd ${prefix}
    rm -rf ${prefix}/build/netcdf
}

build_hdf5
build_netcdf

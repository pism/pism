#!/bin/bash

set -e
set -u
set -x

# Install PPM in ~/local/cdipio using
# parallel HDF5, NetCDF-C, and NetCDF-Fortran in
# ~/local/cdipio/
hdf5=~/local/cdipio/
netcdf=~/local/cdipio/
netcdff=~/local/cdipio/

prefix=$HOME/local/cdipio/
build_dir=~/local/build/cdipio/ppm
url=https://www.dkrz.de/redmine/attachments/download/500/ppm-1.0.6.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar zxf ppm-1.0.6.tar.gz

pushd ppm-1.0.6

export CC=mpicc
export FC=mpifort
./configure --prefix=${prefix} \
            --disable-parmetis \
            --disable-metis \
            --enable-MPI \
            --enable-hdf5 \
            --with-hdf5-root=${hdf5} \
            --enable-netcdf \
            --with-netcdf-root=${netcdf} \
            --with-netcdf-fc-lib=${netcdff}/lib \
            --with-netcdf-fc-mod=${netcdff}/include

make
make install

popd
popd

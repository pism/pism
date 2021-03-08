#!/bin/bash

set -e
set -u
set -x

# Install PPM 1.0.6 using parallel HDF5 in ~/local/hdf5, NetCDF in 
# ~/local/netcdf and NetCDF Fortran in ~/local/netcdff
hdf5=~/local/hdf5
netcdf=~/local/netcdf
netcdff=~/local/netcdff

version=1.0.6
prefix=$HOME/local/PPM
build_dir=~/local/build/PPM
url=https://www.dkrz.de/redmine/attachments/download/500/ppm-1.0.6.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar zxf ppm-${version}.tar.gz

push ppm-${version}

NETCDF_FC_LIB=$(${netcdff}/bin/nf-config --flibs)
NETCDF_FC_MOD=$(${netcdff}/bin/nf-config --fflags)
./configure --prefix=${prefix} --enable-MPI --enable-netcdf \
--enable-hdf5 --disable-parmetis --disable-metis --with-hdf5-root=${hdf5} \
--with-netcdf-include=${netcdf}/include --with-netcdf-lib=${netcdf}/lib

make
make install

popd
popd
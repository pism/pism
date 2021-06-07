#!/bin/bash

set -u
set -e
set -x

# Install NetCDF-Fortran in ~/local/cdipio
# using NetCDF-C in ~/local/cdipio.

netcdf=~/local/cdipio/

version=4.4.3
prefix=$HOME/local/cdipio/
build_dir=~/local/build/cdipio/netcdf-f
url=https://github.com/Unidata/netcdf-fortran/archive/v4.4.3.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar zxf v${version}.tar.gz

pushd netcdf-fortran-${version}

export CPPFLAGS=-I${netcdf}/include
export LDFLAGS=-L${netcdf}/lib

# Compiler flags needed to build this version of NetCDF-Fortran with gfortran from GCC 10:
export FCFLAGS="-w -fallow-argument-mismatch -O2"
export FFLAGS="-w -fallow-argument-mismatch -O2"

./configure --prefix=${prefix} \
            --disable-fortran-type-check | tee netcdf-f-configure.log
make check | tee netcdf-f-make-check.log
make install | tee netcdf-f-make-install.log

popd
popd

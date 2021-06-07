#!/bin/bash

set -u
set -e
set -x

# Install CDI-PIO
netcdf=~/local/cdipio/
ppm=~/local/cdipio/
yaxt=~/local/cdipio/

prefix=$HOME/local/cdipio/
build_dir=~/local/build/cdipio/cdi
url=https://gitlab.dkrz.de/mpim-sw/libcdi.git

mkdir -p ${build_dir}
pushd ${build_dir}

git clone --depth=1 -b cdipio-dev-snapshot-20210223 ${url} || true

pushd libcdi
# Increase the limit of attributes for each variable
sed -i 's/256/4096/g' src/cdi_limits.h

export CC=mpicc
export FC=mpifort
export CFLAGS="-std=gnu99 -O3"
export PPM_CORE_CFLAGS=-I${ppm}/include
export PPM_CORE_LIBS="-L${ppm}/lib -lscalesppmcore -lscalesppm"
export YAXT_CFLAGS="-I${yaxt}/include"
export YAXT_C_LIBS="-L${yaxt}/lib -lyaxt_c"
export YAXT_LIBS="-L${yaxt}/lib -lyaxt"
export ./configure --prefix=${prefix} \
            --with-netcdf=${netcdf} \
            --disable-grib \
            --disable-cgribex \
            --disable-service \
            --disable-extra \
            --disable-ieg \
            --enable-mpi \
            --enable-shared \
            --enable-static

make
make install

popd
popd

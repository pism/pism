#!/bin/bash

set -u
set -e
set -x

# Install CDI-PIO

ppm=~/local/PPM
yaxt=~/local/YAXT
mpi=~/local/MPI

prefix=$HOME/local/CDI
build_dir=~/local/build/CDI
url=https://gitlab.dkrz.de/mpim-sw/libcdi.git

mkdir -p ${build_dir}
pushd ${build_dir}

git clone -b cdipio-dev-snapshot-20210223 '${url}'

pushd cdi-${version}
# Increase the limit of attributes for each variable
sed -i 's/256/4096/g' src/cdi_limits.h

./configure --with-netcdf=${netcdf} --disable-grib --disable-cgribex \
--disable-service --disable-extra --disable-ieg --prefix=${prefix} \
--enable-mpi --enable-shared --enable-static --with-mpi-root=${mpi} \
PPM_CORE_CFLAGS=-I${ppm}/include PPM_CORE_LIBS="-L${ppm}/lib -lscalesppmcore -lscalesppm" \
YAXT_CFLAGS="-I${yaxt}/include" YAXT_C_LIBS="-L${yaxt}/lib -lyaxt_c" \
YAXT_LIBS="-L${yaxt}/lib -lyaxt" CFLAGS="-std=gnu99 -O3"

make
make install

popd
popd
#!/bin/bash

set -u
set -e
set -x

rm -rf yac

# manual-begin
version=3.3.0
prefix=$HOME/local/yac

git clone -b release-${version} \
    https://gitlab.dkrz.de/dkrz-sw/yac.git

cd yac

autoreconf -i

export CC=mpicc FC=mpifort CFLAGS="-O3 -g -march=native"

./configure --prefix=${prefix} \
            --with-yaxt-root=${prefix} \
            --disable-netcdf \
            --disable-examples \
            --disable-tools \
            --disable-deprecated \
            --with-pic

make all && make install
# manual-end

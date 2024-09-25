#!/bin/bash

set -u
set -e
set -x

rm -rf yaxt
rm -rf yac

# manual-begin
prefix=$HOME/local/yac

yaxt_version=0.11.2
git clone -b release-${yaxt_version} \
    https://gitlab.dkrz.de/dkrz-sw/yaxt.git

cd yaxt

autoreconf -i

export CC=mpicc FC=mpifort CFLAGS="-O3 -g -march=native"

./configure --prefix=${prefix} \
            --with-pic

make all && make install

cd -

yac_version=3.3.0
git clone -b release-${yac_version} \
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
